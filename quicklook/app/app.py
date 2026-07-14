import io
import json
import math
import os
import queue
import re
import sqlite3
import sys
import threading
import time
import traceback
from collections import OrderedDict
from pathlib import Path
from threading import Thread, Event

import matplotlib.pyplot as pl
from flask import Flask, render_template, request, jsonify, url_for
from flask_sock import Sock
from loguru import logger
from werkzeug.middleware.proxy_fix import ProxyFix
from quicklook.tql import TessQuickLook
from quicklook.pipelines import ALL_TESS_PIPELINES, HLSP_PIPELINES
from quicklook.cli.ql import sanitize_target_name
from quicklook.exceptions import InvalidInputError, QuickLookError
from quicklook.utils import get_available_pipelines, get_available_sectors

# Directories
BASE_DIR = Path(__file__).parent.resolve()
OUTPUT_DIR = BASE_DIR / "static" / "outputs"
LOG_DIR = OUTPUT_DIR / "logs"
DB_PATH = OUTPUT_DIR / "jobs.db"
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)

app = Flask(
    __name__, static_folder=str(BASE_DIR / "static"), template_folder=str(BASE_DIR / "templates")
)
# Honor only the mount prefix supplied by the trusted muscat-db gateway.  Host,
# scheme, port, and client-address forwarding remain disabled deliberately.
app.wsgi_app = ProxyFix(app.wsgi_app, x_for=0, x_proto=0, x_host=0, x_port=0, x_prefix=1)
sock = Sock(app)

# Original stdout/stderr for fallback
_real_stdout = sys.stdout
_real_stderr = sys.stderr


class _ThreadLocalStream(io.TextIOBase):
    """A sys.stdout/stderr replacement that routes writes per-thread.

    Each job thread registers its log file via set_stream(). Threads without
    a registered stream (e.g. the main Flask thread) fall through to the
    original stream. This avoids the race where concurrent batch jobs stomp
    on each other's global sys.stdout redirect.
    """

    def __init__(self, fallback):
        self._fallback = fallback
        self._local = threading.local()

    def set_stream(self, stream):
        self._local.stream = stream

    def set_fallback(self, stream):
        self._fallback = stream

    def clear_stream(self):
        self._local.stream = None

    @property
    def _current(self):
        return getattr(self._local, "stream", None) or self._fallback

    def write(self, s):
        return self._current.write(s)

    def flush(self):
        return self._current.flush()

    def fileno(self):
        return self._current.fileno()

    @property
    def encoding(self):
        return getattr(self._current, "encoding", "utf-8")


_tls_stdout = _ThreadLocalStream(_real_stdout)
_tls_stderr = _ThreadLocalStream(_real_stderr)


# ---------------------------------------------------------------------------
# SQLite job history (schema + CRUD live in quicklook.app.jobs_db)
# ---------------------------------------------------------------------------
from quicklook.app.jobs_db import JobsDB  # noqa: E402

_jobs_db = JobsDB(DB_PATH)


def _save_job_history(
    name, status, error="", params=None, submitted_at=None, finished_at=None, step_times=None
):
    _jobs_db.save(
        name=name,
        status=status,
        error=error,
        params=params,
        submitted_at=submitted_at,
        finished_at=finished_at,
        step_times=step_times,
    )


_avg_step_cache = {"data": {}, "ts": 0}
_avg_step_lock = threading.Lock()


def _get_avg_step_times():
    """Return average duration (seconds) per pipeline step from recent successful jobs."""
    now = time.time()
    with _avg_step_lock:
        if now - _avg_step_cache["ts"] < 60:
            return dict(_avg_step_cache["data"])
    step_dicts = _jobs_db.recent_step_times(limit=20)
    if not step_dicts:
        return {}
    totals: dict[str, float] = {}
    counts: dict[str, int] = {}
    for st in step_dicts:
        for step, dur in st.items():
            totals[step] = totals.get(step, 0) + dur
            counts[step] = counts.get(step, 0) + 1
    result = {step: totals[step] / counts[step] for step in totals}
    with _avg_step_lock:
        _avg_step_cache["data"] = result
        _avg_step_cache["ts"] = now
    return result


# ---------------------------------------------------------------------------
# Job queue (in-memory for active jobs)
# ---------------------------------------------------------------------------
jobs: OrderedDict = OrderedDict()
jobs_lock = threading.Lock()

# Serial execution queue — lightkurve, astropy, and matplotlib are not
# thread-safe, so pipeline jobs must run one at a time.
_job_queue: queue.Queue = queue.Queue()


def _job_worker():
    """Worker that pulls jobs off the queue and runs them one at a time."""
    while True:
        name, cancel_event, params = _job_queue.get()
        try:
            run_quicklook_background(name=name, cancel_event=cancel_event, **params)
        except Exception:
            pass  # run_quicklook_background handles its own errors
        finally:
            _job_queue.task_done()


_worker_thread = Thread(target=_job_worker, daemon=True)
_worker_thread.start()

# Label for the TGLC local ePSF-extraction step. Kept as a named constant
# because the websocket handler matches on it for fine-grained sub-progress.
EPSF_STEP_LABEL = "Extracting ePSF light curve"

# Pipeline step definitions for progress tracking.
PIPELINE_STEPS = [
    (r"Generating quicklook", "Initializing"),
    (r"Catalog names:|TIC \d+", "Resolving target"),
    (r"Available sectors:|All available lightcurves", "Searching lightcurves"),
    (r"Downloading|search_lightcurve|Using .+ TPF", "Downloading data"),
    (r"ePSF fitting", EPSF_STEP_LABEL),
    (r"Plotting raw light curve|raw lc", "Raw light curve"),
    (r"flatten|biweight|cosine|Flattening", "Flattening light curve"),
    (r"Lomb-Scargle|GLS|Generalized Lomb", "Lomb-Scargle periodogram"),
    (r"Transit Least Squares|TLS|Searching.*periods", "Transit search (TLS)"),
    (r"Appending TLS results", "Processing TLS results"),
    (r"Plotting TLS|Plotting flattened|Plotting phase", "Generating plots"),
    (r"Plotting TPF|tesscut|DSS|archival", "Archival image overlay"),
    (r"Plotting odd-even|secondary eclipse|summary panel", "Final panels"),
    (r"Saved:.*\.png|Runtime:", "Done"),
]


def _detect_step(log_text):
    """Return (current_step_index, total_steps) by scanning log text."""
    total = len(PIPELINE_STEPS)
    last_matched = -1
    for i, (pattern, _) in enumerate(PIPELINE_STEPS):
        if re.search(pattern, log_text, re.IGNORECASE):
            last_matched = i
    return last_matched, total


# ---------------------------------------------------------------------------
# Background job runner
# ---------------------------------------------------------------------------
def run_quicklook_background(name, cancel_event, **kwargs):
    # Atomically: if the job was cancelled (or deleted) while sitting in
    # the queue, finalize its history row and bail. Doing this under the
    # lock closes the race where a cancel arrives between an unlocked
    # `is_set()` check and the `status="running"` write that would
    # otherwise silently overwrite the user's cancellation.
    with jobs_lock:
        info = jobs.get(name)
        if info is None or cancel_event.is_set() or info.get("status") == "cancelled":
            if info is not None:
                info["status"] = "cancelled"
                info.setdefault("finished_at", time.time())
                try:
                    _save_job_history(
                        name=name,
                        status="cancelled",
                        error="",
                        params=info.get("params"),
                        submitted_at=info.get("submitted_at"),
                        finished_at=info["finished_at"],
                        step_times=info.get("step_times", {}),
                    )
                except sqlite3.DatabaseError as e:
                    logger.warning(f"Could not persist cancelled job {name}: {e}")
            return
        info["status"] = "running"

    # Allow overriding the target name passed to the pipeline (e.g. for each-sector jobs
    # where the job name is "TOI1234_s34" but the pipeline target is "TOI1234")
    target_name = kwargs.pop("target_name", name)
    log_file = LOG_DIR / f"{name}.log"
    jobs[name]["log_file"] = str(log_file)
    sink_id = None
    log_fh = None  # explicit file handle — closed in finally, not via `with`

    try:
        # Parse period limits
        period_min = kwargs.get("period_min")
        period_max = kwargs.get("period_max")
        period_limits = None
        if period_min is not None and period_max is not None:
            period_limits = (float(period_min), float(period_max))
        elif period_min is not None:
            period_limits = (float(period_min), None)
        elif period_max is not None:
            period_limits = (None, float(period_max))

        # Parse sigma clip values
        def parse_sigma(val, default_lo=10, default_hi=5):
            if val is None:
                return None
            val = str(val).strip()
            if val == "" or val.lower() == "none":
                return None
            parts = val.split()
            if len(parts) >= 2:
                return (float(parts[0]), float(parts[1]))
            elif len(parts) == 1:
                return (float(parts[0]), default_hi)
            return (default_lo, default_hi)

        sigma_clip_raw = parse_sigma(kwargs.get("sigma_clip_raw", "10 5"))
        sigma_clip_flat = parse_sigma(kwargs.get("sigma_clip_flat"))

        # Parse custom ephemeris
        custom_ephem = kwargs.get("custom_ephem")
        if custom_ephem and custom_ephem.strip():
            parts = custom_ephem.strip().split()
            if len(parts) == 6:
                custom_ephem = [float(x) for x in parts]
            else:
                custom_ephem = None
        else:
            custom_ephem = None

        # Determine outdir and suffix
        outdir = kwargs.get("outdir")
        outdir = str(OUTPUT_DIR) if not outdir else outdir
        suffix = kwargs.get("suffix")
        suffix = None if not suffix else suffix

        # Open the log file and set up streams — file stays open until the
        # outermost finally block so no concurrent write hits a closed FD.
        # Keep stdout/stderr in append mode too. Loguru uses a separate
        # append-mode descriptor below; a normal ``w`` descriptor would keep
        # its own offset and overwrite Loguru messages written in between.
        log_fh = open(log_file, "a", buffering=1, encoding="utf-8")
        log_fh.seek(0)
        log_fh.truncate()
        _tls_stdout.set_stream(log_fh)
        _tls_stderr.set_stream(log_fh)

        # Add loguru sink (uses its own FD to the same path in append mode)
        sink_id = logger.add(
            str(log_file),
            format="{message}",
            level="DEBUG",
            filter="quicklook",
            mode="a",
            encoding="utf-8",
        )

        if cancel_event.is_set():
            raise InterruptedError("Job cancelled before start")
        tql = TessQuickLook(
            target_name=target_name,
            sector=int(kwargs.get("sector", -1)),
            pipeline=kwargs.get("pipeline", "spoc"),
            flux_type=kwargs.get("fluxtype", "pdcsap"),
            # Coerce "" (the form's "auto" option) to None so the pipeline
            # treats it as "let MAST decide" instead of filtering for an
            # empty-string exptime and accidentally excluding every row.
            exptime=(kwargs.get("exptime") or None),
            quality_bitmask=kwargs.get("quality_bitmask", "default"),
            flatten_method=kwargs.get("flatten_method", "biweight"),
            gp_kernel=kwargs.get("gp_kernel", "periodic_auto"),
            window_length=kwargs.get("window_length"),
            pg_method=kwargs.get("pg_method", "gls"),
            edge_cutoff=kwargs.get("edge_cutoff", 0.1),
            sigma_clip_raw=sigma_clip_raw,
            sigma_clip_flat=sigma_clip_flat,
            Porb_limits=period_limits,
            mask_ephem=kwargs.get("mask_ephem", False),
            use_star_priors=kwargs.get("use_priors", False),
            phase_xlim=kwargs.get("phase_xlim"),
            custom_ephem=custom_ephem,
            archival_survey=kwargs.get("survey", "dss1"),
            show_simbad=kwargs.get("show_simbad", False),
            savefig=kwargs.get("save", False),
            savetls=kwargs.get("save", False),
            verbose=kwargs.get("verbose", True),
            overwrite=kwargs.get("overwrite", False),
            outdir=outdir,
            suffix=suffix,
            # Pass the Event's bound is_set so the long TGLC ePSF loop can
            # actually interrupt when the GUI Cancel button is clicked,
            # instead of running to completion and only honouring the
            # cancel at the next phase boundary.
            cancel_check=cancel_event.is_set,
        )
        if cancel_event.is_set():
            raise InterruptedError("Job cancelled")
        tql.plot_tql(return_fig_and_paths=True)
        if cancel_event.is_set():
            jobs[name]["status"] = "cancelled"
        else:
            jobs[name]["status"] = "done"

    except InterruptedError:
        jobs[name]["status"] = "cancelled"
    except QuickLookError as e:
        jobs[name]["status"] = "error"
        jobs[name]["error"] = str(e)
        logger.warning(f"Job failed: {e}")
    except Exception as e:
        tb = traceback.format_exc()
        jobs[name]["status"] = "error"
        # Keep the short message for the GUI badge, but include the traceback
        # so the per-target .log captures the offending file/line.
        jobs[name]["error"] = str(e)
        logger.error(f"Job '{name}' failed unexpectedly: {e}\n{tb}")
    finally:
        # pyplot holds every figure in a global registry, so dropping the local
        # reference is not enough. The CLI exits after each target, but this
        # worker outlives the job and would otherwise accumulate a figure per
        # target across a batch. Jobs run serially on one thread (see
        # _submit_job), so no other thread owns a figure here.
        pl.close("all")
        # Tear down in safe order: streams first, then loguru sink, then file
        _tls_stdout.clear_stream()
        _tls_stderr.clear_stream()
        if sink_id is not None:
            logger.remove(sink_id)
        if log_fh is not None:
            log_fh.close()
        # Persist to SQLite
        jobs[name]["finished_at"] = time.time()
        _save_job_history(
            name=name,
            status=jobs[name]["status"],
            error=jobs[name].get("error", ""),
            params=jobs[name].get("params"),
            submitted_at=jobs[name].get("submitted_at"),
            finished_at=jobs[name]["finished_at"],
            step_times=jobs[name].get("step_times", {}),
        )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _get_recent_results(limit=12):
    """Return the most recent output PNGs as dicts with path and filename."""
    images = sorted(OUTPUT_DIR.glob("*.png"), key=os.path.getmtime, reverse=True)
    return [
        {"path": url_for("static", filename=f"outputs/{img.name}"), "name": img.name}
        for img in images[:limit]
    ]


def _is_truthy(val):
    """Check if a form / JSON field is truthy (handles checkboxes *and* booleans)."""
    if val is None:
        return False
    if isinstance(val, bool):
        return val
    if isinstance(val, str):
        return val.lower() in ("true", "on", "1", "yes")
    return bool(val)


def _parse_params(form):
    """Parse form-data or JSON dict into pipeline kwargs."""

    def safe_int(val, default=-1):
        try:
            return int(val)
        except (ValueError, TypeError):
            return default

    def safe_float(val, default=0.1):
        try:
            return float(val)
        except (ValueError, TypeError):
            return default

    return {
        "sector": safe_int(form.get("sector", -1), -1),
        "pipeline": form.get("pipeline", "spoc"),
        "fluxtype": form.get("fluxtype", "pdcsap"),
        "exptime": safe_int(form.get("exptime") or None, None),
        "quality_bitmask": form.get("quality_bitmask", "default"),
        "flatten_method": form.get("flatten_method", "biweight"),
        "gp_kernel": form.get("gp_kernel", "periodic_auto"),
        "window_length": safe_float(form.get("window_length") or None, None),
        "pg_method": form.get("pg_method", "gls"),
        "edge_cutoff": safe_float(form.get("edge_cutoff", 0.1), 0.1),
        "sigma_clip_raw": form.get("sigma_clip_raw", "10 5"),
        "sigma_clip_flat": form.get("sigma_clip_flat"),
        "period_min": safe_float(form.get("period_min") or None, None),
        "period_max": safe_float(form.get("period_max") or None, None),
        "phase_xlim": safe_float(form.get("phase_xlim") or None, None),
        "mask_ephem": _is_truthy(form.get("mask_ephem")),
        "use_priors": _is_truthy(form.get("use_priors")),
        "custom_ephem": form.get("custom_ephem"),
        "survey": form.get("survey", "dss1"),
        "show_simbad": _is_truthy(form.get("show_simbad")),
        "outdir": form.get("outdir"),
        "suffix": form.get("suffix") if form.get("suffix") else "",
        "save": _is_truthy(form.get("save")),
        "verbose": _is_truthy(form.get("verbose")),
        "overwrite": _is_truthy(form.get("overwrite")),
    }


def _job_dedup_signature(params):
    """Key fields that make a submission semantically distinct.

    Two non-terminal jobs with the same target name AND the same signature
    are considered true duplicates; the same name with a different signature
    (e.g. SPOC vs TGLC, different sector, different flux/photometry) is a
    different job and gets auto-disambiguated rather than rejected.
    """
    return (
        params.get("pipeline", "spoc"),
        params.get("sector", -1),
        params.get("fluxtype", "pdcsap"),
        params.get("exptime", None),
    )


def _submit_job(name, params):
    """Enqueue a job for serial execution. Returns the job name."""
    cancel_event = Event()
    with jobs_lock:
        jobs[name] = {
            "status": "queued",
            "log_file": "",
            "error": "",
            "cancel_event": cancel_event,
            "submitted_at": time.time(),
            "params": params,
            "step_times": {},
        }
        jobs.move_to_end(name)

    _job_queue.put((name, cancel_event, params))
    return name


# ---------------------------------------------------------------------------
# Routes
# ---------------------------------------------------------------------------
@app.route("/")
def index():
    recent_results = _get_recent_results()
    return render_template(
        "index.html",
        recent_results=recent_results,
        pipelines=ALL_TESS_PIPELINES,
    )


@app.errorhandler(InvalidInputError)
def handle_invalid_input(e):
    """Turn a rejected target name into a 400 rather than a 500 + traceback."""
    return jsonify({"ok": False, "reason": str(e)}), 400


@app.route("/submit", methods=["POST"])
def submit_job():
    """AJAX endpoint for single-job submission."""
    form = request.form if request.form else (request.json or {})
    name = (form.get("name") or "").strip()
    if not name:
        return jsonify({"ok": False, "reason": "Target name is required"}), 400

    name = sanitize_target_name(name)
    params = _parse_params(form)

    new_sig = _job_dedup_signature(params)
    with jobs_lock:
        existing = jobs.get(name)
        if existing and existing["status"] not in ("done", "error", "cancelled"):
            if _job_dedup_signature(existing.get("params", {})) == new_sig:
                # True duplicate: same name AND same key settings.
                return (
                    jsonify(
                        {
                            "ok": False,
                            "reason": (f"Job '{name}' with these settings is already running"),
                        }
                    ),
                    409,
                )
            # Same target, different settings (e.g. SPOC vs TGLC). Auto-pick
            # an unused name so both can coexist in the queue.
            base = f"{name}_{params.get('pipeline', 'spoc')}"
            candidate = base
            n = 2
            while True:
                other = jobs.get(candidate)
                if other is None or other["status"] in (
                    "done",
                    "error",
                    "cancelled",
                ):
                    break
                candidate = f"{base}_{n}"
                n += 1
            name = candidate

    _submit_job(name, params)
    return jsonify({"ok": True, "name": name})


@app.route("/batch-submit", methods=["POST"])
def batch_submit():
    """Submit multiple targets sharing the same pipeline parameters."""
    data = request.json or {}
    targets = data.get("targets", [])
    params_template = _parse_params(data.get("params", {}))

    if not targets:
        return jsonify({"ok": False, "reason": "No targets provided"}), 400

    submitted = []
    skipped = []
    for raw_name in targets:
        raw_name = raw_name.strip()
        if not raw_name:
            continue
        try:
            name = sanitize_target_name(raw_name)
        except InvalidInputError:
            # One bad name should not abort the rest of the batch.
            skipped.append(raw_name)
            continue
        with jobs_lock:
            if name in jobs and jobs[name]["status"] not in ("done", "error", "cancelled"):
                skipped.append(name)
                continue
        _submit_job(name, dict(params_template))
        submitted.append(name)

    return jsonify({"ok": True, "submitted": submitted, "skipped": skipped})


@app.route("/available-sectors")
def available_sectors():
    """Return the list of available sectors for a target + pipeline."""
    name = request.args.get("name", "").strip()
    pipeline = request.args.get("pipeline", "spoc").strip()
    if not name:
        return jsonify({"ok": False, "reason": "Target name is required"}), 400
    name = sanitize_target_name(name)
    try:
        sectors = get_available_sectors(name, pipeline=pipeline)
    except Exception as e:
        return jsonify({"ok": False, "reason": str(e)}), 500
    if not sectors:
        return jsonify({"ok": False, "reason": f"No sectors found for {name} with {pipeline}"}), 404
    return jsonify({"ok": True, "sectors": sectors})


@app.route("/each-sector-submit", methods=["POST"])
def each_sector_submit():
    """Submit one job per sector for a single target (--each-sector equivalent)."""
    data = request.json or {}
    name = (data.get("name") or "").strip()
    sectors = data.get("sectors", [])
    params_template = _parse_params(data.get("params", {}))

    if not name:
        return jsonify({"ok": False, "reason": "Target name is required"}), 400
    if not sectors:
        return jsonify({"ok": False, "reason": "No sectors provided"}), 400

    name = sanitize_target_name(name)
    submitted = []
    skipped = []
    for sector in sectors:
        job_name = f"{name}_s{int(sector)}"
        params = dict(params_template)
        params["sector"] = int(sector)
        params["target_name"] = name  # pass original target name to the pipeline
        with jobs_lock:
            if job_name in jobs and jobs[job_name]["status"] not in ("done", "error", "cancelled"):
                skipped.append(job_name)
                continue
        _submit_job(job_name, params)
        submitted.append(job_name)

    return jsonify({"ok": True, "submitted": submitted, "skipped": skipped})


@app.route("/available-pipelines")
def available_pipelines():
    """Return the list of available pipelines for a target (any sector)."""
    name = request.args.get("name", "").strip()
    if not name:
        return jsonify({"ok": False, "reason": "Target name is required"}), 400
    name = sanitize_target_name(name)
    try:
        pipelines = get_available_pipelines(name)
    except Exception as e:
        return jsonify({"ok": False, "reason": str(e)}), 500
    if not pipelines:
        return jsonify({"ok": False, "reason": f"No pipelines found for {name}"}), 404
    return jsonify({"ok": True, "pipelines": pipelines})


@app.route("/each-pipeline-submit", methods=["POST"])
def each_pipeline_submit():
    """Submit one job per pipeline for a single target (--each-pipeline equivalent).

    Sector is forced to -1 (latest available) per pipeline so each pipeline
    runs on whatever sector(s) it actually has. Useful for "fan out by
    pipeline" exploration on a single target.
    """
    data = request.json or {}
    name = (data.get("name") or "").strip()
    pipelines = data.get("pipelines", [])
    params_template = _parse_params(data.get("params", {}))

    if not name:
        return jsonify({"ok": False, "reason": "Target name is required"}), 400
    if not pipelines:
        return jsonify({"ok": False, "reason": "No pipelines provided"}), 400

    name = sanitize_target_name(name)
    submitted = []
    skipped = []
    for pipeline in pipelines:
        pipeline = str(pipeline).lower().strip()
        if not pipeline:
            continue
        job_name = f"{name}_{pipeline}"
        params = dict(params_template)
        params["pipeline"] = pipeline
        params["sector"] = -1  # latest available for this pipeline
        params["target_name"] = name
        with jobs_lock:
            if job_name in jobs and jobs[job_name]["status"] not in ("done", "error", "cancelled"):
                skipped.append(job_name)
                continue
        _submit_job(job_name, params)
        submitted.append(job_name)

    return jsonify({"ok": True, "submitted": submitted, "skipped": skipped})


@app.route("/jobs-json")
def jobs_json():
    """Return status of all jobs (for queue display)."""
    result = []
    with jobs_lock:
        for name, info in jobs.items():
            result.append(
                {
                    "name": name,
                    "status": info["status"],
                    "error": info.get("error", ""),
                    "submitted_at": info.get("submitted_at"),
                }
            )
    return jsonify(result)


@app.route("/cancel/<target>", methods=["POST"])
def cancel(target):
    """Cancel a running job."""
    info = jobs.get(target)
    if info and info["status"] in ("running", "queued"):
        info["cancel_event"].set()
        info["status"] = "cancelled"
        return jsonify({"ok": True})
    return jsonify({"ok": False, "reason": "not running"}), 400


@app.route("/cancel-all", methods=["POST"])
def cancel_all():
    """Cancel all pending (queued or running) jobs."""
    cancelled = []
    for name, info in jobs.items():
        if info["status"] in ("running", "queued"):
            info["cancel_event"].set()
            info["status"] = "cancelled"
            cancelled.append(name)
    return jsonify({"ok": True, "cancelled": cancelled})


@app.route("/delete-job/<target>", methods=["POST"])
def delete_job(target):
    """Remove a job's record from the in-memory queue, SQLite history,
    and the per-target log file.

    Running jobs are refused — they should be cancelled first via
    ``/cancel/<target>``. Queued jobs are accepted: the cancel_event
    is set so the worker skips them once they're picked up.
    """
    with jobs_lock:
        info = jobs.get(target)
        if info and info.get("status") == "running":
            return (
                jsonify(
                    {
                        "ok": False,
                        "reason": "Job is still running; cancel it first.",
                    }
                ),
                400,
            )
        # If queued, mark the cancel_event so the worker doesn't try to
        # touch jobs[target] after we've popped it.
        if info and info.get("status") == "queued":
            evt = info.get("cancel_event")
            if evt is not None:
                evt.set()
        jobs.pop(target, None)

    # Drop the SQLite history row (best-effort).
    try:
        _jobs_db.delete(target)
    except sqlite3.DatabaseError as e:
        logger.warning(f"Failed to delete job_history row for {target}: {e}")

    # Drop the per-target log file (best-effort).
    log_path = LOG_DIR / f"{target}.log"
    try:
        if log_path.exists():
            log_path.unlink()
    except OSError as e:
        logger.warning(f"Failed to delete log file {log_path}: {e}")

    return jsonify({"ok": True})


@app.route("/delete-output", methods=["POST"])
def delete_output():
    """Delete an output PNG and its associated H5 file."""
    filename = request.json.get("filename", "")
    if not filename or ".." in filename or "/" in filename:
        return jsonify({"ok": False, "reason": "invalid filename"}), 400

    png_path = OUTPUT_DIR / filename
    deleted = []
    if png_path.exists() and png_path.suffix == ".png":
        png_path.unlink()
        deleted.append(filename)
        h5_name = png_path.stem + "_tls.h5"
        h5_path = OUTPUT_DIR / h5_name
        if h5_path.exists():
            h5_path.unlink()
            deleted.append(h5_name)

    if deleted:
        return jsonify({"ok": True, "deleted": deleted})
    return jsonify({"ok": False, "reason": "file not found"}), 404


@app.route("/bulk-delete", methods=["POST"])
def bulk_delete():
    """Delete multiple output files at once."""
    filenames = request.json.get("filenames", [])
    if not filenames:
        return jsonify({"ok": False, "reason": "No filenames provided"}), 400

    all_deleted = []
    for filename in filenames:
        if not filename or ".." in filename or "/" in filename:
            continue
        png_path = OUTPUT_DIR / filename
        if png_path.exists() and png_path.suffix == ".png":
            png_path.unlink()
            all_deleted.append(filename)
            h5_name = png_path.stem + "_tls.h5"
            h5_path = OUTPUT_DIR / h5_name
            if h5_path.exists():
                h5_path.unlink()
                all_deleted.append(h5_name)

    return jsonify({"ok": True, "deleted": all_deleted})


@app.route("/results-json/<target>")
def results_json(target):
    """Return paths to output PNG and H5 files for a completed job."""
    # Output files use names without hyphens (e.g. TOI4152 not TOI-4152)
    search_name = target.replace("-", "")
    images = sorted(OUTPUT_DIR.glob(f"*{search_name}*.png"), key=os.path.getmtime, reverse=True)
    h5s = sorted(OUTPUT_DIR.glob(f"*{search_name}*_tls.h5"), key=os.path.getmtime, reverse=True)
    return jsonify(
        {
            "image": url_for("static", filename=f"outputs/{images[0].name}") if images else None,
            "h5": url_for("static", filename=f"outputs/{h5s[0].name}") if h5s else None,
        }
    )


@app.route("/job-log/<target>")
def job_log(target):
    """Return the persisted per-target output log.

    The live WebSocket (``/ws/<target>``) only streams jobs still present in the
    in-memory ``jobs`` queue; once a job finishes and is evicted — or after a
    server restart — its ``.log`` file on disk is the only remaining record.
    This endpoint reads that file directly so the GUI can show a finished job's
    output on demand, without a live connection.
    """
    # Job names never contain path separators or ``..`` (see
    # sanitize_target_name), so reject anything that could escape LOG_DIR.
    if not target or "/" in target or "\\" in target or ".." in target or target.startswith("."):
        return jsonify({"ok": False, "reason": "invalid target"}), 400

    log_path = LOG_DIR / f"{target}.log"
    # Defense in depth: confirm the resolved path really sits inside LOG_DIR.
    try:
        resolved = log_path.resolve()
        resolved.relative_to(Path(LOG_DIR).resolve())
    except (ValueError, OSError):
        return jsonify({"ok": False, "reason": "invalid target"}), 400

    if not resolved.exists():
        return jsonify({"ok": False, "reason": "no log for this job"}), 404

    try:
        content = resolved.read_text(encoding="utf-8", errors="replace")
    except OSError as e:
        logger.warning(f"Failed to read log file {resolved}: {e}")
        return jsonify({"ok": False, "reason": "could not read log"}), 500

    info = jobs.get(target)
    status = info.get("status") if info else None
    return jsonify({"ok": True, "log": content, "status": status})


@app.route("/avg-step-times")
def avg_step_times():
    """Return average pipeline step durations for ETA estimation."""
    return jsonify(_get_avg_step_times())


@sock.route("/ws/<target>")
def ws_log(ws, target):
    """WebSocket endpoint for live log streaming, progress, and ETA."""
    info = jobs.get(target)
    if not info:
        ws.send(json.dumps({"type": "finish", "status": "unknown"}))
        return

    last_size = 0
    last_step_idx = -1
    step_start_time = time.time()

    while True:
        # Check for incoming messages (cancel command)
        try:
            msg = ws.receive(timeout=0)
            if msg:
                data = json.loads(msg)
                if data.get("action") == "cancel":
                    evt = info.get("cancel_event")
                    if evt:
                        evt.set()
                    info["status"] = "cancelled"
        except (json.JSONDecodeError, TypeError, OSError):
            pass

        status = info.get("status", "unknown")

        # Read new log content (re-read path each iteration in case job thread hasn't set it yet)
        new_lines = ""
        full_log = ""
        log_path = Path(info["log_file"]) if info.get("log_file") else None
        if log_path and log_path.exists():
            content = log_path.read_text()
            full_log = content
            if len(content) > last_size:
                new_lines = content[last_size:]
                last_size = len(content)

        # Detect progress & compute ETA
        eta_seconds = None
        if full_log:
            step_idx, step_total = _detect_step(full_log)
            step_label = PIPELINE_STEPS[step_idx][1] if step_idx >= 0 else "Starting"
            pct = int(((step_idx + 1) / step_total) * 100) if step_idx >= 0 else 0

            # Fine-grained sub-progress within the long TGLC ePSF step:
            # interpolate the bar from the "ePSF fitting: N/M (X%)" log lines.
            if step_idx >= 0 and PIPELINE_STEPS[step_idx][1] == EPSF_STEP_LABEL:
                epsf_pcts = re.findall(r"ePSF fitting: \d+/\d+ \((\d+)%\)", full_log)
                if epsf_pcts:
                    frac = int(epsf_pcts[-1]) / 100.0
                    pct = int(((step_idx + frac) / step_total) * 100)
                    step_label = f"{EPSF_STEP_LABEL} ({epsf_pcts[-1]}%)"

            # Record step timing transitions
            if step_idx > last_step_idx:
                now = time.time()
                if last_step_idx >= 0:
                    prev_label = PIPELINE_STEPS[last_step_idx][1]
                    info.setdefault("step_times", {})[prev_label] = now - step_start_time
                step_start_time = now
                last_step_idx = step_idx

            # Estimate remaining time from historical averages
            if step_idx >= 0 and status == "running":
                avg_times = _get_avg_step_times()
                remaining = 0.0
                for i in range(step_idx + 1, step_total):
                    remaining += avg_times.get(PIPELINE_STEPS[i][1], 15)
                elapsed_in_step = time.time() - step_start_time
                avg_current = avg_times.get(PIPELINE_STEPS[step_idx][1], 15)
                remaining += max(0, avg_current - elapsed_in_step)
                eta_seconds = max(0, int(remaining))
        else:
            step_label = "Queued" if status == "queued" else "Starting"
            pct = 0

        # Send update
        try:
            payload = {
                "type": "update",
                "log": new_lines,
                "status": status,
                "progress": pct,
                "step": step_label,
            }
            if eta_seconds is not None:
                payload["eta"] = eta_seconds
            ws.send(json.dumps(payload))
        except (OSError, ConnectionError, BrokenPipeError):
            break  # Client disconnected

        if status in ("done", "error", "cancelled"):
            # Record final step duration
            if last_step_idx >= 0:
                final_label = PIPELINE_STEPS[last_step_idx][1]
                info.setdefault("step_times", {})[final_label] = time.time() - step_start_time

            try:
                ws.send(
                    json.dumps(
                        {
                            "type": "finish",
                            "status": status,
                            "error": info.get("error", ""),
                        }
                    )
                )
            except (OSError, ConnectionError, BrokenPipeError):
                pass
            break

        time.sleep(1.5)


GALLERY_PER_PAGE_CHOICES = [12, 24, 48, 96, 192]
GALLERY_PER_PAGE_DEFAULT = 24


@app.route("/gallery")
def gallery():
    search_name = request.args.get("search", "")
    sort_by = request.args.get("sort", "date_desc")
    page = max(1, int(request.args.get("page", 1)))
    try:
        per_page = int(request.args.get("per_page", GALLERY_PER_PAGE_DEFAULT))
    except (TypeError, ValueError):
        per_page = GALLERY_PER_PAGE_DEFAULT
    if per_page not in GALLERY_PER_PAGE_CHOICES:
        per_page = GALLERY_PER_PAGE_DEFAULT

    images = list(OUTPUT_DIR.glob("*.png"))

    if search_name:
        images = [img for img in images if search_name.lower() in img.stem.lower()]

    # Sort
    if sort_by == "name_asc":
        images.sort(key=lambda p: p.stem.lower())
    elif sort_by == "name_desc":
        images.sort(key=lambda p: p.stem.lower(), reverse=True)
    elif sort_by == "date_asc":
        images.sort(key=os.path.getmtime)
    else:
        images.sort(key=os.path.getmtime, reverse=True)

    total = len(images)
    total_pages = max(1, math.ceil(total / per_page))
    page = min(page, total_pages)
    start = (page - 1) * per_page
    page_images = images[start : start + per_page]

    items = [
        {"path": url_for("static", filename=f"outputs/{img.name}"), "name": img.name}
        for img in page_images
    ]
    return render_template(
        "gallery.html",
        images=items,
        search_name=search_name,
        sort_by=sort_by,
        page=page,
        total_pages=total_pages,
        total=total,
        per_page=per_page,
        per_page_choices=GALLERY_PER_PAGE_CHOICES,
    )


def _parse_tls_filename(stem):
    """Best-effort parse of pipeline, flux_type, cadence from a TLS filename.

    Filenames follow `{name}_s{NN}_{lctype}_{c}c[_mask_ephem][_{suffix}]_tls`
    (see TessQuickLook.check_output_file_exists). ``lctype`` is the flux type
    when pipeline=="spoc" (pdcsap/sap), else the pipeline name itself.
    Returns a dict of fallbacks; values are None when not confidently
    recoverable.
    """
    SPOC_FLUX = {"pdcsap", "sap"}
    # Trailing "_tls" already stripped by Path.stem callers
    if stem.endswith("_tls"):
        stem = stem[: -len("_tls")]
    tokens = stem.split("_")
    pipeline = flux_type = cadence = None
    # Find the sector token (sNN) and the cadence token ((s|l)c).
    sec_idx = next((i for i, t in enumerate(tokens) if re.fullmatch(r"s\d{1,3}", t)), None)
    cad_idx = next((i for i, t in enumerate(tokens) if re.fullmatch(r"[sl]c", t)), None)
    if sec_idx is not None and cad_idx is not None and cad_idx > sec_idx + 1:
        lctype = tokens[sec_idx + 1]
        if lctype in SPOC_FLUX:
            pipeline, flux_type = "spoc", lctype
        elif lctype in HLSP_PIPELINES:
            pipeline = lctype
        cadence = "short" if tokens[cad_idx] == "sc" else "long"
    return {"pipeline": pipeline, "flux_type": flux_type, "cadence": cadence}


def _coerce_scalar(val):
    """Unwrap 0-d/1-elem numpy arrays so Jinja sees a plain Python scalar."""
    if val is None:
        return None
    try:
        import numpy as _np

        if isinstance(val, _np.ndarray):
            if val.ndim == 0:
                return val.item()
            if val.size == 1:
                return val.flat[0].item()
    except Exception:
        pass
    if isinstance(val, bytes):
        try:
            return val.decode()
        except Exception:
            return None
    return val


def _get_allowed_browse_roots():
    """Roots the /list-dirs endpoint is allowed to descend into.

    Defaults to the user's home directory. Override via the
    ``ALLOWED_BROWSE_ROOTS`` env var, colon-separated absolute paths.
    """
    env = os.environ.get("ALLOWED_BROWSE_ROOTS", "")
    if env.strip():
        roots = []
        for part in env.split(":"):
            part = part.strip()
            if not part:
                continue
            try:
                roots.append(Path(part).expanduser().resolve())
            except (OSError, ValueError):
                continue
        if roots:
            return roots
    return [Path.home().resolve()]


def _is_under_allowed_root(p, roots):
    """True if the resolved path is inside any allowed root (or equals one)."""
    try:
        p_resolved = p.resolve()
    except (OSError, ValueError):
        return False
    for root in roots:
        try:
            if p_resolved == root or p_resolved.is_relative_to(root):
                return True
        except (OSError, ValueError):
            continue
    return False


@app.route("/list-dirs")
def list_dirs():
    """List immediate subdirectories of a path. Localhost-only convenience
    for the TLS Summary directory picker — scoped to allowed roots so a
    misconfigured deployment cannot leak arbitrary filesystem layout.
    """
    roots = _get_allowed_browse_roots()
    req = request.args.get("path", "").strip()

    if not req:
        p = roots[0]
    else:
        try:
            p = Path(req).expanduser()
        except (OSError, ValueError):
            return (
                jsonify(
                    {
                        "error": "Invalid path",
                        "entries": [],
                        "path": req,
                        "parent": None,
                        "breadcrumbs": [],
                        "roots": [str(r) for r in roots],
                    }
                ),
                400,
            )

    if not _is_under_allowed_root(p, roots):
        return (
            jsonify(
                {
                    "error": "Path outside allowed roots",
                    "entries": [],
                    "path": str(p),
                    "parent": None,
                    "breadcrumbs": [],
                    "roots": [str(r) for r in roots],
                }
            ),
            403,
        )

    try:
        p = p.resolve()
    except (OSError, ValueError) as e:
        return (
            jsonify(
                {
                    "error": f"Cannot resolve path: {e}",
                    "entries": [],
                    "path": str(p),
                    "parent": None,
                    "breadcrumbs": [],
                    "roots": [str(r) for r in roots],
                }
            ),
            400,
        )

    if not p.is_dir():
        return (
            jsonify(
                {
                    "error": "Not a directory",
                    "entries": [],
                    "path": str(p),
                    "parent": str(p.parent),
                    "breadcrumbs": [],
                    "roots": [str(r) for r in roots],
                }
            ),
            400,
        )

    entries = []
    has_tls_here = False
    try:
        for entry in os.scandir(p):
            try:
                if entry.is_dir(follow_symlinks=False) and not entry.name.startswith("."):
                    sub = Path(entry.path)
                    try:
                        has_tls = any(sub.glob("*_tls.h5"))
                    except (OSError, PermissionError):
                        has_tls = False
                    entries.append(
                        {
                            "name": entry.name,
                            "path": str(sub.resolve()),
                            "has_tls": has_tls,
                        }
                    )
                elif entry.is_file(follow_symlinks=False) and entry.name.endswith("_tls.h5"):
                    has_tls_here = True
            except (OSError, PermissionError):
                continue
    except (OSError, PermissionError) as e:
        return jsonify(
            {
                "error": f"Cannot list directory: {e}",
                "entries": [],
                "path": str(p),
                "parent": str(p.parent),
                "breadcrumbs": [],
                "roots": [str(r) for r in roots],
            }
        )

    entries.sort(key=lambda e: e["name"].lower())

    # Parent only if still inside an allowed root
    parent = None
    if p != p.parent and _is_under_allowed_root(p.parent, roots):
        try:
            parent = str(p.parent.resolve())
        except (OSError, ValueError):
            parent = None

    # Breadcrumb path: walk up while still under an allowed root
    breadcrumbs = []
    current = p
    seen = set()
    while True:
        if str(current) in seen:
            break
        seen.add(str(current))
        breadcrumbs.append({"name": current.name or str(current), "path": str(current)})
        if current == current.parent or not _is_under_allowed_root(current.parent, roots):
            break
        current = current.parent
    breadcrumbs.reverse()

    return jsonify(
        {
            "path": str(p),
            "parent": parent,
            "entries": entries,
            "has_tls_here": has_tls_here,
            "breadcrumbs": breadcrumbs,
            "roots": [str(r) for r in roots],
            "error": None,
        }
    )


@app.route("/tls-summary")
def tls_summary():
    """Sortable, searchable summary table of every saved TLS .h5 result.

    Reads from ``OUTPUT_DIR`` by default; accepts ``?dir=<absolute_path>``
    to point at any other directory containing ``*_tls.h5`` files.
    """
    from quicklook import h5io

    req_dir = request.args.get("dir", "").strip()
    chosen_dir = OUTPUT_DIR
    dir_error = None
    if req_dir:
        try:
            candidate = Path(req_dir).expanduser().resolve()
        except (OSError, ValueError, RuntimeError):
            candidate = None
            dir_error = f"Invalid path: {req_dir}"
        if candidate is not None:
            if candidate.is_dir():
                chosen_dir = candidate
            else:
                dir_error = f"Path not found or not a directory: {req_dir}"

    is_default_dir = chosen_dir.resolve() == OUTPUT_DIR.resolve()
    output_root = OUTPUT_DIR.resolve()

    files = sorted(chosen_dir.glob("*_tls.h5"), key=os.path.getmtime, reverse=True)
    rows = []
    skipped = []
    for fp in files:
        try:
            data = h5io.load(str(fp))
        except Exception as e:
            skipped.append({"name": fp.name, "error": str(e)})
            continue

        stem = fp.stem  # foo_s12_pdcsap_sc_tls
        png_name = stem[: -len("_tls")] + ".png" if stem.endswith("_tls") else stem + ".png"

        # Emit /static/outputs/... URLs only when the file lives inside
        # OUTPUT_DIR; otherwise show the absolute path as plain text so we
        # don't have to add a generic file-streaming endpoint.
        try:
            fp_resolved = fp.resolve()
            external = not fp_resolved.is_relative_to(output_root)
        except (OSError, ValueError):
            external = True
            fp_resolved = fp

        if external:
            png_path = None
            h5_path = None
            disk_path = str(fp_resolved)
        else:
            png_path = (
                url_for("static", filename=f"outputs/{png_name}")
                if (OUTPUT_DIR / png_name).exists()
                else None
            )
            h5_path = url_for("static", filename=f"outputs/{fp.name}")
            disk_path = None

        def first(seq):
            try:
                return seq[0] if seq is not None and len(seq) > 0 else None
            except Exception:
                return None

        ptc = data.get("per_transit_count")
        try:
            ptc_str = ",".join(str(int(x)) for x in ptc) if ptc is not None else ""
        except Exception:
            ptc_str = str(ptc) if ptc is not None else ""

        # Filename-based fallbacks for older h5 files that pre-date the
        # "save pipeline/flux_type/exptime/cadence" fix in tql.py.
        fallback = _parse_tls_filename(stem)
        pipeline = _coerce_scalar(data.get("pipeline")) or fallback["pipeline"]
        flux_type = _coerce_scalar(data.get("flux_type")) or fallback["flux_type"]
        exptime = _coerce_scalar(data.get("exptime"))
        if exptime is None and fallback["cadence"] == "short":
            # Use a hint, not a hard claim; surface in the table so users
            # know it's inferred. 120s is the typical SPOC 2-min cadence.
            exptime_hint = 120 if pipeline == "spoc" else None
            exptime = exptime_hint

        depth = data.get("depth")
        duration_days = _coerce_scalar(data.get("duration"))
        rows.append(
            {
                "filename": fp.name,
                "stem": stem,
                "png": png_path,
                "h5": h5_path,
                "external": external,
                "disk_path": disk_path,
                "mtime": os.path.getmtime(fp),
                "toi": _coerce_scalar(data.get("toiid")),
                "tic": _coerce_scalar(data.get("ticid")),
                "gaiaid": _coerce_scalar(data.get("gaiaid")),
                "pipeline": pipeline,
                "flux_type": flux_type,
                "exptime": exptime,
                "sector": _coerce_scalar(data.get("sector")),
                "Porb": _coerce_scalar(data.get("period")),
                "duration_hr": duration_days * 24.0 if duration_days is not None else None,
                "rp_rs": _coerce_scalar(data.get("rp_rs")),
                "SDE": _coerce_scalar(data.get("SDE")),
                "snr": _coerce_scalar(data.get("snr")),
                "FAP": _coerce_scalar(data.get("FAP")),
                "odd_even": _coerce_scalar(data.get("odd_even_mismatch")),
                "n_tr": _coerce_scalar(data.get("distinct_transit_count")),
                "depth_ppt": (1 - depth) * 1e3 if depth is not None else None,
                "Prot_gls": first(data.get("Prot_gls")),
                "amp_gls": first(data.get("amp_gls")),
                "power_gls": first(data.get("power_gls")),
                "per_transit_count": ptc_str,
                "simbad_obj": _coerce_scalar(data.get("simbad_obj")),
            }
        )

    return render_template(
        "tls_summary.html",
        rows=rows,
        total=len(rows),
        skipped=skipped,
        current_dir=str(chosen_dir),
        is_default_dir=is_default_dir,
        dir_error=dir_error,
    )


@app.route("/compare")
def compare():
    """Side-by-side comparison view."""
    images = sorted(OUTPUT_DIR.glob("*.png"), key=os.path.getmtime, reverse=True)
    items = [
        {"path": url_for("static", filename=f"outputs/{img.name}"), "name": img.name}
        for img in images
    ]
    left = request.args.get("left", "")
    right = request.args.get("right", "")
    return render_template("compare.html", images=items, left=left, right=right)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def run_gui(host="127.0.0.1", port=5000, debug=None):
    """Launch the Flask development server for the QuickLook GUI.

    Parameters
    ----------
    host : str
        Interface to bind to (default ``127.0.0.1``, localhost only).
    port : int
        Port to listen on (default ``5000``).
    debug : bool | None
        Enable the Werkzeug debugger and auto-reloader. ``None`` (the default)
        consults the ``QUICKLOOK_DEBUG`` environment variable. The debugger
        exposes an interactive console (arbitrary code execution) to anyone who
        can reach the port, so it is opt-in.
    """
    if debug is None:
        debug = os.environ.get("QUICKLOOK_DEBUG", "").lower() in ("1", "true", "yes")

    # Route output only while the server is running. Installing these wrappers
    # at module import time captured temporary streams owned by pytest, Typer,
    # and other embedders; once those streams closed, later writes failed with
    # ``ValueError: I/O operation on closed file``.
    previous_stdout = sys.stdout
    previous_stderr = sys.stderr
    _tls_stdout.set_fallback(previous_stdout)
    _tls_stderr.set_fallback(previous_stderr)
    sys.stdout = _tls_stdout
    sys.stderr = _tls_stderr
    try:
        app.run(host=host, port=port, debug=debug, threaded=True)
    finally:
        if sys.stdout is _tls_stdout:
            sys.stdout = previous_stdout
        if sys.stderr is _tls_stderr:
            sys.stderr = previous_stderr


def main():
    """Entry point for the ``ql-gui`` console script."""
    run_gui()


if __name__ == "__main__":
    main()
