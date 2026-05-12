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
from collections import OrderedDict
from pathlib import Path
from threading import Thread, Event

from flask import Flask, render_template, request, jsonify
from flask_sock import Sock
from loguru import logger
from quicklook.tql import TessQuickLook
from quicklook.cli.ql import sanitize_target_name
from quicklook.exceptions import QuickLookError
from quicklook.utils import get_available_sectors

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
sys.stdout = _tls_stdout
sys.stderr = _tls_stderr


# ---------------------------------------------------------------------------
# SQLite job history
# ---------------------------------------------------------------------------
def _init_db():
    conn = sqlite3.connect(str(DB_PATH))
    conn.execute(
        """CREATE TABLE IF NOT EXISTS job_history (
            name         TEXT PRIMARY KEY,
            status       TEXT NOT NULL,
            error        TEXT DEFAULT '',
            params       TEXT DEFAULT '{}',
            submitted_at REAL,
            finished_at  REAL,
            step_times   TEXT DEFAULT '{}'
        )"""
    )
    conn.commit()
    conn.close()


_init_db()


def _save_job_history(
    name, status, error="", params=None, submitted_at=None, finished_at=None, step_times=None
):
    conn = sqlite3.connect(str(DB_PATH))
    conn.execute(
        "INSERT OR REPLACE INTO job_history "
        "(name, status, error, params, submitted_at, finished_at, step_times) "
        "VALUES (?, ?, ?, ?, ?, ?, ?)",
        (
            name,
            status,
            error,
            json.dumps(params or {}),
            submitted_at,
            finished_at,
            json.dumps(step_times or {}),
        ),
    )
    conn.commit()
    conn.close()


_avg_step_cache = {"data": {}, "ts": 0}
_avg_step_lock = threading.Lock()


def _get_avg_step_times():
    """Return average duration (seconds) per pipeline step from recent successful jobs."""
    now = time.time()
    with _avg_step_lock:
        if now - _avg_step_cache["ts"] < 60:
            return dict(_avg_step_cache["data"])
    conn = sqlite3.connect(str(DB_PATH))
    rows = conn.execute(
        "SELECT step_times FROM job_history WHERE status='done' "
        "ORDER BY finished_at DESC LIMIT 20"
    ).fetchall()
    conn.close()
    if not rows:
        return {}
    totals: dict[str, float] = {}
    counts: dict[str, int] = {}
    for (st_json,) in rows:
        st = json.loads(st_json) if st_json else {}
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

# Pipeline step definitions for progress tracking.
PIPELINE_STEPS = [
    (r"Generating quicklook", "Initializing"),
    (r"Catalog names:|TIC \d+", "Resolving target"),
    (r"Available sectors:|All available lightcurves", "Searching lightcurves"),
    (r"Downloading|search_lightcurve|Using .+ TPF", "Downloading data"),
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
    # Skip if cancelled while waiting in queue.
    if cancel_event.is_set():
        return

    # Mark as running now that the worker has picked it up.
    with jobs_lock:
        jobs[name]["status"] = "running"

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
        log_fh = open(log_file, "w", buffering=1, encoding="utf-8")
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
            exptime=kwargs.get("exptime"),
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
            custom_ephem=custom_ephem,
            archival_survey=kwargs.get("survey", "dss1"),
            savefig=kwargs.get("save", False),
            savetls=kwargs.get("save", False),
            verbose=kwargs.get("verbose", True),
            overwrite=kwargs.get("overwrite", False),
            outdir=outdir,
            suffix=suffix,
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
        jobs[name]["status"] = "error"
        jobs[name]["error"] = str(e)
        logger.error(f"Job '{name}' failed unexpectedly: {e}")
    finally:
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
    return [{"path": f"/static/outputs/{img.name}", "name": img.name} for img in images[:limit]]


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
        "mask_ephem": _is_truthy(form.get("mask_ephem")),
        "custom_ephem": form.get("custom_ephem"),
        "survey": form.get("survey", "dss1"),
        "outdir": form.get("outdir"),
        "suffix": form.get("suffix") if form.get("suffix") else "",
        "save": _is_truthy(form.get("save")),
        "verbose": _is_truthy(form.get("verbose")),
        "overwrite": _is_truthy(form.get("overwrite")),
    }


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
    return render_template("index.html", recent_results=recent_results)


@app.route("/submit", methods=["POST"])
def submit_job():
    """AJAX endpoint for single-job submission."""
    form = request.form if request.form else (request.json or {})
    name = (form.get("name") or "").strip()
    if not name:
        return jsonify({"ok": False, "reason": "Target name is required"}), 400

    name = sanitize_target_name(name)
    params = _parse_params(form)

    with jobs_lock:
        if name in jobs and jobs[name]["status"] not in ("done", "error", "cancelled"):
            return jsonify({"ok": False, "reason": f"Job '{name}' is already running"}), 409

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
        name = sanitize_target_name(raw_name)
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
            "image": f"/static/outputs/{images[0].name}" if images else None,
            "h5": f"/static/outputs/{h5s[0].name}" if h5s else None,
        }
    )


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


@app.route("/gallery")
def gallery():
    search_name = request.args.get("search", "")
    sort_by = request.args.get("sort", "date_desc")
    page = max(1, int(request.args.get("page", 1)))
    per_page = 24

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

    items = [{"path": f"/static/outputs/{img.name}", "name": img.name} for img in page_images]
    return render_template(
        "gallery.html",
        images=items,
        search_name=search_name,
        sort_by=sort_by,
        page=page,
        total_pages=total_pages,
        total=total,
    )


@app.route("/compare")
def compare():
    """Side-by-side comparison view."""
    images = sorted(OUTPUT_DIR.glob("*.png"), key=os.path.getmtime, reverse=True)
    items = [{"path": f"/static/outputs/{img.name}", "name": img.name} for img in images]
    left = request.args.get("left", "")
    right = request.args.get("right", "")
    return render_template("compare.html", images=items, left=left, right=right)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    app.run(debug=True, threaded=True)


if __name__ == "__main__":
    main()
