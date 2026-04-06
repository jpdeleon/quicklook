import json
import os
import re
import sys
import time
from collections import OrderedDict
from pathlib import Path
from threading import Thread, Event

from flask import Flask, render_template, request, jsonify
from flask_sock import Sock
from loguru import logger
from quicklook.tql import TessQuickLook
from quicklook.cli.ql import sanitize_target_name

# Directories
BASE_DIR = Path(__file__).parent.resolve()
OUTPUT_DIR = BASE_DIR / "static" / "outputs"
LOG_DIR = OUTPUT_DIR / "logs"
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)

app = Flask(
    __name__, static_folder=str(BASE_DIR / "static"), template_folder=str(BASE_DIR / "templates")
)
sock = Sock(app)

# Job queue: ordered dict preserving submission order.
# Each entry: {status, log_file, error, cancel_event, thread, submitted_at}
jobs = OrderedDict()

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


# --- Background job ---
def run_quicklook_background(name, cancel_event, **kwargs):
    name = sanitize_target_name(name)
    log_file = LOG_DIR / f"{name}.log"
    jobs[name]["log_file"] = str(log_file)

    # Add loguru sink to capture logger output into the log file
    sink_id = logger.add(str(log_file), format="{message}", level="DEBUG", filter="quicklook")

    with open(log_file, "w", buffering=1) as f:
        sys_stdout = sys.stdout
        sys_stderr = sys.stderr
        sys.stdout = f
        sys.stderr = f
        try:
            if cancel_event.is_set():
                raise InterruptedError("Job cancelled before start")
            tql = TessQuickLook(
                target_name=name,
                sector=int(kwargs.get("sector", -1)),
                pipeline=kwargs.get("pipeline", "spoc"),
                flux_type=kwargs.get("fluxtype", "pdcsap"),
                exptime=kwargs.get("exptime"),
                quality_bitmask=kwargs.get("quality_bitmask", "default"),
                flatten_method=kwargs.get("flatten_method", "biweight"),
                pg_method=kwargs.get("pg_method", "gls"),
                edge_cutoff=kwargs.get("edge_cutoff", 0.1),
                savefig=kwargs.get("save", False),
                savetls=kwargs.get("save", False),
                verbose=kwargs.get("verbose", True),
                overwrite=kwargs.get("overwrite", False),
                outdir=str(OUTPUT_DIR),
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
        except SystemExit as e:
            # TessQuickLook calls sys.exit() on certain errors
            jobs[name]["status"] = "error"
            jobs[name]["error"] = f"Pipeline exited (code {e.code})"
        except Exception as e:
            jobs[name]["status"] = "error"
            jobs[name]["error"] = str(e)
            logger.error(f"Job failed: {e}")
        finally:
            sys.stdout = sys_stdout
            sys.stderr = sys_stderr
            logger.remove(sink_id)


def _get_recent_results(limit=12):
    """Return the most recent output PNGs as dicts with path and filename."""
    images = sorted(OUTPUT_DIR.glob("*.png"), key=os.path.getmtime, reverse=True)
    return [{"path": f"/static/outputs/{img.name}", "name": img.name} for img in images[:limit]]


def _submit_job(name, params):
    """Create and start a background job. Returns the sanitized name."""
    name = sanitize_target_name(name)
    cancel_event = Event()
    jobs[name] = {
        "status": "running",
        "log_file": "",
        "error": "",
        "cancel_event": cancel_event,
        "submitted_at": time.time(),
    }
    # Move to end so newest jobs appear last
    jobs.move_to_end(name)

    thread = Thread(
        target=run_quicklook_background,
        kwargs={"name": name, "cancel_event": cancel_event, **params},
        daemon=True,
    )
    jobs[name]["thread"] = thread
    thread.start()
    return name


# --- Routes ---
@app.route("/", methods=["GET", "POST"])
def index():
    recent_results = _get_recent_results()

    if request.method == "POST":
        form = request.form
        name = sanitize_target_name(form.get("name"))

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

        params = {
            "sector": safe_int(form.get("sector", -1), -1),
            "pipeline": form.get("pipeline"),
            "fluxtype": form.get("fluxtype"),
            "exptime": safe_int(form.get("exptime") or None, None),
            "quality_bitmask": form.get("quality_bitmask", "default"),
            "flatten_method": form.get("flatten_method", "biweight"),
            "pg_method": form.get("pg_method", "gls"),
            "edge_cutoff": safe_float(form.get("edge_cutoff", 0.1), 0.1),
            "save": "save" in form,
            "verbose": "verbose" in form,
            "overwrite": "overwrite" in form,
        }

        # Allow re-submission only if previous run finished
        if name not in jobs or jobs[name]["status"] in ("done", "error", "cancelled"):
            _submit_job(name, params)

        return render_template("index.html", target=name, recent_results=recent_results)

    return render_template("index.html", recent_results=recent_results)


@app.route("/jobs-json")
def jobs_json():
    """Return status of all jobs (for queue display)."""
    result = []
    for name, info in jobs.items():
        result.append(
            {
                "name": name,
                "status": info["status"],
                "error": info.get("error", ""),
            }
        )
    return jsonify(result)


@app.route("/cancel/<target>", methods=["POST"])
def cancel(target):
    """Cancel a running job."""
    info = jobs.get(target)
    if info and info["status"] == "running":
        info["cancel_event"].set()
        info["status"] = "cancelled"
        return jsonify({"ok": True})
    return jsonify({"ok": False, "reason": "not running"}), 400


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
        # Also delete matching H5
        h5_name = png_path.stem + "_tls.h5"
        # Handle case where stem already ends with a sector/cadence suffix
        # e.g. TOI-1234_s42_pdcsap_sc.png -> TOI-1234_s42_pdcsap_sc_tls.h5
        h5_path = OUTPUT_DIR / h5_name
        if h5_path.exists():
            h5_path.unlink()
            deleted.append(h5_name)

    if deleted:
        return jsonify({"ok": True, "deleted": deleted})
    return jsonify({"ok": False, "reason": "file not found"}), 404


@app.route("/results-json/<target>")
def results_json(target):
    """Return paths to output PNG and H5 files for a completed job."""
    images = sorted(OUTPUT_DIR.glob(f"*{target}*.png"), key=os.path.getmtime, reverse=True)
    h5s = sorted(OUTPUT_DIR.glob(f"*{target}*_tls.h5"), key=os.path.getmtime, reverse=True)
    return jsonify(
        {
            "image": f"/static/outputs/{images[0].name}" if images else None,
            "h5": f"/static/outputs/{h5s[0].name}" if h5s else None,
        }
    )


@sock.route("/ws/<target>")
def ws_log(ws, target):
    """WebSocket endpoint for live log streaming and cancel commands."""
    info = jobs.get(target)
    if not info:
        ws.send(json.dumps({"type": "finish", "status": "unknown"}))
        return

    log_path = Path(info["log_file"]) if info.get("log_file") else None
    last_size = 0

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
        except (json.JSONDecodeError, TypeError):
            pass
        except Exception:
            # Connection closed or timeout — timeout returns None which is fine
            pass

        status = info.get("status", "unknown")

        # Read new log content
        new_lines = ""
        full_log = ""
        if log_path and log_path.exists():
            content = log_path.read_text()
            full_log = content
            if len(content) > last_size:
                new_lines = content[last_size:]
                last_size = len(content)

        # Detect progress
        if full_log:
            step_idx, step_total = _detect_step(full_log)
            step_label = PIPELINE_STEPS[step_idx][1] if step_idx >= 0 else "Starting"
            pct = int(((step_idx + 1) / step_total) * 100) if step_idx >= 0 else 0
        else:
            step_label = "Starting"
            pct = 0

        # Send update
        try:
            ws.send(
                json.dumps(
                    {
                        "type": "update",
                        "log": new_lines,
                        "status": status,
                        "progress": pct,
                        "step": step_label,
                    }
                )
            )
        except Exception:
            break  # Client disconnected

        if status in ("done", "error", "cancelled"):
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
            except Exception:
                pass
            break

        time.sleep(1.5)


@app.route("/gallery")
def gallery():
    search_name = request.args.get("search", "")
    images = sorted(OUTPUT_DIR.glob("*.png"), key=os.path.getmtime, reverse=True)
    if search_name:
        images = [img for img in images if search_name.lower() in img.stem.lower()]
    items = [{"path": f"/static/outputs/{img.name}", "name": img.name} for img in images]
    return render_template("gallery.html", images=items, search_name=search_name)


# --- Main ---
def main():
    app.run(debug=True, threaded=True)


if __name__ == "__main__":
    main()
