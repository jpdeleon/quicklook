import os
import sys
from pathlib import Path
from flask import Flask, render_template, request, redirect, url_for, jsonify
from threading import Thread
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

jobs = {}  # Tracks running jobs and logs


# --- Background job ---
def run_quicklook_background(name, **kwargs):

    name = sanitize_target_name(name)
    log_file = LOG_DIR / f"{name}.log"
    jobs[name] = {"status": "running", "log_file": str(log_file)}

    with open(log_file, "w", buffering=1) as f:
        sys_stdout = sys.stdout
        sys.stdout = f
        try:
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
                # add other arguments here as needed...
                outdir=str(OUTPUT_DIR),
            )
            tql.plot_tql(return_fig_and_paths=True)
            jobs[name]["status"] = "done"
        except Exception as e:
            print(e, file=sys.stderr)
            jobs[name]["status"] = "error"
        finally:
            sys.stdout = sys_stdout


# --- Routes ---
@app.route("/", methods=["GET", "POST"])
def index():
    running_targets = {k: v for k, v in jobs.items() if v["status"] == "running"}

    if request.method == "POST":

        form = request.form
        name = sanitize_target_name(form.get("name"))
        sector = form.get("sector", -1)
        pipeline = form.get("pipeline")
        fluxtype = form.get("fluxtype")
        exptime = form.get("exptime") or None
        quality_bitmask = form.get("quality_bitmask", "default")
        flatten_method = form.get("flatten_method", "biweight")
        pg_method = form.get("pg_method", "gls")
        edge_cutoff = form.get("edge_cutoff", 0.1)
        save = "save" in form
        verbose = "verbose" in form
        overwrite = "overwrite" in form

        # Convert numeric fields properly
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

        sector = safe_int(sector, -1)
        exptime = safe_int(exptime, None)
        edge_cutoff = safe_float(edge_cutoff, 0.1)

        # Launch background job
        if name not in jobs or jobs[name]["status"] in ["done", "error"]:
            thread = Thread(
                target=run_quicklook_background,
                kwargs={
                    "name": name,
                    "sector": sector,
                    "pipeline": pipeline,
                    "fluxtype": fluxtype,
                    "exptime": exptime,
                    "quality_bitmask": quality_bitmask,
                    "flatten_method": flatten_method,
                    "pg_method": pg_method,
                    "edge_cutoff": edge_cutoff,
                    "save": save,
                    "verbose": verbose,
                    "overwrite": overwrite,
                },
            )
            thread.start()

        return render_template("index.html", running_targets=running_targets, target=name)

    return render_template("index.html", running_targets=running_targets)


@app.route("/log/<target>")
def log(target):
    info = jobs.get(target)
    if info:
        log_path = Path(info["log_file"])
        if log_path.exists():
            return log_path.read_text()
    return ""


@app.route("/status-json/<target>")
def status_json(target):
    info = jobs.get(target)
    if info and isinstance(info, dict):
        return jsonify(info)
    return jsonify({"status": "unknown"})


@app.route("/gallery", methods=["GET", "POST"])
def gallery():
    search_name = request.form.get("search") if request.method == "POST" else None
    images = list(OUTPUT_DIR.glob("*.png"))
    if search_name:
        images = [img for img in images if search_name.lower() in img.stem.lower()]
    images = [f"/static/outputs/{img.name}" for img in images]
    return render_template("gallery.html", images=images, search_name=search_name)


# --- Main ---
def main():
    app.run(debug=True, threaded=True)


if __name__ == "__main__":
    main()
