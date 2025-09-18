import os
import io
import sys
from flask import Flask, render_template, request
from quicklook.tql import TessQuickLook
from pathlib import Path

app = Flask(__name__)
OUTPUT_DIR = os.path.join(app.static_folder, "outputs")


def run_quicklook(name, sector, pipeline, fluxtype):
    """
    Run TessQuickLook and return figure, PNG path, H5 path, and log output.
    """
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    tql = TessQuickLook(
        target_name=name,
        sector=int(sector),
        pipeline=pipeline,
        flux_type=fluxtype,
        savefig=True,
        savetls=True,
        outdir=OUTPUT_DIR,
        overwrite=True,
        verbose=True,
    )

    # Capture log by redirecting stdout temporarily if desired

    log_buffer = io.StringIO()
    sys_stdout = sys.stdout
    sys.stdout = log_buffer

    fig, png_file, h5_file = tql.plot_tql(return_fig_and_paths=True)

    sys.stdout = sys_stdout
    log_output = log_buffer.getvalue()

    return fig, str(png_file), str(h5_file), log_output


@app.route("/", methods=["GET", "POST"])
def index():
    image_url = None
    h5_url = None
    output_log = None

    if request.method == "POST":
        name = request.form.get("name")
        sector = request.form.get("sector")
        pipeline = request.form.get("pipeline")
        fluxtype = request.form.get("fluxtype")

        fig, png_file, h5_file, output_log = run_quicklook(
            name, sector, pipeline, fluxtype
        )

        if Path(png_file).exists():
            image_url = f"/static/outputs/{Path(png_file).name}"
        if Path(h5_file).exists():
            h5_url = f"/static/outputs/{Path(h5_file).name}"

    return render_template(
        "index.html", image=image_url, h5=h5_url, log=output_log
    )


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    app.run(debug=True)


if __name__ == "__main__":
    main()
