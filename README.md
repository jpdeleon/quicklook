# QuickLook

[![Documentation](https://readthedocs.org/projects/quicklook/badge/?version=latest)](https://quicklook.readthedocs.io/en/latest/)
[![CI](https://github.com/jpdeleon/quicklook/actions/workflows/build-and-test.yml/badge.svg)](https://github.com/jpdeleon/quicklook/actions/workflows/build-and-test.yml)
[![PyPI](https://img.shields.io/pypi/v/quicklook-package)](https://pypi.org/project/quicklook-package/)
[![Python](https://img.shields.io/pypi/pyversions/quicklook-package)](https://pypi.org/project/quicklook-package/)

`quicklook` is a Python pipeline that searches for transit signals in [TESS](https://tess.mit.edu/) light curves. Given a target name, it downloads the light curve, estimates the stellar rotation period, searches for transiting companions using the [Transit Least Squares](https://ui.adsabs.harvard.edu/abs/2019A%26A...623A..39H/abstract) (TLS) algorithm, and produces a publication-ready 9-panel diagnostic figure.

Although `quicklook` is optimized to find transiting exoplanets, it can also detect eclipsing binaries, variable stars, and other periodic signals.

## Features

- **Multi-pipeline support** -- SPOC, TESS-SPOC, QLP, CDIPS, PATHOS, TGLC, TASOC
- **Automated detrending** -- biweight, cosine, median, GP, and other [wotan](https://github.com/hippke/wotan) methods
- **Stellar rotation** -- Generalized Lomb-Scargle (GLS) periodogram
- **Transit detection** -- Transit Least Squares (TLS) periodogram
- **Neighbor check** -- Gaia source overlay on archival sky images
- **Batch processing** -- `--each-sector` mode, GNU parallel support, and candidate ranking tools
- **Web GUI** -- Flask-based interface with live progress, job queue, and gallery
- **HDF5 output** -- full TLS results saved for downstream filtering

## Installation

Create a conda environment and install from PyPI:

```bash
conda create -n quicklook python=3.12
conda activate quicklook
pip install -U quicklook-package
```

### Optional extras

```bash
# Web GUI
pip install -U "quicklook-package[gui]"

# Jupyter notebooks
pip install -U "quicklook-package[notebooks]"

# Development (testing, linting, formatting)
pip install -U "quicklook-package[dev]"
```

## Try it on Google Colab

<a href="https://colab.research.google.com/github/jpdeleon/quicklook/blob/main/notebook/examples.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

## Usage

### Command line

```bash
# Basic run on the latest TESS sector
ql --name WASP-21 -save -verbose

# Specific sector and pipeline
ql --name TOI-125.01 --sector 2 --pipeline qlp

# Custom detrending
ql --name TOI-125.01 --flatten_method cosine --window_length 0.3

# Restrict TLS period search range
ql --name TOI-125.01 --period_limits 1 5

# Run on every available sector
ql --name TOI-125.01 --each-sector -save

# Run all sectors with 4 parallel workers
ql --name TOI-125.01 --each-sector -j 4 -save
```

### Python API

```python
from quicklook import TessQuickLook

ql = TessQuickLook(
    target_name="WASP-21",
    sector=56,
    pipeline="SPOC",
    flux_type="pdcsap",
    verbose=True,
)

fig = ql.plot_tql()

# Access results
print(f"Rotation period: {ql.Prot_ls:.2f} days")
print(f"TLS period: {ql.tls_results.period:.4f} days")
print(f"TLS SDE: {ql.tls_results.SDE:.1f}")
```

### Web GUI

```bash
ql-gui
```

Open http://127.0.0.1:5000 in your browser. Enter a target, adjust parameters, and click **Run QuickLook**. Progress is streamed live via WebSocket. Supports single targets, batch submission, and each-sector mode.

![QuickLook Web GUI](docs/img/ql-gui.png)

## Output figure

```bash
ql --name WASP-21 -save -verbose
```

![Example output](tests/WASP-21_s83_pdcsap_sc.png)

The 9-panel figure shows:

| Panel | Content |
|-------|---------|
| 1 | Raw light curve + trend line |
| 2 | GLS periodogram (stellar rotation period) |
| 3 | Phase-folded light curve at rotation period |
| 4 | Flattened light curve + detected transits |
| 5 | TLS periodogram (orbital period) |
| 6 | TESS aperture + Gaia sources on archival image |
| 7 | Phase-folded transit (odd/even comparison) |
| 8 | Secondary eclipse check at phase 0.5 |
| 9 | Summary of stellar and companion parameters |

## CLI tools

| Command | Description |
|---------|-------------|
| `ql` | Run the full QuickLook pipeline on a target |
| `read_tls` | Extract TLS results from a directory of `.h5` files into a CSV |
| `rank_tls` | Filter and rank candidates by SDE from the CSV output |
| `ql-gui` | Launch the web GUI (requires `[gui]` extra) |

## Batch processing

Process a list of TIC IDs:

```bash
# Generate batch script
cat tic_ids.txt | while read tic; do
  echo "ql --name TIC$tic -save --outdir results | tee TIC$tic.log"
done > run_batch.sh

# Run in parallel with GNU parallel
cat run_batch.sh | parallel -j 4

# Extract and rank results
read_tls results/
rank_tls results/ --output_dir ranked

# Combine ranked plots into a PDF
pip install img2pdf
img2pdf ranked/*.png --output ranked.pdf
```

## Documentation

Full documentation is available at [quicklook.readthedocs.io](https://quicklook.readthedocs.io/en/latest/).

## License

See [LICENSE](LICENSE) for details.
