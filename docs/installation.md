# Installation

## Requirements

- Python 3.10 or later
- A working internet connection (for downloading light curves from MAST)

## Basic install

Create a conda environment and install from PyPI:

```bash
conda create -n quicklook python=3.12
conda activate quicklook
pip install -U quicklook-package
```

This installs the `ql`, `read_tls`, and `rank_tls` command-line tools plus the Python API.

## Install with Web GUI

The web GUI requires additional dependencies (Flask, flask-sock). Install them with:

```bash
pip install -U "quicklook-package[gui]"
```

This adds the `ql-gui` command.

## Install for notebooks

To run QuickLook interactively in Jupyter:

```bash
pip install -U "quicklook-package[notebooks]"
```

## Development install

Clone the repository and install in editable mode with all extras:

```bash
git clone https://github.com/jpdeleon/quicklook.git
cd quicklook
pip install -e ".[dev,gui,notebooks]"
```

This installs testing tools (pytest, black, ruff, tox) and enables live code changes without reinstalling.

## Verifying the install

```bash
# Check the CLI
ql --help

# Check imports
python -c "from quicklook import TessQuickLook; print('OK')"
```

## Dependencies

QuickLook depends on these core packages (installed automatically):

| Package | Purpose |
|---------|---------|
| [lightkurve](https://docs.lightkurve.org/) | TESS/Kepler light curve search and download |
| [transitleastsquares](https://github.com/hippke/tls) | Transit detection periodogram |
| [wotan](https://github.com/hippke/wotan) | Light curve detrending/flattening |
| [astropy](https://www.astropy.org/) | Coordinates, units, FITS I/O |
| [matplotlib](https://matplotlib.org/) | Plotting |
| [numpy](https://numpy.org/) | Numerical computation |
| [h5py](https://www.h5py.org/) | HDF5 file I/O for TLS results |
| [reproject](https://reproject.readthedocs.io/) | Archival image reprojection |
| [loguru](https://loguru.readthedocs.io/) | Logging |
