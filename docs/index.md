# QuickLook

**QuickLook** is a Python pipeline for detecting transit signals in [TESS](https://tess.mit.edu/) light curves. Given a target name, it downloads the light curve, runs periodograms to measure the stellar rotation period, searches for transiting companions using the [Transit Least Squares](https://ui.adsabs.harvard.edu/abs/2019A%26A...623A..39H/abstract) (TLS) algorithm, and produces a publication-ready diagnostic figure.

## What it does

- Downloads TESS light curves from [MAST](https://mast.stsci.edu/) via [lightkurve](https://docs.lightkurve.org/)
- Estimates the stellar rotation period using a Generalized Lomb-Scargle (GLS) periodogram
- Searches for periodic transit signals using TLS
- Overlays Gaia sources on archival sky images to check for blended neighbors
- Produces a 9-panel diagnostic figure summarizing all results
- Saves TLS results to HDF5 for downstream ranking and filtering

## Ways to use it

| Method | Best for |
|--------|----------|
| [CLI (`ql`)](cli.md) | Quick single-target analysis from the terminal |
| [Web GUI (`ql-gui`)](gui.md) | Interactive exploration with live progress |
| [Python API](api.md) | Notebooks, scripting, custom workflows |
| [Batch processing](batch.md) | Surveys and candidate vetting at scale |

## Output figure

Running `ql --name WASP-21 -save` produces a 9-panel figure:

![Example output](https://raw.githubusercontent.com/jpdeleon/quicklook/main/tests/WASP-21_s83_pdcsap_sc.png)

See [Pipeline Overview](pipeline.md) for a detailed explanation of each panel.

## Try it now

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jpdeleon/quicklook/blob/main/notebook/examples.ipynb)

## Links

- [GitHub repository](https://github.com/jpdeleon/quicklook)
- [PyPI package](https://pypi.org/project/quicklook-package/)
- [Example notebook](https://github.com/jpdeleon/quicklook/blob/main/notebook/examples.ipynb)
