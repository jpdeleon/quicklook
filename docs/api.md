# API Reference

## TessQuickLook

The main class that runs the full pipeline.

```python
from quicklook import TessQuickLook
```

### Constructor

```python
TessQuickLook(
    target_name: str,
    sector: int = -1,
    pipeline: str = "SPOC",
    flux_type: str = "pdcsap",
    exptime: float = None,
    pg_method: str = "gls",
    flatten_method: str = "biweight",
    gp_kernel: str = "periodic_auto",
    gp_kernel_size: float = 5,
    window_length: float = None,
    edge_cutoff: float = 0.1,
    sigma_clip_raw: tuple = None,
    sigma_clip_flat: tuple = None,
    custom_ephem: list = None,
    mask_ephem: bool = False,
    Porb_limits: tuple = None,
    archival_survey: str = "dss1",
    show_plot: bool = True,
    verbose: bool = True,
    savefig: bool = False,
    savetls: bool = False,
    overwrite: bool = False,
    suffix: str = None,
    quality_bitmask: str = "default",
    outdir: str = ".",
)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `target_name` | str | required | Target name (TOI, TIC, or common name) |
| `sector` | int | `-1` | TESS sector (-1 = latest available) |
| `pipeline` | str | `"SPOC"` | Light curve pipeline |
| `flux_type` | str | `"pdcsap"` | Flux column: `"pdcsap"` or `"sap"` |
| `exptime` | float | `None` | Exposure time in seconds (auto-detected if None) |
| `pg_method` | str | `"gls"` | Periodogram method: `"gls"`, `"ls"`, `"bls"` |
| `flatten_method` | str | `"biweight"` | Detrending method (any wotan method) |
| `gp_kernel` | str | `"periodic_auto"` | GP kernel (if flatten_method is `"gp"`) |
| `gp_kernel_size` | float | `5` | GP kernel size |
| `window_length` | float | `None` | Detrending window in days (auto if None) |
| `edge_cutoff` | float | `0.1` | Days to cut from orbit edges |
| `sigma_clip_raw` | tuple | `None` | `(lo, hi)` sigma clip on raw light curve |
| `sigma_clip_flat` | tuple | `None` | `(lo, hi)` sigma clip on flattened light curve |
| `custom_ephem` | list | `None` | Custom ephemeris: `[Tc, Tc_err, P, P_err, dur, dur_err]` |
| `mask_ephem` | bool | `False` | Mask known transits before detrending |
| `Porb_limits` | tuple | `None` | `(min, max)` TLS period search range in days |
| `archival_survey` | str | `"dss1"` | Sky survey for image overlay |
| `show_plot` | bool | `True` | Display the figure interactively |
| `verbose` | bool | `True` | Print progress messages |
| `savefig` | bool | `False` | Save figure to disk |
| `savetls` | bool | `False` | Save TLS results to HDF5 |
| `overwrite` | bool | `False` | Overwrite existing output files |
| `suffix` | str | `None` | Suffix for output filenames |
| `quality_bitmask` | str | `"default"` | TESS quality bitmask |
| `outdir` | str | `"."` | Output directory |

### Key attributes

After instantiation, these attributes are available:

| Attribute | Type | Description |
|-----------|------|-------------|
| `raw_lc` | `lightkurve.LightCurve` | Raw light curve |
| `flat_lc` | `lightkurve.LightCurve` | Flattened light curve |
| `trend_lc` | `lightkurve.LightCurve` | Trend component |
| `tls_results` | `transitleastsquares.results` | TLS search results |
| `Prot_ls` | float | Stellar rotation period from GLS (days) |
| `tpf` | `lightkurve.TargetPixelFile` | Target pixel file |
| `sector` | int | TESS sector number |
| `pipeline` | str | Pipeline name (lowercase) |
| `cadence` | str | `"short"` or `"long"` |
| `exptime` | float | Exposure time in seconds |
| `gaiaid` | int | Gaia DR2/DR3 source ID |
| `tic` | int | TIC ID |

### Methods

#### `plot_tql()`

Generate the 9-panel diagnostic figure.

```python
fig = ql.plot_tql()
```

Returns a `matplotlib.figure.Figure`. If `savefig=True`, the figure is also saved to disk.

## Exceptions

```python
from quicklook.exceptions import QuickLookError, NoDataError, InvalidInputError, PipelineError
```

| Exception | Raised when |
|-----------|-------------|
| `QuickLookError` | Base exception for all pipeline errors |
| `NoDataError` | No light curves, TPFs, or sectors found for the target |
| `InvalidInputError` | Invalid user input (bad ephemeris, unsupported parameters) |
| `PipelineError` | Unrecoverable processing error during the pipeline |

## Utility functions

```python
from quicklook.utils import get_exofop_json, get_available_sectors
```

### `get_exofop_json(target_name)`

Query ExoFOP-TESS and return the full JSON response as a dict.

```python
data = get_exofop_json("TOI-1234")
```

### `get_available_sectors(target_name, pipeline="SPOC")`

Return a list of TESS sectors with available light curves for the target.

```python
sectors = get_available_sectors("TOI-1234", pipeline="SPOC")
# [27, 56]
```

## GLS periodogram

```python
from quicklook import Gls
```

The Generalized Lomb-Scargle periodogram implementation. Used internally by `TessQuickLook` but also available for standalone use.

```python
gls = Gls(data, Pbeg=0.5, Pend=15.0, verbose=False)
print(f"Best period: {gls.best['P']:.3f} days")
print(f"FAP: {gls.FAP():.2e}")
```
