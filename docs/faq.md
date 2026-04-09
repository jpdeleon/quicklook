# FAQ

## General

### What targets can QuickLook analyze?

Any target observed by TESS that has a light curve on MAST. This includes TOI candidates, known exoplanets, eclipsing binaries, variable stars, and arbitrary TIC IDs.

### What is the SDE and what values are significant?

SDE (Signal Detection Efficiency) is the signal-to-noise metric from the TLS periodogram. Values above ~8 generally indicate a significant periodic transit-like signal. However, astrophysical false positives (eclipsing binaries, contamination) can also produce high SDE values.

### Can QuickLook analyze Kepler or K2 data?

Not currently. QuickLook is designed for TESS data. Kepler/K2 support may be added in the future.

## Errors and troubleshooting

### `NoDataError: No light curves found`

The target has no TESS light curves for the specified pipeline and sector. Try:

- A different pipeline (`--pipeline qlp` or `--pipeline tess-spoc`)
- Removing the sector constraint (`--sector -1`)
- Checking the target name spelling

### `No module named 'distutils'` (Python 3.12)

The `batman-package` dependency imports `distutils`, which was removed in Python 3.12. Install `setuptools`:

```bash
pip install setuptools
```

This is included automatically if you install `quicklook-package >= 1.5`.

### MAST HTTP 500 errors

The MAST archive sometimes experiences outages. These are transient -- wait a few minutes and try again. In CI, network tests automatically mark as `xfail` when MAST is down.

### `I/O operation on closed file` in the web GUI

This was caused by running multiple pipeline jobs concurrently. As of v1.5, the GUI serializes all jobs through a single worker queue. If you see this error, make sure you are running the latest version.

### `GLS power is NaN`

The Generalized Lomb-Scargle periodogram returned NaN, which means the data is unsuitable for GLS (e.g., constant flux, too few data points). QuickLook automatically falls back to astropy's Lomb-Scargle implementation when this happens.

### Plots look garbled or axes overlap

Try a different matplotlib backend:

```bash
export MPLBACKEND=Agg
ql --name TOI-1234 -save
```

## Pipeline parameters

### How do I choose the window length?

The detrending window length controls how aggressive the baseline removal is.

- **Too short**: removes real transit signals along with the trend
- **Too long**: leaves residual systematics in the flattened light curve

The default (`None`) auto-selects 3x the known transit duration, or 0.5 days if no ephemeris is available. Set `--window_length 0` to auto-optimize.

### When should I use GP flattening?

GP (`--flatten_method gp`) is useful for complex stellar variability that biweight filtering cannot capture (e.g., multi-mode pulsators). However, it is computationally expensive for short-cadence (2-min) data and is not recommended for routine use.

### What does `--mask_ephem` do?

It masks known transits (from ExoFOP or custom ephemeris) before detrending. This prevents the detrending algorithm from partially removing the transit signal, which can happen with aggressive window lengths.

## Output files

### How do I read the HDF5 output?

```python
import h5py

with h5py.File("TOI-1234_s56_pdcsap_sc_tls.h5", "r") as f:
    period = f["period"][()]
    power = f["power"][()]
    sde = f.attrs["SDE"]
    depth = f.attrs["depth"]
    print(f"Period: {period:.4f} d, SDE: {sde:.1f}, Depth: {depth:.6f}")
```

Or use the `read_tls` tool to extract a CSV summary from a directory of H5 files.
