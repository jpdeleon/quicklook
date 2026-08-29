# Pipeline Overview

QuickLook produces a 9-panel diagnostic figure. This page explains what each panel shows and how the pipeline works.

## Pipeline steps

```
Target name
    |
    v
1. Resolve target (ExoFOP-TESS -> TIC ID, coordinates, known ephemeris)
    |
    v
2. Download light curve (MAST via lightkurve)
    |
    v
3. Sigma-clip outliers from raw light curve
    |
    v
4. Flatten/detrend (wotan biweight filter by default)
    |
    v
5. Generalized Lomb-Scargle periodogram (stellar rotation)
    |
    v
6. Transit Least Squares search (planet/companion detection)
    |
    v
7. Download target pixel file (TPF)
    |
    v
8. Query Gaia catalog for nearby sources
    |
    v
9. Generate 9-panel diagnostic figure
```

## The 9-panel figure

### Top row -- Stellar variability

**Panel 1: Raw light curve**

The raw PDCSAP (or SAP) flux from the TESS pipeline, shown in blue. The red line is the best-fit trend from the detrending step (wotan biweight filter). This panel reveals stellar variability, systematics, and any obvious transit-like dips.

**Panel 2: Lomb-Scargle periodogram**

The Generalized Lomb-Scargle (GLS) periodogram of the raw light curve (with known transits masked). The peak corresponds to the stellar rotation period. High peaks indicate active or variable stars. The periodogram method can be changed to standard Lomb-Scargle (`ls`) or Box Least Squares (`bls`).

**Panel 3: Phase-folded rotation**

The raw light curve folded at the peak GLS period, showing the rotational modulation pattern. A clean sinusoidal shape suggests starspots rotating in and out of view.

### Middle row -- Transit detection

**Panel 4: Flattened light curve**

The raw light curve divided by the trend (panel 1), producing a normalized, detrended light curve centered on 1.0. Red markers highlight data points inside detected transits (from the TLS search in panel 5). This panel is where transits are most visible.

**Panel 5: TLS periodogram**

The [Transit Least Squares](https://ui.adsabs.harvard.edu/abs/2019A%26A...623A..39H/abstract) periodogram, which searches for periodic box-shaped (transit-like) dips. The y-axis is the Signal Detection Efficiency (SDE). The peak period is the best-fit orbital period. An SDE above ~8 is generally considered a significant detection.

**Panel 6: Archival image + Gaia overlay**

The TESS aperture (blue polygon) overlaid on an archival sky survey image (DSS by default). Orange and red circles show nearby Gaia sources, scaled by brightness. With `-show_simbad` (CLI) or the GUI checkbox, nearby non-stellar SIMBAD objects within the field of view are also labelled by their condensed object type (e.g. `EclBin`), which helps flag known eclipsing binaries or other variables near the target. DSS cutouts are cached under `~/.astropy/cache` so repeated runs of the same field do not re-download the image. This panel helps identify:

- Blended neighbors that could be the true source of the signal
- Background eclipsing binaries contaminating the aperture
- Crowded fields where photometric dilution is significant

### Bottom row -- Transit characterization

**Panel 7: Phase-folded transit (odd/even)**

The flattened light curve folded at the TLS best-fit period and zoomed around the transit. Odd transits (red) and even transits (blue) are shown separately. The black line is the best-fit transit model. If odd and even depths differ significantly, the signal is likely an eclipsing binary rather than a planet.

**Panel 8: Secondary eclipse check**

The phase-folded light curve centered at phase 0.5 (halfway between transits). A dip here indicates a secondary eclipse, which is a strong sign of a self-luminous companion (eclipsing binary or high-albedo brown dwarf). The dashed horizontal line shows the primary transit depth for reference.

**Panel 9: Summary panel**

A text summary with:

- Target identifiers (TIC, TOI, Gaia ID, coordinates)
- Stellar parameters from ExoFOP (radius, mass, Teff)
- TLS results (period, depth, duration, SDE)
- Derived companion parameters (radius ratio, estimated radius)
- SIMBAD object classification
- Pipeline metadata (sector, cadence, runtime)

## Detrending methods

QuickLook uses [wotan](https://github.com/hippke/wotan) for most detrending methods. The default is `biweight`, but any wotan method can be used:

| Method | Description | Best for |
|--------|-------------|----------|
| `biweight` | Tukey's biweight (robust) | General use (default) |
| `cosine` | Cosine filtering | Smooth trends |
| `median` | Sliding median | Simple systematics |
| `mean` | Sliding mean | Simple systematics |
| `gp` | Gaussian process | Complex variability (slow for short cadence) |

!!! warning
    GP flattening (`--flatten_method gp`) is not recommended for short-cadence data (exposure time <= 120s) due to computational cost.

### Notch and LOCoR

Two additional methods are implemented directly in QuickLook (`quicklook/notch_locor.py`), adapted from the Notch filter and Locally Optimized Combination of Rotations (LOCoR) algorithms of [Rizzuto et al. 2017](https://arxiv.org/abs/1709.09670) ([arizzuto/Notch_and_LOCoR](https://github.com/arizzuto/Notch_and_LOCoR)):

| Method | Description | `--window-length` meaning |
|--------|-------------|----------------------------|
| `notch` | Sliding-window quadratic + notch-box fit, BIC-selected against a plain quadratic per window, so a real transit dip is preserved rather than divided out | Notch window size in days (default 0.5) |
| `locor` | Phase-folds on the rotation period and fits each cycle as an optimal linear combination of neighboring cycles | Rotation period in days (auto-estimated via GLS if left unset) |

!!! note
    The reference implementation pins its per-window fits to `mpyfit`, a C-backed fitter that isn't on PyPI and needs a local compiler. QuickLook's port solves the same windowed model with `scipy.optimize.least_squares` instead, so results follow the same algorithm but aren't bit-for-bit identical to the upstream repository.

## Output files

When `-save` is used, QuickLook produces:

| File | Description |
|------|-------------|
| `TARGET_sSECTOR_FLUX_CADENCE.png` | The 9-panel diagnostic figure |
| `TARGET_sSECTOR_FLUX_CADENCE_tls.h5` | HDF5 file with full TLS results |

Example: `WASP-21_s56_pdcsap_sc.png` and `WASP-21_s56_pdcsap_sc_tls.h5`

For TGLC runs, the flux token encodes the photometry method: `--fluxtype aperture` writes `tglc_aper` and `--fluxtype psf` writes `tglc_psf` (e.g. `TIC123_s11_tglc_psf_lc.png`). This keeps aperture and PSF runs of the same target/sector from overwriting each other on disk. SPOC stems and default (auto) TGLC stems are unchanged.

The HDF5 file contains the TLS periodogram, best-fit parameters, stellar parameters, and metadata. Use `ql read-tls` to extract a summary CSV from a directory of these files.
