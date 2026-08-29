# CLI Reference

QuickLook provides a unified Typer-based CLI (`ql`) with subcommands.

## Global help

```bash
ql --help
```

Shows available subcommands: `run`, `read-tls`, `rank-tls`.

---

## `ql run` -- Run QuickLook analysis

```bash
ql run --name TARGET [options]
```

### Required

| Flag | Description |
|------|-------------|
| `--name NAME` | Target name (TOI, TIC, or common name) |

### Observation selection

| Flag | Default | Description |
|------|---------|-------------|
| `--sector SECTOR` | `-1` (latest) | TESS sector number |
| `--pipeline PIPELINE` | `spoc` | Light curve pipeline: `spoc`, `tess-spoc`, `qlp`, `cdips`, `pathos`, `tglc`, `tasoc`, `t16` |
| `--fluxtype TYPE` | `pdcsap` | Light-curve type. SPOC: `pdcsap` or `sap`. TGLC: `aperture`, `psf`, or `auto` |
| `--exptime SECONDS` | auto | Exposure time in seconds |
| `--quality-bitmask MASK` | `default` | Quality mask: `none`, `default`, `hard`, `hardest` |

### Detrending

| Flag | Default | Description |
|------|---------|-------------|
| `--flatten-method METHOD` | `biweight` | Detrending method: any [wotan](https://github.com/hippke/wotan) method, or `notch`/`locor` (see [Notch and LOCoR](pipeline.md#notch-and-locor)) |
| `--window-length DAYS` | auto | Detrending window length in days (0 = auto-optimize); for `locor` this is the rotation period instead (0/unset = auto-estimate via GLS) |
| `--edge-cutoff DAYS` | `0.1` | Cut this many days from each orbit edge |
| `--sigma-clip-raw LO HI` | `10 5` | Sigma clipping on raw light curve before flattening |
| `--sigma-clip-flat LO HI` | none | Sigma clipping on flattened light curve |
| `--gp-kernel KERNEL` | `periodic_auto` | GP kernel: `periodic_auto`, `periodic`, `squared_exp` |
| `--gp-kernel-size N` | `5` | GP kernel size |

### Period search

| Flag | Default | Description |
|------|---------|-------------|
| `--pg-method METHOD` | `gls` | Periodogram method: `gls`, `ls`, or `bls` |
| `--period-limits MIN MAX` | auto | TLS search range in days |
| `--phase-xlim HALF_WIDTH` | auto | Phase half-width for transit/eclipse panels |

### Ephemeris

| Flag | Default | Description |
|------|---------|-------------|
| `--custom-ephem Tc Tcerr P Perr Tdur Tdurerr` | none | Custom ephemeris (6 values in days) |
| `--mask-ephem` | off | Mask known transits before detrending |
| `--use-priors` | off | Use ExoFOP stellar params as TLS priors |

### Output

| Flag | Default | Description |
|------|---------|-------------|
| `--outdir DIR` | `.` | Output directory |
| `--suffix TEXT` | none | Suffix appended to output filenames |
| `--survey NAME` | `dss1` | Archival image survey for overlay |
| `--show-simbad` | off | Overplot nearby SIMBAD objects on the archival image |
| `--save` | off | Save figure (.png) and TLS results (.h5) |
| `--verbose` | off | Print detailed progress |
| `--overwrite` | off | Overwrite existing output files |

### Multi-sector / multi-pipeline mode

| Flag | Default | Description |
|------|---------|-------------|
| `--each-sector` | off | Run on every available sector individually |
| `--each-pipeline` | off | Run on every available pipeline for the latest sector (mutually exclusive with `--each-sector`) |
| `-j, --jobs N` | `1` | Number of parallel jobs for `--each-sector` |
| `--nice INC` | unchanged | Lower CPU priority by this increment (POSIX; e.g. `19` = lowest) |
| `--cores N` | auto | CPU cores used by TLS per run |

### Examples

```bash
# Basic run on latest sector
ql run --name TOI-1234 --save --verbose

# Specific sector and pipeline
ql run --name TOI-125.01 --sector 2 --pipeline qlp

# Custom detrending
ql run --name HD209458 --flatten-method cosine --window-length 0.3

# Restrict period search range
ql run --name TOI-125.01 --period-limits 1 5

# Run all sectors, save results
ql run --name TOI-125.01 --each-sector --save

# Run all sectors with 4 parallel workers
ql run --name TOI-125.01 --each-sector -j 4 --save

# Run every available pipeline on the latest sector
ql run --name TOI-125.01 --each-pipeline --save

# TGLC PSF photometry with nearby SIMBAD objects overlaid
ql run --name TOI-125.01 --pipeline tglc --fluxtype psf --show-simbad --save

# Pipe output to a log file
ql run --name TIC-52368076 --verbose --save | tee output.log
```

> Underscore-style flags (`--flatten_method`) are also accepted as equivalents
> to hyphen-style flags (`--flatten-method`).

---

## `ql read-tls` -- Read TLS results

Scan a directory of `.h5` TLS output files and produce a summary CSV:

```bash
ql read-tls OUTPUT_DIR
```

This creates `OUTPUT_DIR_tls.csv` with columns for period, depth, SDE, and other
TLS metrics. Useful for ranking candidates after batch processing.

| Flag | Default | Description |
|------|---------|-------------|
| `--param, -p PARAM` | `SDE_tls` | Parameter to sort by |

---

## `ql rank-tls` -- Rank candidates

Filter and rank candidates from the CSV produced by `read-tls`:

```bash
ql rank-tls OUTPUT_DIR --output-dir ranked
```

This copies the plots of candidates that pass the filters into the `ranked/`
directory, sorted by Signal Detection Efficiency (SDE).

| Flag | Default | Description |
|------|---------|-------------|
| `--csv-path PATH` | auto | Path to the existing CSV (auto-generated if omitted) |
| `--output-dir DIR` | `{input_dir}/temp` | Output directory for ranked plots |
| `--column COL` | `SDE_tls` | Column to sort by |
| `--remove-filter` | off | Skip filtering (rank all candidates) |
| `--use-ascending` | off | Sort ascending instead of descending |
| `--min-sde N` | `5.0` | Minimum SDE threshold |
| `--min-depth N` | `1.0` | Minimum transit depth (ppt) |
| `--max-depth N` | `100.0` | Maximum transit depth (ppt) |
