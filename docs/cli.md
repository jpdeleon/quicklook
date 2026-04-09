# CLI Reference

QuickLook provides three command-line tools: `ql`, `read_tls`, and `rank_tls`.

## `ql` -- Run QuickLook analysis

```bash
ql --name TARGET [options]
```

### Required

| Flag | Description |
|------|-------------|
| `--name NAME` | Target name (TOI, TIC, or common name) |

### Observation selection

| Flag | Default | Description |
|------|---------|-------------|
| `--sector SECTOR` | `-1` (latest) | TESS sector number |
| `--pipeline PIPELINE` | `spoc` | Light curve pipeline: `spoc`, `tess-spoc`, `qlp`, `cdips`, `pathos`, `tglc`, `tasoc` |
| `--fluxtype TYPE` | `pdcsap` | Flux column: `pdcsap` or `sap` |
| `--exptime SECONDS` | auto | Exposure time in seconds |
| `--quality_bitmask MASK` | `default` | Quality mask: `none`, `default`, `hard`, `hardest` |

### Detrending

| Flag | Default | Description |
|------|---------|-------------|
| `--flatten_method METHOD` | `biweight` | Detrending method (any [wotan](https://github.com/hippke/wotan) method) |
| `--window_length DAYS` | auto | Detrending window length in days (0 = auto-optimize) |
| `--edge_cutoff DAYS` | `0.1` | Cut this many days from each orbit edge |
| `--sigma_clip_raw LO HI` | `10 5` | Sigma clipping on raw light curve before flattening |
| `--sigma_clip_flat LO HI` | none | Sigma clipping on flattened light curve |

### Period search

| Flag | Default | Description |
|------|---------|-------------|
| `--pg_method METHOD` | `gls` | Periodogram method: `gls`, `ls`, or `bls` |
| `--period_limits MIN MAX` | auto | TLS search range in days (default: 0.5 to baseline/2) |

### Ephemeris

| Flag | Default | Description |
|------|---------|-------------|
| `--custom_ephem Tc Tcerr P Perr Tdur Tdurerr` | none | Custom ephemeris (6 values in days) |
| `-mask_ephem` | off | Mask known transits before detrending |

### Output

| Flag | Default | Description |
|------|---------|-------------|
| `--outdir DIR` | `.` | Output directory |
| `--suffix TEXT` | none | Suffix appended to output filenames |
| `--survey NAME` | `dss1` | Archival image survey for overlay |
| `-save` | off | Save figure (.png) and TLS results (.h5) |
| `-verbose` | off | Print detailed progress |
| `-overwrite` | off | Overwrite existing output files |

### Multi-sector mode

| Flag | Default | Description |
|------|---------|-------------|
| `--each-sector` | off | Run on every available sector individually |
| `-j, --jobs N` | `1` | Number of parallel jobs for `--each-sector` |

### Examples

```bash
# Basic run on latest sector
ql --name TOI-1234 -save -verbose

# Specific sector and pipeline
ql --name TOI-125.01 --sector 2 --pipeline qlp

# Custom detrending
ql --name HD209458 --flatten_method cosine --window_length 0.3

# Restrict period search range
ql --name TOI-125.01 --period_limits 1 5

# Run all sectors, save results
ql --name TOI-125.01 --each-sector -save

# Run all sectors with 4 parallel workers
ql --name TOI-125.01 --each-sector -j 4 -save

# Pipe output to a log file
ql --name TIC-52368076 -verbose -save | tee output.log
```

## `read_tls` -- Read TLS results

Scan a directory of `.h5` TLS output files and produce a summary CSV:

```bash
read_tls OUTPUT_DIR
```

This creates `OUTPUT_DIR_tls.csv` with columns for period, depth, SDE, and other TLS metrics. Useful for ranking candidates after batch processing.

## `rank_tls` -- Rank candidates

Filter and rank candidates from the CSV produced by `read_tls`:

```bash
rank_tls OUTPUT_DIR --output_dir ranked
```

This copies the plots of candidates that pass the filters into the `ranked/` directory, sorted by Signal Detection Efficiency (SDE).
