# Getting Started

## Quick example (CLI)

Run QuickLook on a known exoplanet host:

```bash
ql --name WASP-21 -save -verbose
```

This will:

1. Resolve "WASP-21" to a TIC ID via ExoFOP
2. Download the most recent TESS sector light curve from MAST
3. Run the full pipeline (detrend, GLS periodogram, TLS transit search)
4. Display the 9-panel diagnostic figure
5. Save the figure (`.png`) and TLS results (`.h5`) to the current directory

## Quick example (Python)

```python
from quicklook import TessQuickLook

ql = TessQuickLook(
    target_name="WASP-21",
    sector=56,
    pipeline="SPOC",
    flux_type="pdcsap",
    verbose=True,
)

# Generate the diagnostic figure
fig = ql.plot_tql()

# Access results
print(f"Rotation period: {ql.Prot_ls:.2f} days")
print(f"TLS period: {ql.tls_results.period:.4f} days")
print(f"TLS SDE: {ql.tls_results.SDE:.1f}")
print(f"Transit depth: {ql.tls_results.depth:.6f}")
```

## Quick example (Web GUI)

```bash
ql-gui
```

Open [http://127.0.0.1:5000](http://127.0.0.1:5000) in your browser. Enter a target name, select parameters, and click **Run QuickLook**. Progress is streamed live via WebSocket.

## Choosing a target

QuickLook accepts several target name formats:

| Format | Example |
|--------|---------|
| TOI number | `TOI-1234` or `TOI 1234` |
| TIC ID | `TIC 12345678` |
| Common name | `WASP-21`, `K2-100`, `HD 209458` |

The target name is resolved to a TIC ID using [ExoFOP-TESS](https://exofop.ipac.caltech.edu/tess/).

## Choosing a sector

- `--sector -1` (default): uses the most recent available sector
- `--sector 56`: uses a specific sector
- `--each-sector`: runs the pipeline on every available sector individually

## Choosing a pipeline

TESS light curves are produced by different pipelines:

| Pipeline | Description |
|----------|-------------|
| `spoc` | Science Processing Operations Center (default, 2-min cadence) |
| `tess-spoc` | TESS-SPOC FFI light curves |
| `qlp` | Quick-Look Pipeline (30-min and 10-min cadence) |
| `cdips` | Cluster Difference Imaging Photometric Survey |
| `pathos` | PATHOS pipeline |
| `tglc` | TESS-Gaia Light Curve |
| `tasoc` | TESS Asteroseismic Science Operations Center |

## Next steps

- [CLI Reference](cli.md) -- all command-line options
- [Web GUI](gui.md) -- interactive web interface
- [Pipeline Overview](pipeline.md) -- what each panel means
- [Batch Processing](batch.md) -- analyzing many targets at scale
