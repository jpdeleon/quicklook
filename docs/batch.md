# Batch Processing

QuickLook supports processing many targets at scale using the CLI, shell scripts, or the web GUI.

## CLI: `--each-sector`

Run the pipeline on every available sector for a single target:

```bash
ql run --name TOI-125.01 --each-sector --save
```

Add `-j N` to parallelize across N workers:

```bash
ql run --name TOI-125.01 --each-sector -j 4 --save
```

Every other `ql run` option (`--custom-ephem`, `--window-length`,
`--sigma-clip-raw/--sigma-clip-flat`, `--gp-kernel-size`, `--use-priors`,
`--survey`, ...) applies to each sector exactly as in a single run; the same
holds for `--each-pipeline`. Only the TLS cores per run differ: without
`--cores`, the CPUs are split across the `-j` workers.

## Shell: batch script from a target list

Given a file `tic_ids.txt` with one TIC ID per line:

```
12345678
23456789
34567890
```

### Generate a batch script

```bash
cat tic_ids.txt | while read tic; do
  echo "ql run --name TIC$tic --save --outdir results | tee TIC$tic.log"
done > run_batch.sh
```

### Test a single line

```bash
sed -n '1p' run_batch.sh | sh
```

### Run all targets in parallel

Using [GNU parallel](https://www.gnu.org/software/parallel/):

```bash
cat run_batch.sh | parallel -j 4
```

## Post-processing: rank candidates

After batch processing, use `ql read-tls` and `ql rank-tls` to find the best candidates.

### Step 1: Extract TLS results

```bash
ql read-tls results/
```

This scans all `*_tls.h5` files in `results/` and creates `results_tls.csv` with one row per target.

### Step 2: Rank and filter

```bash
ql rank-tls results/ --output-dir ranked
```

This filters candidates by SDE and other metrics, then copies the plots of the top candidates into `ranked/`.

### Optional: reject shared single events

With `--single-event-search`, scattered light near orbit gaps and momentum
dumps can produce a "single transit" on many stars at once. Compare event
times across the batch:

```bash
ql shared-epochs results/ [--tol 0.25] [--min-targets 3]
```

An event is *shared* when at least `--min-targets` distinct targets (itself
included) have an event within `--tol` days of it. The command writes
`results_shared_events.csv` (one row per event, with `n_targets_nearby` and
`shared`) and `results_shared_targets.csv`. In the targets table,
`event_class_vetted` is `systematic` for mono/duo targets whose events are all
shared. In a 263-star young-host run, it re-labelled 10 sector-14 Kepler-field
monos and duos as systematic and kept every event that coincides with a
catalogue transit.

### Step 3: Visual inspection

Combine ranked plots into a single PDF for quick review:

```bash
pip install img2pdf
img2pdf ranked/*.png --output ranked.pdf
```

## Web GUI: batch submission

The web GUI supports batch submission:

1. Enter multiple target names (one per line) in the batch input field
2. Set shared pipeline parameters
3. Click **Submit Batch**

Jobs are queued and executed one at a time. The queue panel shows progress for all jobs. When a job finishes, the GUI automatically starts watching the next one.

## Web GUI: each-sector mode

1. Enter a single target name
2. Check **Run each available sector**
3. Click **Run QuickLook**

The GUI queries available sectors, then submits one job per sector. All jobs share the same pipeline parameters.

## Tips for large runs

- **Save outputs** (`--save`): always save to disk so results persist across sessions.
- **Use an output directory** (`--outdir results/`): keep results organized.
- **Check for existing outputs**: QuickLook skips targets that already have output files unless `--overwrite` is used.
- **Monitor MAST availability**: the MAST archive occasionally has outages. If many targets fail with HTTP errors, wait and retry.
- **Disk space**: each target produces a ~1 MB PNG and a ~1 MB HDF5 file.
