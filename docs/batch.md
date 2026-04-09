# Batch Processing

QuickLook supports processing many targets at scale using the CLI, shell scripts, or the web GUI.

## CLI: `--each-sector`

Run the pipeline on every available sector for a single target:

```bash
ql --name TOI-125.01 --each-sector -save
```

Add `-j N` to parallelize across N workers:

```bash
ql --name TOI-125.01 --each-sector -j 4 -save
```

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
  echo "ql --name TIC$tic -save --outdir results | tee TIC$tic.log"
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

After batch processing, use `read_tls` and `rank_tls` to find the best candidates.

### Step 1: Extract TLS results

```bash
read_tls results/
```

This scans all `*_tls.h5` files in `results/` and creates `results_tls.csv` with one row per target.

### Step 2: Rank and filter

```bash
rank_tls results/ --output_dir ranked
```

This filters candidates by SDE and other metrics, then copies the plots of the top candidates into `ranked/`.

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

- **Save outputs** (`-save`): always save to disk so results persist across sessions.
- **Use an output directory** (`--outdir results/`): keep results organized.
- **Check for existing outputs**: QuickLook skips targets that already have output files unless `-overwrite` is used.
- **Monitor MAST availability**: the MAST archive occasionally has outages. If many targets fail with HTTP errors, wait and retry.
- **Disk space**: each target produces a ~1 MB PNG and a ~1 MB HDF5 file.
