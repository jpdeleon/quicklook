"""Unified CLI for TESS quick-look analysis using Typer.

All heavy imports (numpy, pandas, matplotlib, lightkurve, etc.) are deferred
inside the functions that need them so that ``--help`` on any subcommand or
``read-tls`` / ``rank-tls`` execution is fast.
"""

import os
from typing import Optional, Tuple

import typer

from quicklook.exceptions import InvalidInputError


# --- validation helpers (stdlib / typer only) ---

_FLUXTYPE_CHOICES = {"pdcsap", "sap", "aperture", "psf", "auto"}
_PIPELINE_CHOICES = {"spoc", "tess-spoc", "tasoc", "cdips", "pathos", "qlp", "tglc", "t16"}
_BITMASK_CHOICES = {"none", "default", "hard", "hardest"}
_FLATTEN_CHOICES = {
    "biweight",
    "gp",
    "medfilt",
    "median",
    "rspline",
    "hspline",
    "pspline",
    "lowess",
    "cofiam",
    "supersmoother",
    "savgol",
    "notch",
    "locor",
}
_GP_KERNEL_CHOICES = {"periodic_auto", "periodic", "squared_exp"}
_PG_CHOICES = {"gls", "ls", "bls"}


def _validate_choice(value: str, choices: set[str], name: str) -> str:
    low = value.lower()
    if low not in choices:
        valid = ", ".join(sorted(choices))
        raise typer.BadParameter(f"{name} must be one of: {valid}")
    return low


# --- batch helper (heavy imports inside) ---


def _run_ql_for_sector(
    name: str,
    sector: int,
    pipeline: str,
    fluxtype: str,
    quality_bitmask: str,
    flatten_method: str,
    gp_kernel: str,
    pg_method: str,
    edge_cutoff: float,
    outdir: str,
    exptime: Optional[int],
    save: bool,
    overwrite: bool,
    verbose: bool,
    mask_ephem: bool,
    suffix: Optional[str],
    period_limits: Optional[Tuple[float, float]],
    tls_use_threads: int,
    phase_xlim: Optional[float],
    show_simbad: bool = False,
    iterative_search: bool = False,
    min_sde_iterative: float = 5.0,
    max_planets: Optional[int] = None,
):
    import matplotlib.pyplot as pl
    from quicklook.tql import TessQuickLook

    ql = TessQuickLook(
        target_name=name,
        sector=sector,
        pipeline=pipeline,
        exptime=exptime,
        flux_type=fluxtype,
        pg_method=pg_method,
        custom_ephem=None,
        mask_ephem=mask_ephem,
        quality_bitmask=quality_bitmask,
        flatten_method=flatten_method,
        gp_kernel=gp_kernel,
        gp_kernel_size=5,
        window_length=None,
        sigma_clip_raw=(10, 5),
        sigma_clip_flat=None,
        Porb_limits=period_limits,
        edge_cutoff=edge_cutoff,
        phase_xlim=phase_xlim,
        iterative_search=iterative_search,
        min_sde_iterative=min_sde_iterative,
        max_planets=max_planets,
        archival_survey="dss1",
        show_simbad=show_simbad,
        savefig=save,
        savetls=save,
        outdir=outdir,
        verbose=verbose,
        overwrite=overwrite,
        suffix=suffix,
        tls_use_threads=tls_use_threads,
    )
    ql.plot_tql()
    if not save:
        pl.show()
    pl.close()


# --- rank-tls helpers (heavy imports inside) ---


def _is_near_integer(series, tolerance):
    import numpy as np

    fractional = series % 1
    return np.isclose(fractional, 0, atol=tolerance) | np.isclose(fractional, 1, atol=tolerance)


def _apply_rank_filters(
    df,
    tolerance: float = 0.2,
    min_depth: float = 1,
    max_depth: float = 100,
    min_SDE: float = 5,
    max_odd_even_sigma: float = 3.0,
    max_secondary_sde: float = 5.0,
    max_depth_variance_ratio: float = 0.5,
    duration_ratio_bounds: tuple = (0.3, 3.0),
    max_ruwe: float = 1.4,
):
    import numpy as np
    import pandas as pd

    df2 = df.copy()
    depth_range = (df2["depth"] >= min_depth) & (df2["depth"] < max_depth)
    errmsg = f"No candidates satisfy `{min_depth}<depth_range<{max_depth}` ppt."
    assert depth_range.sum() > 0, errmsg

    ratio_forward = _is_near_integer(df2["Prot_gls"] / df2["Porb_tls"], tolerance)
    ratio_inverse = _is_near_integer(df2["Porb_tls"] / df2["Prot_gls"], tolerance)
    not_near = ~(ratio_forward | ratio_inverse)
    errmsg = "No candidates satisfy `not_near_period_ratio`."
    assert not_near.sum() > 0, errmsg

    not_eb = ~df2["simbad_object"].fillna("").str.lower().isin(["eclipsing binary", "eclbin"])
    errmsg = "No candidates satisfy `not_EB`."
    assert not_eb.sum() > 0, errmsg

    strong_signal = df["SDE_tls"] > min_SDE
    errmsg = f"No candidates satisfy `SDE>{min_SDE}`."
    assert strong_signal.sum() > 0, errmsg

    # Advanced vetting metrics (odd-even depth mismatch, secondary eclipse,
    # per-transit depth/duration consistency, Gaia astrometric blend flag).
    # These are computed by TessQuickLook per candidate but were previously
    # dropped before reaching this filter. A metric that is NaN (not
    # computable for a given candidate, or absent from an older CSV) is
    # treated as "not flagged" rather than excluded, matching the
    # rotation-harmonic check above.
    def _numeric_column(name, default=np.nan):
        return pd.to_numeric(df2.get(name, pd.Series(default, index=df2.index)), errors="coerce")

    odd_even = _numeric_column("odd_even_mismatch")
    not_odd_even_eb = ~(odd_even.abs() > max_odd_even_sigma)
    errmsg = f"No candidates satisfy `|odd_even_mismatch|<={max_odd_even_sigma}sigma`."
    assert not_odd_even_eb.sum() > 0, errmsg

    secondary_sde = _numeric_column("secondary_sde")
    no_significant_secondary = ~(secondary_sde > max_secondary_sde)
    errmsg = f"No candidates satisfy `secondary_sde<={max_secondary_sde}`."
    assert no_significant_secondary.sum() > 0, errmsg

    depth_variance_ratio = _numeric_column("depth_variance_ratio")
    consistent_depth = ~(depth_variance_ratio > max_depth_variance_ratio)
    errmsg = f"No candidates satisfy `depth_variance_ratio<={max_depth_variance_ratio}`."
    assert consistent_depth.sum() > 0, errmsg

    duration_ratio = _numeric_column("duration_consistency_ratio")
    lo, hi = duration_ratio_bounds
    consistent_duration = duration_ratio.isna() | duration_ratio.between(lo, hi)
    errmsg = f"No candidates satisfy `duration_consistency_ratio` in {duration_ratio_bounds}."
    assert consistent_duration.sum() > 0, errmsg

    ruwe = _numeric_column("gaia_ruwe")
    not_astrometric_blend = ~(ruwe > max_ruwe)
    errmsg = f"No candidates satisfy `gaia_ruwe<={max_ruwe}`."
    assert not_astrometric_blend.sum() > 0, errmsg

    not_toi = df2["TOI"].isna() | (df2["TOI"] == "")

    exp = pd.to_numeric(
        df2.get("exptime", pd.Series(9999, index=df2.index)), errors="coerce"
    ).fillna(9999)
    count = df2.get("per_transit_count", pd.Series(0, index=df2.index))
    if count.dtype == object:
        import ast

        def _min_count(val):
            if pd.isna(val) or val == "" or val == "0":
                return 0
            try:
                parsed = ast.literal_eval(str(val))
                if isinstance(parsed, (tuple, list)) and len(parsed) > 0:
                    return min(int(x) for x in parsed)
                return int(parsed)
            except (ValueError, SyntaxError):
                return 0

        count = count.apply(_min_count)

    short_valid = (exp <= 120) & (count > 10)
    long_valid = (exp > 120) & (count > 5)
    not_sparse = short_valid | long_valid | (count == 0)

    return df2[
        depth_range
        & not_near
        & not_eb
        & strong_signal
        & not_odd_even_eb
        & no_significant_secondary
        & consistent_depth
        & consistent_duration
        & not_astrometric_blend
        & not_toi
        & not_sparse
    ]


# --- read-tls helper (heavy imports inside) ---


def _summarize_tls_to_csv(input_dir: str, param: str = "SDE_tls") -> str:
    from glob import glob
    from tqdm import tqdm
    from quicklook import h5io
    import pandas as pd

    files = glob(input_dir + "/*.h5")
    if not files:
        raise typer.BadParameter(f"No *.h5 files found in {input_dir}")

    records = []
    errors = []
    for file in tqdm(files, desc="Reading *.h5 files"):
        try:
            data = h5io.load(file)
            d = {
                "TOI": data.get("toiid"),
                "TIC": data.get("ticid"),
                "gaiaid": data.get("gaiaid"),
                "power_gls": data.get("power_gls")[0] if data.get("power_gls") else None,
                "SDE_tls": data.get("SDE"),
                "Porb_tls": data.get("period"),
                "Prot_gls": data.get("Prot_gls")[0] if data.get("Prot_gls") else None,
                "amp_gls": data.get("amp_gls")[0] if data.get("amp_gls") else None,
                "depth": (1 - data.get("depth")) * 1e3 if data.get("depth") is not None else None,
                "per_transit_count": tuple(int(x) for x in data.get("per_transit_count")),
                "n_planets": data.get("n_planets", 1),
                "exptime": data.get("exptime"),
                "pipeline": data.get("pipeline"),
                "flux_type": data.get("flux_type"),
                "simbad_object": data.get("simbad_obj"),
                "odd_even_mismatch": data.get("odd_even_mismatch"),
                "secondary_depth": data.get("secondary_depth"),
                "secondary_sde": data.get("secondary_sde"),
                "depth_variance_ratio": data.get("depth_variance_ratio"),
                "duration_consistency_ratio": data.get("duration_consistency_ratio"),
                "gaia_ruwe": data.get("gaia_ruwe"),
                "filename": file,
            }
            records.append(pd.Series(d))
        except Exception:
            errors.append(data.get("ticid") if data else file)

    df = (
        pd.DataFrame(records)
        .sort_values(by=param, ascending=False)
        .reset_index(drop=True)
        .drop_duplicates()
    )
    fp = input_dir.rstrip("/") + "_tls.csv"
    df.to_csv(fp, index=False)
    if errors:
        typer.echo(f"Errors reading files: {errors}", err=True)
    return fp


# --- commands ---

app = typer.Typer(
    name="ql",
    help="Quick look analysis of TESS lightcurves.",
    no_args_is_help=True,
)


@app.command()
def run(
    name: str = typer.Option(..., "--name", help="Target name"),
    sector: int = typer.Option(-1, "--sector", help="TESS sector (default=-1 for last available)"),
    each_sector: bool = typer.Option(
        False, "--each-sector", help="Run on each available sector individually"
    ),
    each_pipeline: bool = typer.Option(
        False, "--each-pipeline", help="Run on each available pipeline for the latest sector"
    ),
    jobs: int = typer.Option(1, "--jobs", "-j", help="Number of parallel jobs"),
    nice: Optional[int] = typer.Option(
        None, "--nice", help="Lower CPU priority (POSIX; e.g. 19 = lowest)"
    ),
    cores: Optional[int] = typer.Option(None, "--cores", help="CPU cores used by TLS per run"),
    fluxtype: str = typer.Option(
        "pdcsap", "--fluxtype", help=f"Flux type ({'/'.join(_FLUXTYPE_CHOICES)})"
    ),
    pipeline: str = typer.Option(
        "SPOC", "--pipeline", help=f"Pipeline ({'/'.join(_PIPELINE_CHOICES)})"
    ),
    exptime: Optional[int] = typer.Option(None, "--exptime", help="Exposure time (minutes)"),
    quality_bitmask: str = typer.Option(
        "default", "--quality-bitmask", help=f"Quality bitmask ({'/'.join(_BITMASK_CHOICES)})"
    ),
    flatten_method: str = typer.Option(
        "biweight",
        "--flatten-method",
        help=f"Flatten method ({'/'.join(sorted(_FLATTEN_CHOICES))})",
    ),
    gp_kernel: str = typer.Option(
        "periodic_auto", "--gp-kernel", help=f"Wotan GP kernel ({'/'.join(_GP_KERNEL_CHOICES)})"
    ),
    gp_kernel_size: int = typer.Option(5, "--gp-kernel-size", help="Wotan GP kernel size"),
    pg_method: str = typer.Option(
        "gls", "--pg-method", help=f"Periodogram method ({'/'.join(_PG_CHOICES)})"
    ),
    window_length: Optional[float] = typer.Option(
        None, "--window-length", help="Flatten window length (days); 0 = auto-optimize"
    ),
    edge_cutoff: float = typer.Option(0.1, "--edge-cutoff", help="Edge cutoff (days)"),
    sigma_clip_raw: Optional[Tuple[float, float]] = typer.Option(
        (10.0, 5.0), "--sigma-clip-raw", help="(sigma_lo, sigma_hi) for raw LC outlier rejection"
    ),
    sigma_clip_flat: Optional[Tuple[float, float]] = typer.Option(
        None, "--sigma-clip-flat", help="(sigma_lo, sigma_hi) for flattened LC outlier rejection"
    ),
    period_limits: Optional[Tuple[float, float]] = typer.Option(
        None, "--period-limits", help="Period limits for TLS search (low high) in days"
    ),
    phase_xlim: Optional[float] = typer.Option(
        None, "--phase-xlim", help="Phase half-width for transit/eclipse panels"
    ),
    use_priors: bool = typer.Option(
        False, "--use-priors", help="Use ExoFOP stellar params as TLS priors"
    ),
    survey: str = typer.Option("dss1", "--survey", help="Archival image survey name"),
    custom_ephem: Optional[Tuple[float, float, float, float, float, float]] = typer.Option(
        None, "--custom-ephem", help="Custom ephemeris: Tc Tcerr P Perr Tdur Tdurerr"
    ),
    outdir: str = typer.Option(".", "--outdir", help="Output directory"),
    save: bool = typer.Option(False, "--save", help="Save figure and TLS"),
    verbose: bool = typer.Option(False, "--verbose", help="Show details"),
    overwrite: bool = typer.Option(False, "--overwrite", help="Overwrite existing files"),
    show_simbad: bool = typer.Option(
        False, "--show-simbad", help="Overplot nearby SIMBAD objects on archival image"
    ),
    mask_ephem: bool = typer.Option(False, "--mask-ephem", help="Mask transits using ephemerides"),
    suffix: Optional[str] = typer.Option(None, "--suffix", help="Suffix for output filenames"),
    iterative: bool = typer.Option(
        False,
        "--iterative",
        help="Iteratively mask detected transits and re-run TLS to find additional planets",
    ),
    min_sde_iterative: float = typer.Option(
        5.0,
        "--min-sde-iterative",
        help="Minimum SDE for additional iterative candidates (default 5)",
    ),
    max_planets: Optional[int] = typer.Option(
        None,
        "--max-planets",
        help="Maximum number of planets to search for (default: until SDE < --min-sde-iterative)",
    ),
):
    """Run a quick-look analysis of a TESS lightcurve.

    Examples:

        ql run --name TOI-1234 --sector 27 --save --verbose
    """
    from concurrent.futures import ProcessPoolExecutor, as_completed
    from quicklook.cli.ql import sanitize_target_name
    from quicklook.utils import get_available_pipelines, get_available_sectors
    from quicklook.plot import dss_description
    import matplotlib.pyplot as pl
    from quicklook.tql import TessQuickLook

    if nice is not None:
        try:
            new_nice = os.nice(nice)
            if verbose:
                typer.echo(f"Process niceness set to {new_nice}")
        except (AttributeError, OSError) as e:
            typer.echo(f"Warning: could not apply --nice {nice}: {e}", err=True)

    try:
        target_name = sanitize_target_name(name)
    except InvalidInputError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(2)

    _validate_choice(fluxtype, _FLUXTYPE_CHOICES, "fluxtype")
    pipeline_lower = _validate_choice(pipeline, _PIPELINE_CHOICES, "pipeline")
    _validate_choice(quality_bitmask, _BITMASK_CHOICES, "quality_bitmask")
    _validate_choice(flatten_method, _FLATTEN_CHOICES, "flatten_method")
    _validate_choice(gp_kernel, _GP_KERNEL_CHOICES, "gp_kernel")
    _validate_choice(pg_method, _PG_CHOICES, "pg_method")

    if survey not in dss_description:
        valid = ", ".join(dss_description.keys())
        raise typer.BadParameter(f"survey must be one of: {valid}")

    if each_sector and each_pipeline:
        typer.echo("Error: --each-sector and --each-pipeline are mutually exclusive.", err=True)
        raise typer.Exit(2)

    if not save:
        import sys

        if not os.environ.get("DISPLAY") and sys.platform == "linux":
            typer.echo(
                "No display detected on this headless system.\n"
                "Use --save to save the figure and results to disk.",
                err=True,
            )
            raise typer.Exit(1)

    cpu_total = os.cpu_count() or 1
    if each_sector:
        tls_use_threads = cores if cores is not None else max(1, cpu_total // max(1, jobs))
    else:
        tls_use_threads = cores if cores is not None else max(1, cpu_total // 2)
    tls_use_threads = max(1, min(tls_use_threads, cpu_total))

    if each_sector and jobs * tls_use_threads > cpu_total:
        typer.echo(
            f"Warning: jobs ({jobs}) * cores ({tls_use_threads}) = "
            f"{jobs * tls_use_threads} exceeds cpu_count ({cpu_total}); "
            "expect oversubscription.",
            err=True,
        )

    if each_sector:
        sectors = get_available_sectors(target_name, pipeline=pipeline_lower)
        if not sectors:
            typer.echo(f"No sectors found for {target_name} with pipeline {pipeline_lower}")
            raise typer.Exit(1)

        typer.echo(f"Found {len(sectors)} sectors for {target_name}: {sectors}")
        typer.echo(f"Running with {jobs} parallel job(s), {tls_use_threads} TLS core(s) each...")

        failed = 0
        with ProcessPoolExecutor(max_workers=jobs) as executor:
            futures = {
                executor.submit(
                    _run_ql_for_sector,
                    target_name,
                    sector,
                    pipeline_lower,
                    fluxtype,
                    quality_bitmask,
                    flatten_method,
                    gp_kernel,
                    pg_method,
                    edge_cutoff,
                    outdir,
                    exptime,
                    save,
                    overwrite,
                    verbose,
                    mask_ephem,
                    suffix,
                    period_limits,
                    tls_use_threads,
                    phase_xlim,
                    show_simbad,
                    iterative,
                    min_sde_iterative,
                    max_planets,
                ): sector
                for sector in sectors
            }
            for i, future in enumerate(as_completed(futures), 1):
                s = futures[future]
                try:
                    future.result()
                    typer.echo(f"[{i}/{len(sectors)}] sector {s}: done")
                except Exception as e:
                    typer.echo(f"[{i}/{len(sectors)}] sector {s}: FAILED - {e}")
                    failed += 1

        typer.echo(f"\nCompleted: {len(sectors) - failed}/{len(sectors)} sectors successful")
        if failed:
            raise typer.Exit(1)
        raise typer.Exit()

    if each_pipeline:
        pipelines = get_available_pipelines(target_name)
        if not pipelines:
            typer.echo(f"No pipelines found for {target_name}")
            raise typer.Exit(1)

        typer.echo(f"Found {len(pipelines)} pipelines for {target_name}: {pipelines}")
        typer.echo(
            f"Running with {jobs} parallel job(s), {tls_use_threads} TLS core(s) each, "
            "sector=-1 (latest available per pipeline)..."
        )

        failed = 0
        with ProcessPoolExecutor(max_workers=jobs) as executor:
            futures = {
                executor.submit(
                    _run_ql_for_sector,
                    target_name,
                    -1,
                    p,
                    fluxtype,
                    quality_bitmask,
                    flatten_method,
                    gp_kernel,
                    pg_method,
                    edge_cutoff,
                    outdir,
                    exptime,
                    save,
                    overwrite,
                    verbose,
                    mask_ephem,
                    suffix,
                    period_limits,
                    tls_use_threads,
                    phase_xlim,
                    show_simbad,
                    iterative,
                    min_sde_iterative,
                    max_planets,
                ): p
                for p in pipelines
            }
            for i, future in enumerate(as_completed(futures), 1):
                p = futures[future]
                try:
                    future.result()
                    typer.echo(f"[{i}/{len(pipelines)}] pipeline {p}: done")
                except Exception as e:
                    typer.echo(f"[{i}/{len(pipelines)}] pipeline {p}: FAILED - {e}")
                    failed += 1

        typer.echo(f"\nCompleted: {len(pipelines) - failed}/{len(pipelines)} pipelines successful")
        if failed:
            raise typer.Exit(1)
        raise typer.Exit()

    ql = TessQuickLook(
        target_name=target_name,
        sector=sector,
        pipeline=pipeline_lower,
        exptime=exptime,
        flux_type=fluxtype,
        pg_method=pg_method,
        custom_ephem=custom_ephem,
        mask_ephem=mask_ephem,
        quality_bitmask=quality_bitmask,
        flatten_method=flatten_method,
        gp_kernel=gp_kernel,
        gp_kernel_size=gp_kernel_size,
        window_length=window_length,
        sigma_clip_raw=sigma_clip_raw,
        sigma_clip_flat=sigma_clip_flat,
        Porb_limits=period_limits,
        use_star_priors=use_priors,
        phase_xlim=phase_xlim,
        edge_cutoff=edge_cutoff,
        archival_survey=survey,
        show_simbad=show_simbad,
        savefig=save,
        savetls=save,
        outdir=outdir,
        verbose=verbose,
        overwrite=overwrite,
        suffix=suffix,
        tls_use_threads=tls_use_threads,
        iterative_search=iterative,
        min_sde_iterative=min_sde_iterative,
        max_planets=max_planets,
    )
    _ = ql.plot_tql()
    if not save:
        pl.show()
    pl.close()


@app.command(name="read-tls")
def read_tls(
    input_dir: str = typer.Argument(..., help="Input directory containing *.h5 files"),
    param: str = typer.Option("SDE_tls", "--param", "-p", help="Parameter to sort by"),
):
    """Summarize TLS results from *.h5 files into a single CSV file."""
    csv_path = _summarize_tls_to_csv(input_dir, param)
    typer.echo(f"Saved: {csv_path}")


@app.command(name="rank-tls")
def rank_tls(
    input_dir: str = typer.Argument(..., help="Directory containing the original files"),
    csv_path: Optional[str] = typer.Option(None, "--csv-path", help="Path to the CSV file"),
    output_dir: Optional[str] = typer.Option(
        None, "--output-dir", help="Output directory (default: {input_dir}/temp)"
    ),
    column: Optional[str] = typer.Option(
        None, "--column", help="Column to rank by (default: SDE_tls)"
    ),
    remove_filter: bool = typer.Option(
        False, "--remove-filter", help="Remove default filters (rank all)"
    ),
    use_ascending: bool = typer.Option(False, "--use-ascending", help="Ascending sort order"),
    min_sde: float = typer.Option(5.0, "--min-sde", help="Minimum SDE threshold"),
    min_depth: float = typer.Option(1.0, "--min-depth", help="Minimum transit depth (ppt)"),
    max_depth: float = typer.Option(100.0, "--max-depth", help="Maximum transit depth (ppt)"),
    max_odd_even_sigma: float = typer.Option(
        3.0, "--max-odd-even-sigma", help="Max odd/even transit-depth mismatch (sigma)"
    ),
    max_secondary_sde: float = typer.Option(
        5.0, "--max-secondary-sde", help="Max significance of a secondary-eclipse-like dip"
    ),
    max_depth_variance_ratio: float = typer.Option(
        0.5, "--max-depth-variance-ratio", help="Max relative scatter across per-transit depths"
    ),
    duration_ratio_bounds: Tuple[float, float] = typer.Option(
        (0.3, 3.0),
        "--duration-ratio-bounds",
        help="Allowed (min, max) observed/expected transit duration ratio",
    ),
    max_ruwe: float = typer.Option(
        1.4, "--max-ruwe", help="Max Gaia DR3 RUWE (flags likely unresolved binaries/blends)"
    ),
):
    """Copy and prefix PNG files sorted by SDE_tls or a custom column.

    Uses the CSV produced by ``read-tls``. The CSV is auto-generated from
    ``input_dir`` if not provided via ``--csv-path``.
    """
    import pandas as pd
    import shutil
    from pathlib import Path

    if csv_path is None:
        csv_path = f"{input_dir.rstrip('/')}_tls.csv"
        if not os.path.exists(csv_path):
            _summarize_tls_to_csv(input_dir)

    if output_dir is None:
        output_dir = f"{input_dir.rstrip('/')}/temp"
    else:
        output_dir = f"{input_dir.rstrip('/')}/{output_dir}"

    Path(output_dir).mkdir(exist_ok=True)

    df = pd.read_csv(csv_path)
    if column is not None:
        if column not in df.columns:
            raise typer.BadParameter(f"Column '{column}' not found in CSV.")
        df = df.sort_values(by=column, ascending=use_ascending)

    if remove_filter:
        typer.echo("Ranking all files...")
        df2 = df.copy()
    else:
        typer.echo("Applying filters to files...")
        df2 = _apply_rank_filters(
            df,
            tolerance=0.2,
            min_depth=min_depth,
            max_depth=max_depth,
            min_SDE=min_sde,
            max_odd_even_sigma=max_odd_even_sigma,
            max_secondary_sde=max_secondary_sde,
            max_depth_variance_ratio=max_depth_variance_ratio,
            duration_ratio_bounds=duration_ratio_bounds,
            max_ruwe=max_ruwe,
        )
        if len(df2) == 0:
            typer.echo(
                "No good candidates left after applying filters. Try --remove-filter.", err=True
            )
            raise typer.Exit(1)

    filenames = df2.filename.apply(lambda x: x.split("/")[-1].split("_tls.")[0] + ".png").values

    num_digits = len(str(len(filenames)))
    for idx, filename in enumerate(filenames, start=1):
        src_path = Path(input_dir, filename)
        if not src_path.exists():
            typer.echo(f"Warning: File not found: {src_path}", err=True)
            continue
        prefix = str(idx).zfill(num_digits)
        new_name = f"{prefix}_{filename}"
        dst_path = Path(output_dir, new_name)
        shutil.copy2(src_path, dst_path)
        typer.echo(f"Copied: {src_path} -> {dst_path}")


@app.command()
def gui(
    host: str = typer.Option("127.0.0.1", "--host", help="Interface to bind (default: localhost)"),
    port: int = typer.Option(5000, "--port", help="Port to listen on"),
    debug: bool = typer.Option(
        False, "--debug", help="Enable the Werkzeug debugger and auto-reloader (dev only)"
    ),
):
    """Launch the QuickLook web GUI (Flask). Needs the optional gui extra.

    Open http://<host>:<port> in a browser to run analyses interactively.
    The Werkzeug debugger exposes an interactive console, so ``--debug`` is
    opt-in; leave it off on any host other users can reach.

    Examples:

        quicklook gui
        quicklook gui --host 0.0.0.0 --port 8080
    """
    from quicklook.app.app import run_gui

    # An explicit --debug forces the debugger on; without it, defer to the
    # QUICKLOOK_DEBUG environment variable handled inside run_gui.
    run_gui(host=host, port=port, debug=True if debug else None)


if __name__ == "__main__":
    app()
