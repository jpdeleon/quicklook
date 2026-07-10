#!/usr/bin/env python
import os
import sys
import argparse
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# import logging
import matplotlib.pyplot as pl
from quicklook.tql import TessQuickLook
from quicklook.utils import get_available_pipelines, get_available_sectors
from quicklook.plot import dss_description


def run_ql_for_sector(
    name,
    sector,
    pipeline,
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
    show_simbad=False,
):
    """Run ql for a single sector. Used by --each-sector."""
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


long_decription = """Run a quick look analysis of a TESS lightcurve.
Notes:
* use single hyphen (-flag) if no value is needed.
* use double hyphen (--flag value) if value is needed.

Example: ql --name TOI-1234 --sector 27 -save -verbose
"""


def sanitize_target_name(name: str) -> str:
    """
    Sanitize the user input target name for TessQuickLook.

    Rules:
    - Strip leading/trailing spaces
    - If the target starts with 'TOI', ensure it is formatted as 'TOI-xxxx'
    - Otherwise, remove spaces
    """
    name = name.strip()
    if name[:3].lower() == "toi":
        # Replace spaces with hyphen
        name = name.replace(" ", "-")
        # Ensure single hyphen after TOI
        if not name.upper().startswith("TOI-"):
            name = "TOI-" + name[3:].lstrip("-")
    else:
        # Remove all spaces for non-TOI targets
        name = name.replace(" ", "")
    return name


def main():
    parser = argparse.ArgumentParser(
        description=long_decription, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--name", type=str, help="target name", default=None)
    parser.add_argument(
        "--sector",
        # type=int,
        help="TESS sector (default=-1 (last available sector))",
        default=-1,
    )
    parser.add_argument(
        "--each-sector",
        action="store_true",
        help="run on each available sector for the given pipeline individually",
        default=False,
    )
    parser.add_argument(
        "--each-pipeline",
        action="store_true",
        help="run on each available pipeline for the latest sector "
        "(mutually exclusive with --each-sector)",
        default=False,
    )
    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=1,
        help="number of parallel jobs (default=1)",
    )
    parser.add_argument(
        "--nice",
        type=int,
        default=None,
        help="lower CPU priority by this increment (POSIX; e.g. 19 = lowest). Default: unchanged.",
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=None,
        help="CPU cores used by TLS per run. "
        "Default: cpu_count//2 for single run; cpu_count//jobs for --each-sector.",
    )
    # parser.add_argument(
    #     "-c",
    #     "--cadence",
    #     type=str,
    #     choices=["long", "short"],
    #     help="30-min long or 2-min short (default)",
    #     default="short",
    # )
    # parser.add_argument(
    #     "-sr",
    #     "--search_radius",
    #     type=float,
    #     help="search radius in arcsec (default=3)",
    #     default=3,
    # )
    parser.add_argument(
        "--fluxtype",
        type=str,
        help="type of lightcurve (SPOC: pdcsap/sap; TGLC: aperture/psf/auto)",
        choices=["pdcsap", "sap", "aperture", "psf", "auto"],
        default="pdcsap",
    )
    parser.add_argument(
        "--pipeline",
        type=str,
        help="lightcurve produced from which pipeline (default=SPOC)",
        choices=["spoc", "tess-spoc", "tasoc", "cdips", "pathos", "qlp", "tglc", "t16"],
        default="SPOC",
    )
    parser.add_argument(
        "--exptime",
        type=int,
        help="exposure time (default is whatever is used in available sector)",
        default=None,
    )
    # parser.add_argument(
    #     "--binsize",
    #     type=float,
    #     help="size (in minutes) to bin the raw lc",
    #     default=None,
    # )
    # parser.add_argument(
    #     "-m",
    #     "--mission",
    #     type=str,
    #     help="TESS or K2 or Kepler",
    #     default='TESS',
    # )
    # parser.add_argument(
    #     "-a",
    #     "--aper_mask",
    #     type=str,
    #     help="aperture mask type",
    #     choices=["pipeline", "round", "square", "percentile", "threshold"],
    #     default=None,
    # )
    # parser.add_argument(
    #     "-t", "--threshold", type=float, help="mask threshold in sigma", default=5
    # )
    # parser.add_argument(
    #     "-r", "--aper_radius", type=int, help="mask radius in pix", default=1
    # )
    # parser.add_argument(
    #     "-perc", "--percentile", type=float, help="mask percentile", default=90
    # )
    parser.add_argument(
        "--quality_bitmask",
        help="remove specific data points identified in TESS data release notes",
        type=str,
        choices=["none", "default", "hard", "hardest"],
        default="default",
    )
    # parser.add_argument(
    #     "-size",
    #     "--cutout_size",
    #     nargs=2,
    #     type=float,
    #     help="FFI cutout size for long cadence (default=[12,12] pix)",
    #     default=(12, 12),
    # )
    parser.add_argument(
        "--flatten_method",
        type=str,
        help="wotan flatten method (default=biweight); see https://wotan.readthedocs.io/en/latest/Usage.html",
        default="biweight",
        choices=[
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
        ],
    )
    parser.add_argument(
        "--gp_kernel",
        type=str,
        help="wotan gp kernel (default=periodic_auto)",
        default="periodic_auto",
        choices=["periodic_auto", "periodic", "squared_exp"],
    )
    parser.add_argument(
        "--gp_kernel_size",
        type=int,
        help="wotan gp kernel size (default=5)",
        default=5,
    )
    parser.add_argument(
        "--pg_method",
        type=str,
        help="periodogram method (default=gls)",
        default="gls",
        choices=["gls", "ls", "bls"],
    )
    parser.add_argument(
        "--window_length",
        type=float,
        help="flatten method window length (default=0.5 days); window length is optimized when set to 0",
        default=None,
    )
    parser.add_argument(
        "--edge_cutoff",
        type=int,
        help="cut each edges (default=0.1 days)",
        default=0.1,
    )
    parser.add_argument(
        "--sigma_clip_raw",
        type=float,
        help="(sigma_lo,sigma_hi) for outlier rejection of raw lc before flattening/detrending",
        nargs=2,
        default=(10, 5),
    )
    parser.add_argument(
        "--sigma_clip_flat",
        help="(sigma_lo,sigma_hi) for outlier rejection of flattened/detrended lc",
        nargs=2,
        type=float,
        default=None,
    )
    parser.add_argument(
        "--period_limits",
        help="period limits in TLS search; default=(0.1, baseline/2) d",
        nargs=2,
        type=float,
        default=None,
    )
    parser.add_argument(
        "--phase_xlim",
        type=float,
        default=None,
        help=(
            "phase half-width for panels 7 and 8. Example: 0.1 plots "
            "transit phase -0.1..0.1 and eclipse phase 0.4..0.6. "
            "Default: automatic zoom from transit duration."
        ),
    )
    parser.add_argument(
        "-u",
        "--use_priors",
        action="store_true",
        help="use ExoFOP stellar params (R_star, M_star) as TLS priors",
        default=False,
    )
    parser.add_argument(
        "--survey",
        help="archival image survey name if using img option (default=dss1)",
        choices=list(dss_description.keys()),
        default="dss1",
    )
    parser.add_argument(
        "--custom_ephem",
        help="custom ephemeris in days. Example: --custom_ephem Tc Tcerr P Perr Tdur Tdurerr",
        nargs=6,
        type=float,
        default=None,
    )
    parser.add_argument("--outdir", type=str, help="output directory", default=".")
    parser.add_argument(
        "-save",
        action="store_true",
        help="save figure and tls",
        default=False,
    )
    parser.add_argument("-verbose", action="store_true", help="show details", default=False)
    parser.add_argument("-overwrite", action="store_true", help="overwrite files", default=False)
    parser.add_argument(
        "-show_simbad",
        action="store_true",
        help="overplot nearby SIMBAD objects on the archival image (default=False)",
        default=False,
    )
    # parser.add_argument(
    #     "-use_tpf_image",
    #     action="store_true",
    #     help="plot gaia sources on tpf image instead of archival image (default=False)",
    #     default=False,
    # )
    parser.add_argument(
        "-mask_ephem",
        help="mask transits either using TOI or custom ephemerides if available (default=False)",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--suffix",
        help="add suffix to filename if -save flag is used",
        type=str,
        default=None,
    )

    # prints help if no argument supplied
    args = parser.parse_args(None if sys.argv[1:] else ["-h"])

    if args.nice is not None:
        try:
            new_nice = os.nice(args.nice)
            if args.verbose:
                print(f"Process niceness set to {new_nice}")
        except (AttributeError, OSError) as e:
            print(
                f"Warning: could not apply --nice {args.nice}: {e}",
                file=sys.stderr,
            )

    target_name = sanitize_target_name(args.name)

    cpu_total = os.cpu_count() or 1
    if args.each_sector:
        tls_use_threads = (
            args.cores if args.cores is not None else max(1, cpu_total // max(1, args.jobs))
        )
    else:
        tls_use_threads = args.cores if args.cores is not None else max(1, cpu_total // 2)
    tls_use_threads = max(1, min(tls_use_threads, cpu_total))
    if args.each_sector and args.each_pipeline:
        print("Error: --each-sector and --each-pipeline are mutually exclusive.", file=sys.stderr)
        sys.exit(2)

    if args.each_sector and args.jobs * tls_use_threads > cpu_total:
        print(
            f"Warning: jobs ({args.jobs}) * cores ({tls_use_threads}) = "
            f"{args.jobs * tls_use_threads} exceeds cpu_count ({cpu_total}); "
            "expect oversubscription.",
            file=sys.stderr,
        )

    if args.each_sector:
        sectors = get_available_sectors(target_name, pipeline=args.pipeline)
        if not sectors:
            print(f"No sectors found for {target_name} with pipeline {args.pipeline}")
            sys.exit(1)

        print(f"Found {len(sectors)} sectors for {target_name}: {sectors}")
        print(f"Running with {args.jobs} parallel job(s), {tls_use_threads} TLS core(s) each...")

        failed = 0
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            futures = {}
            for sector in sectors:
                future = executor.submit(
                    run_ql_for_sector,
                    target_name,
                    sector,
                    args.pipeline,
                    args.fluxtype,
                    args.quality_bitmask,
                    args.flatten_method,
                    args.gp_kernel,
                    args.pg_method,
                    args.edge_cutoff,
                    args.outdir,
                    args.exptime,
                    args.save,
                    args.overwrite,
                    args.verbose,
                    args.mask_ephem,
                    args.suffix,
                    args.period_limits,
                    tls_use_threads,
                    args.phase_xlim,
                    args.show_simbad,
                )
                futures[future] = sector

            for i, future in enumerate(as_completed(futures), 1):
                sector = futures[future]
                try:
                    future.result()
                    print(f"[{i}/{len(sectors)}] sector {sector}: done")
                except Exception as e:
                    print(f"[{i}/{len(sectors)}] sector {sector}: FAILED - {e}")
                    failed += 1

        print(f"\nCompleted: {len(sectors) - failed}/{len(sectors)} sectors successful")
        if failed:
            sys.exit(1)
        sys.exit(0)

    if args.each_pipeline:
        pipelines = get_available_pipelines(target_name)
        if not pipelines:
            print(f"No pipelines found for {target_name}")
            sys.exit(1)

        print(f"Found {len(pipelines)} pipelines for {target_name}: {pipelines}")
        print(
            f"Running with {args.jobs} parallel job(s), {tls_use_threads} TLS core(s) each, "
            "sector=-1 (latest available per pipeline)..."
        )

        failed = 0
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            futures = {}
            for pipeline in pipelines:
                future = executor.submit(
                    run_ql_for_sector,
                    target_name,
                    -1,  # latest sector for each pipeline
                    pipeline,
                    args.fluxtype,
                    args.quality_bitmask,
                    args.flatten_method,
                    args.gp_kernel,
                    args.pg_method,
                    args.edge_cutoff,
                    args.outdir,
                    args.exptime,
                    args.save,
                    args.overwrite,
                    args.verbose,
                    args.mask_ephem,
                    args.suffix,
                    args.period_limits,
                    tls_use_threads,
                    args.phase_xlim,
                    args.show_simbad,
                )
                futures[future] = pipeline

            for i, future in enumerate(as_completed(futures), 1):
                pipeline = futures[future]
                try:
                    future.result()
                    print(f"[{i}/{len(pipelines)}] pipeline {pipeline}: done")
                except Exception as e:
                    print(f"[{i}/{len(pipelines)}] pipeline {pipeline}: FAILED - {e}")
                    failed += 1

        print(f"\nCompleted: {len(pipelines) - failed}/{len(pipelines)} pipelines successful")
        if failed:
            sys.exit(1)
        sys.exit(0)

    ql = TessQuickLook(
        # gaia2id=args.gaia2id,
        # gaia3id=args.gaia3id,
        # toiid=args.toi,
        # ticid=args.tic,
        # coords=args.coords,
        target_name=target_name,
        # search_radius=args.search_radius,
        sector=args.sector,
        # cadence=args.cadence,
        pipeline=args.pipeline,
        exptime=args.exptime,
        flux_type=args.fluxtype,
        pg_method=args.pg_method,
        custom_ephem=args.custom_ephem,
        mask_ephem=args.mask_ephem,
        # sap_mask=args.aper_mask,
        # aper_radius=args.aper_radius,
        # threshold_sigma=args.threshold,
        # percentile=args.percentile,
        quality_bitmask=args.quality_bitmask,
        # apply_data_quality_mask=args.quality_mask,
        flatten_method=args.flatten_method,
        gp_kernel=args.gp_kernel,
        gp_kernel_size=args.gp_kernel_size,
        window_length=args.window_length,
        sigma_clip_raw=args.sigma_clip_raw,
        sigma_clip_flat=args.sigma_clip_flat,
        # cutout_size=args.cutout_size,
        # bin_hr=args.bin_hr,
        Porb_limits=args.period_limits,
        use_star_priors=args.use_priors,
        phase_xlim=args.phase_xlim,
        edge_cutoff=args.edge_cutoff,
        # find_cluster=args.find_cluster,
        # nearby_gaia_radius=args.nearby_gaia_radius,
        # run_gls=args.gls,
        archival_survey=args.survey,
        show_simbad=args.show_simbad,
        savefig=args.save,
        savetls=args.save,
        outdir=args.outdir,
        verbose=args.verbose,
        overwrite=args.overwrite,
        suffix=args.suffix,
        tls_use_threads=tls_use_threads,
        # check_if_variable=args.check_if_variable,
        # estimate_spec_type=args.spec_type,
        # estimate_gyro_age=args.gyro_age,
    )
    _ = ql.plot_tql()
    if not args.save:
        pl.show()
    pl.close()


if __name__ == "__main__":
    main()
