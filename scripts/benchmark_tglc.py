"""Benchmark: does the local TGLC pipeline reproduce the MAST TGLC HLSP?

Example
-------
    uv run scripts/benchmark_tglc.py "TIC 259377017" --sector 3 --size 99

Writes ``<outdir>/<target>_s<sector>_benchmark.json`` and a comparison plot,
and exits non-zero if the calibrated light curves fail the agreement criteria.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys

from quicklook.benchmark_tglc import (
    DEFAULT_MAX_RESID_TO_SCATTER,
    DEFAULT_MIN_OVERLAP_FRAC,
    DEFAULT_MIN_PEARSON,
    plot_comparison,
    run_benchmark,
)

# Non-saturated, uncrowded default targets with MAST TGLC coverage.
DEFAULT_TARGETS = (("TIC 259377017", 3),)


def _parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("target", nargs="?", help="target name, e.g. 'TIC 259377017'")
    p.add_argument("--sector", type=int, help="TESS sector (required with target)")
    p.add_argument("--size", type=int, default=99, help="TESScut size in px (50-99)")
    p.add_argument("--limit-mag", type=float, default=16)
    p.add_argument("--outdir", default="benchmark_tglc_output")
    p.add_argument("--no-cache", action="store_true", help="force fresh FFI download + ePSF fit")
    p.add_argument("--no-plot", action="store_true")
    p.add_argument("--min-pearson", type=float, default=DEFAULT_MIN_PEARSON)
    p.add_argument("--max-resid-to-scatter", type=float, default=DEFAULT_MAX_RESID_TO_SCATTER)
    p.add_argument("--min-overlap-frac", type=float, default=DEFAULT_MIN_OVERLAP_FRAC)
    args = p.parse_args(argv)
    if args.target and args.sector is None:
        p.error("--sector is required when a target is given")
    return args


def _format_table(result) -> str:
    head = f"{'column':<15}{'N':>6}{'r':>8}{'med_ratio':>11}{'resid/scat':>12}{'scat_ratio':>12}"
    rows = [head, "-" * len(head)]
    for col, m in result.columns.items():
        rows.append(
            f"{col:<15}{m['n_compared']:>6}{m['pearson_r']:>8.3f}{m['median_ratio']:>11.3f}"
            f"{m['resid_to_scatter']:>12.3f}{m['scatter_ratio']:>12.3f}"
        )
    return "\n".join(rows)


def main(argv=None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    args = _parse_args(argv)
    targets = [(args.target, args.sector)] if args.target else list(DEFAULT_TARGETS)
    os.makedirs(args.outdir, exist_ok=True)

    all_passed = True
    for target, sector in targets:
        result, local_lc, mast_lc = run_benchmark(
            target,
            sector,
            size=args.size,
            limit_mag=args.limit_mag,
            use_cache=not args.no_cache,
            min_pearson=args.min_pearson,
            max_resid_to_scatter=args.max_resid_to_scatter,
            min_overlap_frac=args.min_overlap_frac,
        )
        stem = os.path.join(args.outdir, f"{target.replace(' ', '')}_s{sector}_benchmark")
        with open(f"{stem}.json", "w") as fh:
            json.dump(result.to_dict(), fh, indent=2, default=float)
        if not args.no_plot:
            plot_comparison(local_lc, mast_lc, result, f"{stem}.png")

        status = "PASS" if result.passed else "FAIL"
        print(f"\n{target} sector {sector}: {status}")
        print(
            f"cadences: local={result.n_local} mast(good)={result.n_mast} "
            f"matched={result.n_matched} ({result.overlap_frac:.0%})"
        )
        print(_format_table(result))
        for reason in result.failures:
            print(f"  - {reason}")
        all_passed &= result.passed

    return 0 if all_passed else 1


if __name__ == "__main__":
    sys.exit(main())
