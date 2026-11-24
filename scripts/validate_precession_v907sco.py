#!/usr/bin/env python
"""Real-data check of the EB precession branch on V907 Sco (needs network).

Internal dev tool, not part of the installed package. V907 Sco has only FFI
light curves, at a different cadence in each sector (S13: 30 min, S39: 10 min,
S66: 200 s), so ``ql run --sector all`` (one ``--exptime``) cannot stitch them.
This script downloads the TESS-SPOC light curves, flattens each sector with
both eclipses masked, and calls ``search_eb_precession`` directly.

Usage:

    uv run python scripts/validate_precession_v907sco.py [--flux pdcsap_flux]
        [--third-light 0.0] [--out docs/img/precession-v907sco.png]

Reference: Zasche et al. (2023), arXiv:2302.04613, Table 2 (i1 = 73.37 deg in
S13, 76.57/76.61 deg in the halves of S39).
"""

from __future__ import annotations

import argparse
import os

os.environ.setdefault("MPLBACKEND", "Agg")

import lightkurve as lk
import numpy as np

from quicklook.plot_precession import plot_eb_precession
from quicklook.precession import sector_span_from_lc
from quicklook.precession_eb import search_eb_precession

PERIOD = 3.7762049  # d, Zasche et al. 2023, Table 1
T0_BTJD = 2459362.3272 - 2457000.0  # primary eclipse
DURATION = 5.3 / 24  # d
MASK_DURATIONS = 1.6
FLATTEN_WINDOW_DAYS = 1.5


def _flatten(lc):
    lc = lc.remove_nans().normalize()
    width = MASK_DURATIONS * DURATION
    mask = lc.create_transit_mask(
        period=[PERIOD, PERIOD],
        transit_time=[T0_BTJD, T0_BTJD + PERIOD / 2],
        duration=[width, width],
    )
    cadence = float(np.median(np.diff(lc.time.value)))
    window = max(int(FLATTEN_WINDOW_DAYS / cadence) | 1, 7)  # odd, at least 7 cadences
    return lc, lc.flatten(window_length=window, mask=mask)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--flux", default="pdcsap_flux", help="flux column to download")
    parser.add_argument("--third-light", type=float, default=0.0)
    parser.add_argument("--out", default="precession-v907sco.png")
    args = parser.parse_args()

    search = lk.search_lightcurve("V907 Sco", mission="TESS", author="TESS-SPOC")
    times, fluxes, errs, spans = [], [], [], []
    for row in search:
        raw, flat = _flatten(row.download(flux_column=args.flux, quality_bitmask="default"))
        span = sector_span_from_lc(raw)
        cadence = np.median(np.diff(raw.time.value)) * 1440
        print(f"S{span.sector}: {len(flat)} pts, {cadence:.1f} min, CROWDSAP={span.crowdsap:.4f}")
        times.append(flat.time.value)
        fluxes.append(flat.flux.value)
        errs.append(flat.flux_err.value)
        spans.append(span)

    result = search_eb_precession(
        np.concatenate(times),
        np.concatenate(fluxes),
        np.concatenate(errs),
        spans,
        PERIOD,
        T0_BTJD,
        DURATION,
        third_light=args.third_light,
    )
    b = result.binary
    if b is None:
        raise SystemExit(result.reason)
    print(
        f"binary: r_sum={b.r_sum:.4f} k={b.k:.3f} J={b.sb_ratio:.3f} "
        f"dt_sec={b.dt_sec * 1440:.1f} min"
    )
    for f in result.sector_fits:
        print(
            f"S{f.sector}: i={f.inc:.3f}±{f.inc_err:.3f} deg, T14={f.t14 * 24:.2f} h, "
            f"depth={f.depth_data * 1e3:.1f}±{f.depth_data_err * 1e3:.1f} ppt"
        )
    print(
        f"di/dt = {result.di_dt:.3f} ± {result.di_dt_err:.3f} deg/yr "
        f"({result.significance:.1f} sigma); P_nodal <= {result.max_nodal_period_yr:.0f} yr"
    )
    print(f"candidate={result.is_candidate} flags={result.flags}")
    fig = plot_eb_precession(result, "V907 Sco (TESS-SPOC FFI)")
    fig.savefig(args.out, dpi=120, bbox_inches="tight")
    print(f"Saved: {args.out}")


if __name__ == "__main__":
    main()
