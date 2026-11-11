"""Helpers for long-period targets.

* ``default_period_ceiling`` - TLS period ceiling. TLS' default stops at half
  the baseline, which hides planets with 2-3 transits (e.g. a 20 d planet in a
  single ~27 d sector). The extended ceiling reaches nearly the full baseline.
* ``epoch_check`` - is there a dip at the catalogue-predicted transit times?
  Independent of the TLS period search, so it still answers "is the signal in
  the data" when TLS picks an alias or another signal.
"""

from __future__ import annotations

import numpy as np

EPOCH_SNR_MIN = 5.0  # stacked SNR for epoch_detected
EXTENDED_CEILING_FRACTION = 0.98
_N_CONTROL = 30
_IN_TRANSIT_FRACTION = 0.4  # half-width of the in-transit window, in durations
_MIN_COVERAGE = 0.5  # fraction of expected in-transit cadences required
_MIN_OOT_POINTS = 10


def default_period_ceiling(time_span: float, extend: bool = True) -> float:
    """TLS ``period_max`` for a light curve of the given baseline (days)."""
    return EXTENDED_CEILING_FRACTION * time_span if extend else time_span / 2


def needs_standard_fallback(n_transits: int, period_max: float, time_span: float) -> bool:
    """True if an extended-range TLS winner is a single-transit period.

    With the extended ceiling TLS may be asked for ``n_transits_min=1``; a
    single event cannot define a period, so such a winner is discarded in
    favour of the standard (<= half baseline) search.
    """
    return period_max > time_span / 2 and n_transits < 2


def _stack(time, flux, period, epoch, duration, shift=0.0):
    """Stacked dip SNR over covered predicted transits -> (snr, n_covered)."""
    half = _IN_TRANSIT_FRACTION * duration
    k0 = int(np.ceil((time[0] - epoch - shift) / period))
    k1 = int(np.floor((time[-1] - epoch - shift) / period))
    inn, out, covered = [], [], 0
    for k in range(k0, k1 + 1):
        tc = epoch + shift + k * period
        lo, hi = np.searchsorted(time, [tc - 3 * duration, tc + 3 * duration])
        tt, ff = time[lo:hi], flux[lo:hi]
        if len(tt) < _MIN_OOT_POINTS:
            continue
        in_win = np.abs(tt - tc) < half
        cadence = np.median(np.diff(tt))
        if in_win.sum() < _MIN_COVERAGE * (2 * half / cadence):
            continue
        oot = np.abs(tt - tc) > 0.75 * duration
        if oot.sum() < _MIN_OOT_POINTS:
            continue
        covered += 1
        inn.append(ff[in_win])
        out.append(ff[oot])
    if covered == 0:
        return np.nan, 0
    inn, out = np.concatenate(inn), np.concatenate(out)
    sigma = 1.4826 * np.median(np.abs(out - np.median(out)))
    if not sigma > 0:
        return np.nan, covered
    return float((np.median(out) - inn.mean()) / (sigma / np.sqrt(len(inn)))), covered


def epoch_check(
    time, flux, period, epoch, duration, n_control: int = _N_CONTROL, seed: int = 0
) -> dict:
    """Stack the flattened flux at the catalogue ephemeris.

    Parameters are in days (``epoch`` on the same time system as ``time``).
    Returns ``epoch_snr``, ``epoch_n_covered``, ``epoch_control_snr_p95`` (the
    95th percentile of the same statistic at random phase offsets, a
    per-target false-alarm yardstick) and ``epoch_detected``.
    """
    empty = {
        "epoch_snr": np.nan,
        "epoch_n_covered": 0,
        "epoch_control_snr_p95": np.nan,
        "epoch_detected": False,
    }
    vals = np.array([period, epoch, duration], dtype=float)
    if not (np.all(np.isfinite(vals)) and period > 0 and duration > 0):
        return empty
    time, flux = np.asarray(time, dtype=float), np.asarray(flux, dtype=float)
    ok = np.isfinite(time) & np.isfinite(flux)
    if ok.sum() < 2 * _MIN_OOT_POINTS:
        return empty
    order = np.argsort(time[ok])
    time, flux = time[ok][order], flux[ok][order]

    snr, covered = _stack(time, flux, period, epoch, duration)
    if covered == 0:
        return empty
    rng = np.random.default_rng(seed)
    control = [
        c
        for s in rng.uniform(0.15 * period, 0.85 * period, n_control)
        for c, n in [_stack(time, flux, period, epoch, duration, shift=s)]
        if n and np.isfinite(c)
    ]
    return {
        "epoch_snr": snr,
        "epoch_n_covered": covered,
        "epoch_control_snr_p95": float(np.percentile(control, 95)) if control else np.nan,
        "epoch_detected": bool(np.isfinite(snr) and snr >= EPOCH_SNR_MIN),
    }
