"""Detrending dedicated to the single-event search.

The pipeline's detrending window is tuned for TLS: short enough to remove
stellar variability, often 3x a catalogue duration or 0.5-1 d. Long single
transits need more room, and even a window 3x the transit duration sags into
the dip near a data edge (TOI-2180 b in TESS S48). So the single-event search
detrends the raw light curve itself:

1. ``choose_window``: at least ``window_floor`` (3x the longest duration the
   search tries for the host star, ~2.9 d for a Sun-like star), unless the user
   set the window or the star rotates fast (a wide window would leave the
   variability in).
2. Search with trial durations capped at window / ``WINDOW_TO_DURATION``.
3. Re-detrend with the detected events masked (trend interpolated across
   them), as in the wotan documentation, and search again.

TLS and the default figure keep using the pipeline's own flattened light curve.
"""

from __future__ import annotations

import numpy as np
from astropy.timeseries import LombScargle
from wotan import flatten

from quicklook.single_event import DEFAULT_P_MAX, MAX_DURATION, _ECC_ALLOWANCE, central_duration

WINDOW_TO_DURATION = 3.0  # Hippke et al. 2019: windows < ~3x duration distort transits
FAST_ROTATION_DAYS = 2.0
ROTATION_POWER_MIN = 0.1  # Lomb-Scargle power for a real rotation signal
MASK_FACTOR = 1.0  # extra masked half-width, in units of half the event duration
DETREND_METHOD = "biweight"  # robust, and well defined with a masked trend
_ROTATION_P_MIN = 0.1  # days
_ROTATION_BIN = 30 / 1440  # days; bin before the periodogram for speed
_MIN_POINTS = 20


def longest_trial_duration(rho_star: float = 1.0) -> float:
    """Longest duration (days) the single-event search tries for this star."""
    return min(MAX_DURATION, _ECC_ALLOWANCE * central_duration(DEFAULT_P_MAX, rho_star))


def window_floor(rho_star: float = 1.0) -> float:
    """Smallest window (days) that preserves every duration the search tries.

    A 1.6 d biweight (3x a 1-yr central duration) removes a 22 h transit
    entirely, so the floor follows the longest trial duration instead.
    """
    return WINDOW_TO_DURATION * longest_trial_duration(rho_star)


def choose_window(
    base_window: float, rho_star: float, user_set: bool, fast_rotator: bool
) -> tuple[float, str]:
    """Window for the single-event detrending and why it was chosen.

    ``reason`` is ``"user"`` (window given explicitly), ``"fast rotator"``
    (kept the pipeline window), ``"pipeline"`` (pipeline window already wide
    enough) or ``"floor"`` (raised to ``window_floor``).
    """
    if user_set:
        return float(base_window), "user"
    if fast_rotator:
        return float(base_window), "fast rotator"
    floor = window_floor(rho_star)
    if base_window >= floor:
        return float(base_window), "pipeline"
    return float(floor), "floor"


def is_fast_rotator(time, flux) -> bool:
    """True if a significant Lomb-Scargle peak lies below ``FAST_ROTATION_DAYS``."""
    t, f = np.asarray(time, dtype=float), np.asarray(flux, dtype=float)
    ok = np.isfinite(t) & np.isfinite(f)
    t, f = t[ok], f[ok]
    if t.size < _MIN_POINTS:
        return False
    edges = np.arange(t.min(), t.max() + _ROTATION_BIN, _ROTATION_BIN)
    idx = np.digitize(t, edges)
    counts = np.bincount(idx)
    keep = counts > 0
    tb = (np.bincount(idx, weights=t)[keep]) / counts[keep]
    fb = (np.bincount(idx, weights=f)[keep]) / counts[keep]
    span = tb.max() - tb.min()
    if span <= 2 * _ROTATION_P_MIN:
        return False
    freq = np.linspace(1 / (span / 2), 1 / _ROTATION_P_MIN, 5000)
    power = LombScargle(tb, fb).power(freq)
    best = int(np.argmax(power))
    return bool(power[best] >= ROTATION_POWER_MIN and 1 / freq[best] < FAST_ROTATION_DAYS)


def event_mask(time, events) -> np.ndarray:
    """Cadences within each event plus a margin of ``MASK_FACTOR`` half-durations."""
    t = np.asarray(time, dtype=float)
    mask = np.zeros(t.size, dtype=bool)
    for ev in events:
        half = ev.duration / 2 * (1 + MASK_FACTOR)
        mask |= np.abs(t - ev.t0) < half
    return mask


def detrend(time, flux, window: float, mask=None, edge_cutoff: float = 0.0):
    """Biweight-detrend ``flux``; masked cadences do not shape the trend.

    Returns ``(flat, trend)`` on the input cadences. The trend under the mask
    is linearly interpolated from the unmasked trend (as in the wotan docs);
    non-finite input stays NaN. Variability on timescales within a few times
    the masked span cannot be separated from a long transit this way (or any
    other): a 22 h, 4 ppt transit on a 2 ppt, P = 6 d signal comes out ~25%
    shallow, while on a P = 20 d signal it is recovered.
    """
    t, f = np.asarray(time, dtype=float), np.asarray(flux, dtype=float)
    use = np.isfinite(t) & np.isfinite(f)
    if mask is not None:
        use &= ~np.asarray(mask, dtype=bool)
    trend = np.full(t.size, np.nan)
    if use.sum() < _MIN_POINTS:
        return np.full(t.size, np.nan), trend
    _, trend_use = flatten(
        t[use],
        f[use],
        method=DETREND_METHOD,
        window_length=window,
        edge_cutoff=edge_cutoff,
        break_tolerance=window / 2,
        return_trend=True,
    )
    good = np.isfinite(trend_use)
    t_ref, y_ref = t[use][good], trend_use[good]
    finite = np.isfinite(t) & np.isfinite(f)
    trend[finite] = np.interp(t[finite], t_ref, y_ref)
    return f / trend, trend


__all__ = ["choose_window", "detrend", "event_mask", "is_fast_rotator", "window_floor"]
