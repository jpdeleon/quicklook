"""Transit-timing variations: per-transit fits, O-C and event chaining.

Step 5 of issue #42. TLS folds on a strictly linear ephemeris, so transits
displaced by more than about half a duration are smeared out, and masking at
the linear ephemeris leaves them in the light curve for the next search. Here
each transit is fitted on its own:

1. ``fit_transit_times``: around each predicted time ``t0 + n*P`` fit a
   trapezoid with free mid-time (within +-``delta``), depth and duration on a
   local linear baseline. Epochs without enough data on the transit and on
   both sides of it are skipped.
2. ``ttv_statistics``: weighted linear ephemeris through the fitted times,
   O-C, and the chi-square p-value of the linear model. A small p-value means
   the times are not consistent with a constant period.
3. ``cycle_count_solutions``: across long gaps (e.g. between sectors years
   apart) the number of orbits can be off by one; every combination of +-1
   cycle per gap that still fits is reported.
4. ``chain_events``: three or more single events with consistent depth and
   duration, each within a tolerance of ``t0 + n*P``, give a starting ephemeris
   when TLS finds nothing (it needs the transits to line up exactly).

Quasi-periodic search (QATS; Carter & Agol 2013) is a later option for
transits too shallow to fit one by one.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import chi2 as chi2_dist

from quicklook.duo import events_consistent
from quicklook.single_event import SingleEvent

INGRESS_FRACTION = 0.15  # ingress/egress length as a fraction of the total duration
TTV_DELTA_DURATIONS = 1.0  # search +-delta = this many durations around the prediction ...
TTV_DELTA_PERIOD_FRACTION = 0.25  # ... but never more than this fraction of the period
TTV_SNR_MIN = 3.0  # per-transit depth / depth error
TTV_PVALUE_MAX = 1e-3  # linear ephemeris rejected below this p-value
TTV_MIN_TRANSITS = 3
MIN_COVERAGE = 0.5  # fraction of expected cadences in transit and on each side
DURATION_RANGE = (0.5, 2.5)  # fitted duration bounds, in units of the input duration
_EDGE_FRACTION = 0.98  # a fit pinned at the edge of +-delta is not trusted
_GRID_STEPS_PER_DURATION = 30
CYCLE_GAP_MIN = 30.0  # days between fitted transits that make the cycle count uncertain
CYCLE_DCHI2_MAX = 9.0  # alternative cycle counts within this chi2 of the best are kept
CYCLE_MAX_GAPS = 4  # 3**4 = 81 combinations at most
CHAIN_MIN_EVENTS = 3
CHAIN_TOL_DURATIONS = 1.0  # each event within this many durations of t0 + n*P
CHAIN_P_MIN = 0.5  # days
CHAIN_MAX_CYCLES = 500


@dataclass(frozen=True)
class TransitTime:
    """One fitted transit. Times and durations in days, depth as a fraction."""

    epoch: int
    tc: float
    tc_err: float
    duration: float
    depth: float
    snr: float


@dataclass(frozen=True)
class Ephemeris:
    t0: float  # mid-time at epoch 0
    period: float
    t0_err: float
    period_err: float
    chi2: float
    dof: int


@dataclass(frozen=True)
class TTVStats:
    ephemeris: Ephemeris | None
    epochs: np.ndarray = field(repr=False)
    oc: np.ndarray = field(repr=False)  # observed - calculated (days)
    oc_err: np.ndarray = field(repr=False)
    amplitude: float  # half the O-C peak-to-peak (days)
    rms: float  # O-C root mean square (days)
    pvalue: float  # chi-square p-value of the linear ephemeris
    pvalue_robust: float  # same, without the transit with the largest |O-C|/error
    detected: bool

    @property
    def n(self) -> int:
        return int(self.epochs.size)


@dataclass(frozen=True)
class CycleSolution:
    """Linear ephemeris for one choice of cycle-count shifts across long gaps."""

    shifts: tuple[int, ...]  # cumulative epoch shift added after each gap
    period: float
    t0: float
    chi2: float


# --- model ----------------------------------------------------------------------------


def trapezoid(time, tc: float, duration: float) -> np.ndarray:
    """Unit-depth trapezoid: 1 at the centre, 0 out of transit."""
    x = np.abs(np.asarray(time, dtype=float) - tc)
    half = duration / 2
    ingress = max(INGRESS_FRACTION * duration, 1e-6)
    return np.clip((half - x) / ingress, 0.0, 1.0)


def _cadence(time) -> float:
    return float(np.median(np.diff(time)))


def _scatter(flux) -> float:
    d = np.diff(flux)
    return float(1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2))


def _covered(time, tc, duration, cadence) -> bool:
    """Enough cadences in transit and within one duration on each side."""
    half = duration / 2
    need = MIN_COVERAGE * duration / cadence
    edges = [tc - half - duration, tc - half, tc + half, tc + half + duration]
    a, b, c, d = np.searchsorted(time, edges)
    return (c - b) >= need and (b - a) >= need and (d - c) >= need


def _linear_depth(t, f, tc, duration, t_ref):
    """Least-squares depth for a box-free trapezoid at ``tc`` on a linear baseline."""
    design = np.column_stack([np.ones_like(t), t - t_ref, -trapezoid(t, tc, duration)])
    coef, *_ = np.linalg.lstsq(design, f, rcond=None)
    resid = f - design @ coef
    return coef, float(resid @ resid)


def _coarse_tc(t, f, tc_pred, duration, delta, cadence):
    step = max(cadence, duration / _GRID_STEPS_PER_DURATION)
    best = (np.inf, tc_pred, None)
    for tc in np.arange(tc_pred - delta, tc_pred + delta + step / 2, step):
        coef, ss = _linear_depth(t, f, tc, duration, tc_pred)
        if coef[2] > 0 and ss < best[0]:
            best = (ss, tc, coef)
    return best[1], best[2]


def _refine(t, f, tc_pred, tc0, duration0, coef0, delta, sigma):
    def model(p):
        tc, dur, depth, a, b = p
        return (a + b * (t - tc_pred)) * (1 - depth * trapezoid(t, tc, dur))

    x0 = np.array([tc0, duration0, coef0[2] / max(coef0[0], 1e-9), coef0[0], coef0[1]])
    lo = [tc_pred - delta, DURATION_RANGE[0] * duration0, 0.0, -np.inf, -np.inf]
    hi = [tc_pred + delta, DURATION_RANGE[1] * duration0, 1.0, np.inf, np.inf]
    x0 = np.clip(x0, np.array(lo) + 1e-12, np.array(hi) - 1e-12)
    sol = least_squares(lambda p: (f - model(p)) / sigma, x0, bounds=(lo, hi), x_scale="jac")
    dof = max(t.size - x0.size, 1)
    red = max(1.0, float(sol.fun @ sol.fun) / dof)
    try:
        cov = np.linalg.inv(sol.jac.T @ sol.jac) * red
        err = np.sqrt(np.clip(np.diag(cov), 0, np.inf))
    except np.linalg.LinAlgError:
        err = np.full(x0.size, np.inf)
    return sol.x, err


def fit_transit_time(
    time, flux, tc_pred: float, duration: float, delta: float, epoch: int = 0
) -> TransitTime | None:
    """Fit one transit near ``tc_pred``; None if not covered or not significant.

    The mid-time is free within ``tc_pred +- delta``, the duration within
    ``DURATION_RANGE`` x ``duration``. Errors come from the Jacobian, scaled up by the
    reduced chi-square when the local scatter exceeds the white noise.
    """
    t, f = np.asarray(time, dtype=float), np.asarray(flux, dtype=float)
    ok = np.isfinite(t) & np.isfinite(f)
    t, f = t[ok], f[ok]
    if t.size < 10:
        return None
    cadence = _cadence(t)
    # Room for the longest allowed duration plus a baseline as long on each side.
    half_window = delta + 1.5 * DURATION_RANGE[1] * duration
    sel = np.abs(t - tc_pred) < half_window
    t, f = t[sel], f[sel]
    if t.size < 10 or not _covered(t, tc_pred, duration, cadence):
        return None
    tc0, coef = _coarse_tc(t, f, tc_pred, duration, delta, cadence)
    if coef is None:
        return None
    sigma = max(_scatter(f), 1e-9)
    (tc, dur, depth, *_), err = _refine(t, f, tc_pred, tc0, duration, coef, delta, sigma)
    tc_err, depth_err = float(err[0]), float(err[2])
    snr = depth / depth_err if depth_err > 0 else 0.0
    pinned = abs(tc - tc_pred) > _EDGE_FRACTION * delta
    if snr < TTV_SNR_MIN or pinned or not np.isfinite(tc_err):
        return None
    if not _covered(t, tc, dur, cadence):
        return None
    return TransitTime(int(epoch), float(tc), tc_err, float(dur), float(depth), float(snr))


def default_delta(period: float, duration: float) -> float:
    """Half-width (days) of the window searched around each predicted time."""
    return float(min(TTV_DELTA_DURATIONS * duration, TTV_DELTA_PERIOD_FRACTION * period))


def fit_transit_times(
    time, flux, period: float, t0: float, duration: float, delta: float | None = None
) -> tuple[TransitTime, ...]:
    """Fit every covered transit of the ephemeris ``t0 + n*period``.

    Epoch numbers are relative to ``t0`` (which may be negative).
    """
    t = np.asarray(time, dtype=float)
    finite = t[np.isfinite(t)]
    if finite.size < 10 or not (period > 0 and duration > 0):
        return ()
    delta = default_delta(period, duration) if delta is None else float(delta)
    n_lo = int(np.ceil((finite.min() - t0) / period))
    n_hi = int(np.floor((finite.max() - t0) / period))
    fits = (
        fit_transit_time(time, flux, t0 + n * period, duration, delta, epoch=n)
        for n in range(n_lo, n_hi + 1)
    )
    return tuple(tt for tt in fits if tt is not None)


# --- ephemeris and O-C ----------------------------------------------------------------


def _weighted_line(epochs, tcs, errs):
    w = 1 / np.asarray(errs, dtype=float) ** 2
    design = np.column_stack([np.ones_like(epochs, dtype=float), epochs])
    cov = np.linalg.inv(design.T @ (design * w[:, None]))
    coef = cov @ (design.T @ (w * tcs))
    resid = tcs - design @ coef
    return coef, cov, float(np.sum(w * resid**2))


def linear_ephemeris(times) -> Ephemeris:
    """Weighted least-squares ``tc = t0 + period * epoch``."""
    epochs = np.array([tt.epoch for tt in times], dtype=float)
    if np.unique(epochs).size < 2:
        raise ValueError("A linear ephemeris needs at least two distinct epochs.")
    tcs = np.array([tt.tc for tt in times])
    errs = np.array([tt.tc_err for tt in times])
    (t0, period), cov, chi2 = _weighted_line(epochs, tcs, errs)
    return Ephemeris(
        float(t0),
        float(period),
        float(np.sqrt(cov[0, 0])),
        float(np.sqrt(cov[1, 1])),
        chi2,
        int(epochs.size - 2),
    )


def ttv_statistics(times) -> TTVStats:
    """O-C against the best linear ephemeris and its chi-square p-value."""
    epochs = np.array([tt.epoch for tt in times], dtype=int)
    errs = np.array([tt.tc_err for tt in times], dtype=float)
    if epochs.size < TTV_MIN_TRANSITS or np.unique(epochs).size < TTV_MIN_TRANSITS:
        empty = np.array([])
        return TTVStats(None, epochs, empty, errs, np.nan, np.nan, np.nan, np.nan, False)
    eph = linear_ephemeris(times)
    tcs = np.array([tt.tc for tt in times])
    oc = tcs - (eph.t0 + eph.period * epochs)
    pvalue = float(chi2_dist.sf(eph.chi2, eph.dof))
    robust = _pvalue_without_worst(times, oc / errs)
    return TTVStats(
        eph,
        epochs,
        oc,
        errs,
        float((oc.max() - oc.min()) / 2),
        float(np.sqrt(np.mean(oc**2))),
        pvalue,
        robust,
        bool(pvalue < TTV_PVALUE_MAX and robust < TTV_PVALUE_MAX),
    )


def _pvalue_without_worst(times, pulls) -> float:
    """Linear-ephemeris p-value after dropping the largest |pull|.

    Per-transit timing errors have heavier tails than a Gaussian (a noise
    spike in ingress can drag one fit), so a single outlier must not be enough
    to claim TTVs. With only three transits there is nothing to drop.
    """
    if len(times) <= TTV_MIN_TRANSITS:
        return float(chi2_dist.sf(linear_ephemeris(times).chi2, len(times) - 2))
    worst = int(np.argmax(np.abs(pulls)))
    rest = [tt for i, tt in enumerate(times) if i != worst]
    eph = linear_ephemeris(rest)
    return float(chi2_dist.sf(eph.chi2, eph.dof))


def ttv_mask(time, times, margin: float = 1.2) -> np.ndarray:
    """Cadences within ``margin`` x half the fitted duration of each fitted mid-time."""
    t = np.asarray(time, dtype=float)
    mask = np.zeros(t.size, dtype=bool)
    for tt in times:
        mask |= np.abs(t - tt.tc) < tt.duration * margin / 2
    return mask


# --- cycle-count ambiguity ------------------------------------------------------------


def _segment_starts(tcs, epochs, gap_min):
    """Indices (into time-sorted transits) that start a new segment after a long gap.

    A gap only makes the cycle count uncertain if epochs are missing across it:
    consecutive transits of a long-period planet are far apart but unambiguous.
    """
    return [
        i + 1
        for i, (dt, de) in enumerate(zip(np.diff(tcs), np.diff(epochs)))
        if dt > gap_min and de > 1
    ]


def cycle_count_solutions(times, gap_min: float = CYCLE_GAP_MIN) -> list[CycleSolution]:
    """Every +-1 cycle shift per long gap whose linear fit is within ``CYCLE_DCHI2_MAX``.

    The first entry is the best fit. With no long gaps the only solution is the
    one given (``shifts == ()``).
    """
    order = np.argsort([tt.tc for tt in times])
    tcs = np.array([times[i].tc for i in order])
    epochs = np.array([times[i].epoch for i in order], dtype=float)
    errs = np.array([times[i].tc_err for i in order])
    starts = _segment_starts(tcs, epochs, gap_min)[:CYCLE_MAX_GAPS]
    solutions = []
    for shifts in itertools.product((0, -1, 1), repeat=len(starts)):
        shifted = epochs.copy()
        for start, shift in zip(starts, shifts):
            shifted[start:] += shift
        if np.any(np.diff(shifted) <= 0):
            continue
        try:
            (t0, period), _, chi2 = _weighted_line(shifted, tcs, errs)
        except np.linalg.LinAlgError:  # zero/inf weights (WASP-93)
            continue
        if period > 0:
            solutions.append(CycleSolution(tuple(shifts), float(period), float(t0), chi2))
    if not solutions:
        return []
    best = min(s.chi2 for s in solutions)
    kept = [s for s in solutions if s.chi2 - best <= CYCLE_DCHI2_MAX]
    return sorted(kept, key=lambda s: (s.chi2, any(s.shifts)))


# --- chaining single events -----------------------------------------------------------


def _consistent_members(events):
    ref = max(events, key=lambda e: e.snr)
    return [e for e in events if e is ref or events_consistent(ref, e)]


def _candidate_periods(t0s, p_floor):
    for i, j in itertools.combinations(range(t0s.size), 2):
        span = abs(t0s[j] - t0s[i])
        n_max = min(int(span // p_floor), CHAIN_MAX_CYCLES)
        for n in range(1, n_max + 1):
            yield span / n, t0s[min(i, j)]


def _members_at(events, period, t_ref):
    matched = {}
    for ev in events:
        n = int(np.round((ev.t0 - t_ref) / period))
        resid = abs(ev.t0 - t_ref - n * period)
        if resid <= CHAIN_TOL_DURATIONS * ev.duration and n not in matched:
            matched[n] = ev
    return matched


def chain_events(events) -> tuple[float, float, tuple[SingleEvent, ...]] | None:
    """Longest period that places >= 3 consistent events within tolerance of ``t0 + n*P``.

    Returns ``(period, t0, members)`` from a linear fit through the members, or
    None. Like a duotransit, an integer fraction of the period also fits; the
    longest one is returned, so check for missing transits in covered data.
    """
    events = sorted(events, key=lambda e: e.t0)
    if len(events) < CHAIN_MIN_EVENTS:
        return None
    members = sorted(_consistent_members(events), key=lambda e: e.t0)
    if len(members) < CHAIN_MIN_EVENTS:
        return None
    t0s = np.array([e.t0 for e in members])
    p_floor = max(CHAIN_P_MIN, 3 * max(e.duration for e in members))
    best = (0, 0.0, None)
    for period, t_ref in _candidate_periods(t0s, p_floor):
        matched = _members_at(members, period, t_ref)
        if (len(matched), period) > best[:2]:
            best = (len(matched), period, (t_ref, matched))
    count, period, payload = best
    if count < CHAIN_MIN_EVENTS:
        return None
    t_ref, matched = payload
    epochs = np.array(sorted(matched), dtype=float)
    tcs = np.array([matched[int(n)].t0 for n in epochs])
    (t0, p_fit), _, _ = _weighted_line(epochs - epochs[0], tcs, np.ones_like(tcs))
    return float(p_fit), float(t0), tuple(matched[int(n)] for n in epochs)


__all__ = [
    "CycleSolution",
    "Ephemeris",
    "TTVStats",
    "TransitTime",
    "chain_events",
    "cycle_count_solutions",
    "default_delta",
    "fit_transit_time",
    "fit_transit_times",
    "linear_ephemeris",
    "trapezoid",
    "ttv_mask",
    "ttv_statistics",
]
