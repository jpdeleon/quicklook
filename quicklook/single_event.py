"""Single-event (mono-transit) search on a flattened light curve.

Periodic searches (TLS, BLS) need about three transits, so a planet with one
transit in the data is never found by them. This module looks for individual
dips instead, in the spirit of Foreman-Mackey et al. (2016, AJ 152, 206):

1. Slide a box over the light curve for a grid of durations set by the star
   (``duration_grid``), giving a detection SNR at every trial mid-time.
2. Keep the strongest non-overlapping peaks above ``snr_min``.
3. Compare each peak locally against systematics models (a flux step and a
   quadratic ramp) by BIC, and keep it only if a transit wins clearly.

Events at a data edge or gap edge (no baseline on one side) are excluded.
``min_period`` turns an event into a period lower limit: the smallest period
whose other transits all fall outside the covered data.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

SNR_MIN = 7.0
DELTA_BIC_MIN = 10.0
MAX_EVENTS = 5
N_DURATIONS = 12
MAX_DURATION = 1.5  # days; longer dips do not survive typical detrending
MIN_DURATION = 0.02  # days (~30 min)
DEFAULT_P_MAX = 1000.0  # days, sets the longest trial duration
CENTRAL_DURATION_1YR = 13.0 / 24.0  # days; central transit, P=1 yr, solar density
_ECC_ALLOWANCE = 1.3  # eccentric orbits can transit longer than circular
_GRAZING_FRACTION = 0.3  # high-impact transits are shorter than central ones
_GRID_STEP_FRACTION = 1 / 3  # trial-t0 spacing, in units of the shortest duration
_MIN_BASELINE = 0.2  # days of out-of-transit window on each side (at least)
_MIN_SIDE = 0.1  # days of data required on each side of a trial ...
_MIN_SIDE_FRACTION = 0.25  # ... or this fraction of the duration, if longer
_N_SHIFT_TRIALS = 9  # transit-centre shifts tried when refining an event
_REFINE_RANGE = (0.5, 3.0)  # durations tried when refining, relative to the grid peak
_MAX_REFINE_ITER = 4
_N_REFINE_DURATIONS = 15  # finer than the search grid
_MIN_COVERAGE = 0.5  # fraction of expected cadences required in each window
_N_STEP_TRIALS = 21
_GAP_CADENCES = 10  # a gap is a jump longer than this many cadences ...
_GAP_MIN = 0.1  # ... and at least this many days


@dataclass(frozen=True)
class SingleEvent:
    """One transit-like dip. Times and durations in days, depth as a fraction."""

    t0: float
    duration: float
    depth: float
    snr: float
    delta_bic_null: float  # BIC(flat baseline) - BIC(transit); > 0 favours transit
    delta_bic_systematics: float  # BIC(best step/ramp) - BIC(transit)
    dist_to_gap: float  # to the nearest data gap or edge
    n_in: int


@dataclass(frozen=True)
class SingleEventResult:
    events: tuple[SingleEvent, ...]
    grid_t0: np.ndarray = field(repr=False)
    grid_snr: np.ndarray = field(repr=False)
    grid_duration: np.ndarray = field(repr=False)
    durations: np.ndarray = field(repr=False)
    snr_min: float = SNR_MIN


def central_duration(period: float, rho_star: float = 1.0) -> float:
    """Central-transit duration (days) for a circular orbit; rho_star in solar units."""
    return CENTRAL_DURATION_1YR * (period / 365.25) ** (1 / 3) * rho_star ** (-1 / 3)


def duration_grid(
    span: float,
    cadence: float,
    rho_star: float = 1.0,
    p_min: float | None = None,
    p_max: float = DEFAULT_P_MAX,
    n: int = N_DURATIONS,
) -> np.ndarray:
    """Log-spaced trial durations (days) from grazing short-period to eccentric long-period.

    Talens, Aigrain et al. (2025, RNAAS) show long-period transits can be much
    shorter or longer than the usual 1-12 h, so the range follows the star.
    """
    rho = rho_star if np.isfinite(rho_star) and rho_star > 0 else 1.0
    p_lo = p_min if p_min is not None else max(span / 2, 1.0)
    d_min = max(3 * cadence, MIN_DURATION, _GRAZING_FRACTION * central_duration(p_lo, rho))
    d_max = min(MAX_DURATION, _ECC_ALLOWANCE * central_duration(p_max, rho))
    d_max = max(d_max, 2 * d_min)
    return np.geomspace(d_min, d_max, n)


def _prepare(time, flux):
    time, flux = np.asarray(time, dtype=float), np.asarray(flux, dtype=float)
    ok = np.isfinite(time) & np.isfinite(flux)
    order = np.argsort(time[ok])
    return time[ok][order], flux[ok][order]


def _cadence(time) -> float:
    return float(np.median(np.diff(time)))


def _point_scatter(flux) -> float:
    """White-noise level from point-to-point differences (robust to dips)."""
    d = np.diff(flux)
    return float(1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2))


def _trial_grid(time, step, windows):
    grid = np.arange(time[0], time[-1], step)
    if windows:
        keep = np.zeros(grid.size, dtype=bool)
        for lo, hi in windows:
            keep |= (grid >= lo) & (grid <= hi)
        grid = grid[keep]
    return grid


def _box_depth(time, csum, grid, dur, cadence):
    """Box depth and its white-noise weight at every trial t0 (NaN where not testable).

    Returns ``(depth, inv_n)`` where ``inv_n = 1/n_in + 1/n_out``.
    """
    base = max(2 * dur, _MIN_BASELINE)
    edges = [
        grid - dur / 2 - base,
        grid - dur / 2,
        grid + dur / 2,
        grid + dur / 2 + base,
    ]
    i0, i1, i2, i3 = (np.searchsorted(time, e) for e in edges)
    n_left, n_in, n_right = i1 - i0, i2 - i1, i3 - i2
    s_left, s_in, s_right = csum[i1] - csum[i0], csum[i2] - csum[i1], csum[i3] - csum[i2]
    need_in = max(3, _MIN_COVERAGE * dur / cadence)
    # A fixed minimum, not a fraction of the (up to 2x duration) window, so a
    # long transit near a data edge stays testable while a dip with no
    # baseline at all on one side is still rejected.
    need_side = max(_MIN_SIDE, _MIN_SIDE_FRACTION * dur) / cadence
    ok = (n_in >= need_in) & (n_left >= need_side) & (n_right >= need_side)
    depth = np.full(grid.size, np.nan)
    inv_n = np.full(grid.size, np.nan)
    n_out = n_left + n_right
    depth[ok] = (s_left[ok] + s_right[ok]) / n_out[ok] - s_in[ok] / n_in[ok]
    inv_n[ok] = 1 / n_in[ok] + 1 / n_out[ok]
    return depth, inv_n


def _red_scale(time, csum, dur, cadence) -> float:
    """Robust scatter of box depths over the whole light curve at this duration.

    Captures correlated (red) noise on the transit timescale, which a
    point-to-point estimate misses (cf. Pont et al. 2006). Always measured on
    the full light curve, so a windowed search is not normalised by its own
    transit.
    """
    full = np.arange(time[0], time[-1], max(dur / 2, cadence))
    depth, _ = _box_depth(time, csum, full, dur, cadence)
    depth = depth[np.isfinite(depth)]
    if depth.size < 10:
        return 0.0
    return float(1.4826 * np.median(np.abs(depth - np.median(depth))))


def _box_snr(time, csum, grid, dur, cadence, sigma):
    """Box-search SNR at every trial t0 for one duration (NaN where not testable).

    The noise is the larger of the white-noise expectation and the red-noise
    scale from ``_red_scale``.
    """
    depth, inv_n = _box_depth(time, csum, grid, dur, cadence)
    noise = np.maximum(sigma * np.sqrt(inv_n), _red_scale(time, csum, dur, cadence))
    return depth / noise


def _bic(design, y):
    coef, *_ = np.linalg.lstsq(design, y, rcond=None)
    rss = float(np.sum((y - design @ coef) ** 2))
    n, k = design.shape
    return n * np.log(max(rss, 1e-300) / n) + k * np.log(n), coef


@dataclass(frozen=True)
class LocalModels:
    """Local fits around one event; ``time`` is absolute, models share its length."""

    time: np.ndarray
    flux: np.ndarray
    transit: np.ndarray
    systematics: np.ndarray  # the better of step and ramp
    systematics_kind: str  # "step" or "ramp"
    bic: dict
    depth: float
    n_in: int


def local_models(time, flux, t0: float, duration: float) -> LocalModels:
    """Fit flat, transit, step and quadratic-ramp models in a window around ``t0``.

    All models share a linear baseline. The step position is scanned across the
    window so a genuine discontinuity is not penalised for being off-centre.
    """
    time, flux = _prepare(time, flux)
    half = duration / 2 + max(2 * duration, _MIN_BASELINE)
    sel = np.abs(time - t0) < half
    t, y = time[sel] - t0, flux[sel]
    one = np.ones_like(t)
    box = (np.abs(t) < duration / 2).astype(float)
    designs = {
        "null": np.column_stack([one, t]),
        "transit": np.column_stack([one, t, box]),
        "ramp": np.column_stack([one, t, t**2]),
    }
    fits = {k: _bic(d, y) for k, d in designs.items()}
    steps = [
        np.column_stack([one, t, (t > s).astype(float)])
        for s in np.linspace(-half * 0.9, half * 0.9, _N_STEP_TRIALS)
    ]
    step_fits = [(_bic(d, y), d) for d in steps]
    (bic_step, step_coef), step_design = min(step_fits, key=lambda x: x[0][0])
    bic = {k: v[0] for k, v in fits.items()} | {"step": bic_step}
    kind = "step" if bic_step < bic["ramp"] else "ramp"
    sys_model = step_design @ step_coef if kind == "step" else designs["ramp"] @ fits["ramp"][1]
    transit_coef = fits["transit"][1]
    return LocalModels(
        time=time[sel],
        flux=y,
        transit=designs["transit"] @ transit_coef,
        systematics=sys_model,
        systematics_kind=kind,
        bic=bic,
        depth=-float(transit_coef[2]),
        n_in=int(box.sum()),
    )


def _has_baseline(time, center: float, dur: float, cadence: float) -> bool:
    """Same both-sides rule as the grid search, for one box."""
    need = max(_MIN_SIDE, _MIN_SIDE_FRACTION * dur)
    lo, hi = center - dur / 2, center + dur / 2
    i0, i1, i2, i3 = np.searchsorted(time, [lo - need, lo, hi, hi + need])
    n_need = _MIN_COVERAGE * need / cadence
    return (i1 - i0) >= n_need and (i3 - i2) >= n_need


def _refine_once(time, flux, t0, duration, durations, cadence):
    lo, hi = _REFINE_RANGE
    d_lo = max(lo * duration, durations[0])
    d_hi = min(hi * duration, durations[-1])
    trial_durs = np.geomspace(d_lo, d_hi, _N_REFINE_DURATIONS) if d_hi > d_lo else [duration]
    half = max(trial_durs) / 2 + max(2 * max(trial_durs), _MIN_BASELINE)
    sel = np.abs(time - t0) < half
    t, y = time[sel] - t0, flux[sel]
    one = np.ones_like(t)
    best = (np.inf, t0, duration)
    for dur in trial_durs:
        for shift in np.linspace(-duration / 2, duration / 2, _N_SHIFT_TRIALS):
            box = (np.abs(t - shift) < dur / 2).astype(float)
            if box.sum() < 3 or not _has_baseline(time, t0 + shift, dur, cadence):
                continue
            design = np.column_stack([one, t, box])
            coef, *_ = np.linalg.lstsq(design, y, rcond=None)
            rss = float(np.sum((y - design @ coef) ** 2))
            if rss < best[0] and coef[2] < 0:
                best = (rss, t0 + shift, dur)
    return best[1], best[2]


def refine_event(time, flux, t0: float, duration: float, durations) -> tuple[float, float]:
    """Best transit centre and duration near a grid peak, by least squares.

    The step model's position is scanned in ``local_models``; scanning the
    transit's centre and duration too keeps the comparison fair. The grid peak
    can be a sub-optimal box (e.g. 7 h inside a 24 h dip), which a step fits
    better than that box. Repeats until the solution stops moving, since one
    pass only reaches a few times the starting duration. Boxes without data on
    both sides are skipped, as in the grid search.
    """
    cadence = _cadence(time)
    for _ in range(_MAX_REFINE_ITER):
        new_t0, new_dur = _refine_once(time, flux, t0, duration, durations, cadence)
        if (
            np.isclose(new_dur, duration, rtol=0.05)
            and abs(new_t0 - t0) < duration / _N_SHIFT_TRIALS
        ):
            return new_t0, new_dur
        t0, duration = new_t0, new_dur
    return t0, duration


def _compare_models(time, flux, t0, dur):
    """Local BIC comparison: transit vs flat, step and quadratic-ramp baselines."""
    m = local_models(time, flux, t0, dur)
    b = m.bic
    return m.depth, b["null"] - b["transit"], min(b["step"], b["ramp"]) - b["transit"], m.n_in


def longest_segment(time) -> float:
    """Length (days) of the longest stretch of data without a gap.

    This, not the total span, bounds the shortest period a single event can
    have, so it sets the shortest trial duration for stitched multi-sector data.
    """
    b = np.sort(gap_boundaries(time))
    return float(np.max(b[1::2] - b[0::2])) if b.size >= 2 else float(np.ptp(time))


def gap_boundaries(time) -> np.ndarray:
    """Data edges plus the start and end of every gap."""
    cad = _cadence(time)
    jumps = np.flatnonzero(np.diff(time) > max(_GAP_CADENCES * cad, _GAP_MIN))
    return np.concatenate([[time[0]], time[jumps], time[jumps + 1], [time[-1]]])


def _pick_peaks(grid, snr, best_dur, snr_min, max_events):
    order = np.argsort(np.where(np.isfinite(snr), snr, -np.inf))[::-1]
    peaks = []
    for i in order:
        if not snr[i] >= snr_min or len(peaks) >= 4 * max_events:
            break
        if all(abs(grid[i] - grid[j]) > 2 * max(best_dur[i], best_dur[j]) for j in peaks):
            peaks.append(i)
    return peaks


def single_event_search(
    time,
    flux,
    rho_star: float = 1.0,
    durations=None,
    windows=None,
    snr_min: float = SNR_MIN,
    delta_bic_min: float = DELTA_BIC_MIN,
    max_events: int = MAX_EVENTS,
) -> SingleEventResult:
    """Find individual transit-like dips in a flattened light curve.

    Parameters
    ----------
    time, flux : array-like
        Flattened light curve (flux ~ 1). Non-finite points are dropped.
    rho_star : float
        Stellar density in solar units, used for the duration grid.
    durations : array-like, optional
        Trial durations in days; defaults to ``duration_grid``.
    windows : list of (lo, hi), optional
        Only test mid-times inside these intervals (e.g. around predicted
        transits).
    snr_min, delta_bic_min : float
        Detection SNR threshold and the BIC margin a transit must have over the
        best systematics model.
    max_events : int
        Maximum number of events returned, strongest first.
    """
    time, flux = _prepare(time, flux)
    empty = np.array([])
    if time.size < 20:
        return SingleEventResult((), empty, empty, empty, empty, snr_min)
    cadence = _cadence(time)
    durations = (
        duration_grid(longest_segment(time), cadence, rho_star)
        if durations is None
        else np.asarray(durations, dtype=float)
    )
    sigma = _point_scatter(flux)
    grid = _trial_grid(time, durations[0] * _GRID_STEP_FRACTION, windows)
    if grid.size == 0 or not sigma > 0:
        return SingleEventResult((), grid, np.full(grid.size, np.nan), grid, durations, snr_min)

    csum = np.concatenate([[0.0], np.cumsum(flux)])
    snr_all = np.vstack([_box_snr(time, csum, grid, d, cadence, sigma) for d in durations])
    filled = np.where(np.isfinite(snr_all), snr_all, -np.inf)
    best = np.argmax(filled, axis=0)
    grid_snr = snr_all[best, np.arange(grid.size)]
    best_dur = durations[best]

    bounds = gap_boundaries(time)
    events = []
    for i in _pick_peaks(grid, grid_snr, best_dur, snr_min, max_events):
        t0, dur = refine_event(time, flux, grid[i], best_dur[i], durations)
        depth, dbic_null, dbic_sys, n_in = _compare_models(time, flux, t0, dur)
        if depth <= 0 or dbic_sys < delta_bic_min:
            continue
        events.append(
            SingleEvent(
                t0=float(t0),
                duration=float(dur),
                depth=depth,
                snr=float(grid_snr[i]),
                delta_bic_null=float(dbic_null),
                delta_bic_systematics=float(dbic_sys),
                dist_to_gap=float(np.min(np.abs(bounds - t0))),
                n_in=n_in,
            )
        )
        if len(events) >= max_events:
            break
    return SingleEventResult(tuple(events), grid, grid_snr, best_dur, durations, snr_min)


def _covered(time, tc, duration, cadence) -> bool:
    lo, hi = np.searchsorted(time, [tc - duration / 2, tc + duration / 2])
    return (hi - lo) >= _MIN_COVERAGE * duration / cadence


def period_is_allowed(time, t0, duration, period, cadence=None, observed=(0,)) -> bool:
    """True if no transit at ``t0 + k*period`` lands in covered data, except the
    ``observed`` epochs ``k`` (the event itself, and for a duotransit its partner).
    """
    cadence = _cadence(time) if cadence is None else cadence
    k_lo = int(np.ceil((time[0] - duration - t0) / period))
    k_hi = int(np.floor((time[-1] + duration - t0) / period))
    return not any(
        _covered(time, t0 + k * period, duration, cadence)
        for k in range(k_lo, k_hi + 1)
        if k not in observed
    )


def min_period(time, t0: float, duration: float) -> float:
    """Smallest period (days) whose other transits all miss the covered data.

    Resolution is a quarter of the transit duration. A period just longer than
    the data span is always allowed, so the result is at most about the span.
    """
    time, _ = _prepare(time, np.zeros_like(np.asarray(time, dtype=float)))
    cadence = _cadence(time)
    span = time[-1] - time[0]
    for p in np.arange(duration, span + 2 * duration, duration / 4):
        if period_is_allowed(time, t0, duration, p, cadence):
            return float(p)
    return float(span + duration)


def classify_events(
    n_events: int, tls_periodic: bool, n_unexplained: int, has_duo: bool = False
) -> str:
    """``periodic`` if TLS found a real periodic signal, else ``duo``/``mono``/``none``.

    ``n_unexplained`` counts events that do not coincide with a TLS transit;
    ``has_duo`` is True when two of them pair into a duotransit. The classes
    are mutually exclusive by construction.
    """
    if tls_periodic:
        return "periodic"
    if has_duo:
        return "duo"
    if n_unexplained > 0 or n_events > 0:
        return "mono"
    return "none"
