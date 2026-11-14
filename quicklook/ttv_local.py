"""Local-baseline transit timing for active stars (TTV diagnostics, follow-up to #42).

The per-transit fits in ``quicklook.ttv`` work on the flattened light curve. On
young, spotted stars the global detrending leaves residuals larger than the
transit, so here every transit is fitted on the *raw* (normalised) flux with a
trapezoid plus its own cubic baseline over ``+-HALF_WINDOW`` days, all linear in
the depth and baseline coefficients. A grid over the mid-time offset then gives:

* a duration x offset map of all transits stacked (is the ephemeris right?),
* one mean offset per visit (sector or stitched block of sectors),
* one offset per transit,
* timing errors calibrated by injecting transits at random phases into the same
  light curve, because formal errors on spotted stars are several times too small.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.stats import chi2 as chi2_dist

from quicklook.ttv import trapezoid

MINUTES_PER_DAY = 1440.0
HALF_WINDOW = 0.30  # days of data on each side of a transit for its baseline
BASELINE_ORDER = 3
OFFSET_MAX = 0.25  # days (6 h); also capped at a quarter period
OFFSET_STEP = 2.0 / MINUTES_PER_DAY
DURATION_GRID_H = np.arange(0.75, 6.01, 0.25)
GRID_STEP = 4.0 / MINUTES_PER_DAY
SNR_MIN = 3.0
MIN_IN_TRANSIT = 10  # cadences within ~1 h of the predicted mid-time
MIN_WINDOW_FRACTION = 0.6
VISIT_GAP = 5.0  # days; longer gaps start a new visit
N_INJECTIONS = 10  # per visit; 0 disables calibration
INJ_SHIFT_MAX = 90.0 / MINUTES_PER_DAY
CATASTROPHIC = 60.0 / MINUTES_PER_DAY


@dataclass(frozen=True)
class OffsetFit:
    """Best common offset of one or more transits from their predicted times."""

    offset: float  # days
    err: float  # days, Delta-chi2 contour inflated by sqrt(red. chi2) * beta
    snr: float
    depth: float
    beta: float
    edge: bool
    second_snr: float
    offsets: np.ndarray = field(repr=False)
    snr_curve: np.ndarray = field(repr=False)

    @property
    def usable(self) -> bool:
        return self.snr >= SNR_MIN and not self.edge


@dataclass(frozen=True)
class VisitTiming:
    visit: int
    t_mid: float
    epochs: tuple[int, ...]
    fit: OffsetFit | None
    err_calibrated: float = np.nan  # days, from injections; NaN if not run

    @property
    def err(self) -> float:
        cal = self.err_calibrated
        return max(self.fit.err, cal) if np.isfinite(cal) else self.fit.err


@dataclass(frozen=True)
class SingleTiming:
    visit: int
    epoch: int
    fit: OffsetFit
    err_calibrated: float = np.nan

    @property
    def err(self) -> float:
        cal = self.err_calibrated
        return max(self.fit.err, cal) if np.isfinite(cal) else self.fit.err


@dataclass(frozen=True)
class Calibration:
    single_sigma: float  # days (68th percentile of |recovered - injected|)
    stack_sigma: float
    single_catastrophic: float  # fraction off by more than CATASTROPHIC
    n_single: int
    n_stack: int


@dataclass
class LocalTTVResult:
    period: float
    t0: float
    duration: float
    source: str
    grid_offsets: np.ndarray
    grid_durations: np.ndarray
    grid_snr: np.ndarray
    visits: list[VisitTiming]
    singles: list[SingleTiming]
    calibration: dict[int, Calibration]
    visit_chi2: float = np.nan
    visit_dof: int = 0
    visit_pvalue: float = np.nan

    @property
    def best_common_offset(self) -> float:
        i, j = np.unravel_index(np.nanargmax(self.grid_snr), self.grid_snr.shape)
        return float(self.grid_offsets[j])

    def usable_visits(self) -> list[VisitTiming]:
        return [v for v in self.visits if v.fit is not None and v.fit.usable]

    def usable_singles(self) -> list[SingleTiming]:
        return [s for s in self.singles if s.fit.usable]


# --- light-curve helpers -------------------------------------------------------------


def split_visits(time, gap: float = VISIT_GAP) -> list[np.ndarray]:
    """Index arrays of contiguous blocks separated by more than ``gap`` days."""
    t = np.asarray(time, dtype=float)
    order = np.argsort(t)
    cuts = np.where(np.diff(t[order]) > gap)[0] + 1
    return [np.sort(ix) for ix in np.split(order, cuts)]


def covered_epochs(time, period: float, t0: float) -> list[int]:
    """Epochs with enough cadences in the baseline window and near mid-transit."""
    t = np.asarray(time, dtype=float)
    if t.size == 0 or not (np.isfinite(period) and np.isfinite(t0) and period > 0):
        return []
    cadence = np.median(np.diff(np.sort(t))) if t.size > 1 else 2 / MINUTES_PER_DAY
    expected = 2 * HALF_WINDOW / cadence
    lo = int(np.ceil((t.min() - t0) / period))
    hi = int(np.floor((t.max() - t0) / period))
    out = []
    for n in range(lo, hi + 1):
        tc = t0 + n * period
        win = np.abs(t - tc) < HALF_WINDOW
        near = np.abs(t - tc) < 0.04
        if win.sum() > MIN_WINDOW_FRACTION * expected and near.sum() >= MIN_IN_TRANSIT:
            out.append(n)
    return out


def _robust_sigma(flux) -> float:
    f = np.asarray(flux, dtype=float)
    d = np.diff(f[np.isfinite(f)])
    return float(1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2)) if d.size else np.nan


# --- joint linear fit ---------------------------------------------------------------


def _segment(t, f, tc, duration, keep):
    """Baseline-projected flux and transit model for one transit window (or None)."""
    lo, hi = np.searchsorted(t, [tc - HALF_WINDOW, tc + HALF_WINDOW])
    tt, ff, kk = t[lo:hi], f[lo:hi], keep[lo:hi]
    tt, ff = tt[kk], ff[kk]
    if tt.size < BASELINE_ORDER + 3:
        return None
    q, _ = np.linalg.qr(np.vander((tt - tc) / HALF_WINDOW, BASELINE_ORDER + 1))
    y = ff - q @ (q.T @ ff)
    m = -trapezoid(tt, tc, duration)
    m = m - q @ (q.T @ m)
    return y, m


def joint_fit(t, f, tcs, duration, keep=None):
    """Common depth + a cubic baseline per transit; returns (depth, dchi2_raw, resid, dof).

    Each baseline is projected out of its own window, which is equivalent to the
    joint linear least-squares fit but costs O(N) instead of O(N * n_transits^2).
    ``t`` must be sorted.
    """
    keep = np.ones(t.size, bool) if keep is None else keep
    segs = [s for s in (_segment(t, f, tc, duration, keep) for tc in tcs) if s is not None]
    if not segs:
        return None
    y = np.concatenate([s[0] for s in segs])
    m = np.concatenate([s[1] for s in segs])
    if y.size < 30:
        return None
    mm = float(m @ m)
    if mm <= 0:
        return None
    my = float(m @ y)
    depth = my / mm
    resid = y - depth * m
    return depth, my * my / mm, resid, y.size - len(segs) * (BASELINE_ORDER + 1) - 1


def _beta(resid, duration, cadence):
    """Pont et al. (2006) time-averaging factor at the transit-duration scale."""
    n = max(int(duration / cadence), 2)
    m = resid.size // n
    if m < 4:
        return 1.0
    binned = resid[: m * n].reshape(m, n).mean(axis=1)
    expected = np.std(resid) / np.sqrt(n) * np.sqrt(m / (m - 1))
    return float(max(np.std(binned) / expected, 1.0)) if expected > 0 else 1.0


def scan_offsets(t, f, tcs, duration, sigma, offsets, keep=None):
    """Signed Delta-chi2 (in units of sigma^2) and depth at each common offset."""
    dchi = np.full(offsets.size, np.nan)
    depth = np.full(offsets.size, np.nan)
    for i, o in enumerate(offsets):
        r = joint_fit(t, f, [tc + o for tc in tcs], duration, keep)
        if r is None:
            continue
        depth[i] = r[0]
        dchi[i] = np.sign(r[0]) * r[1] / sigma**2
    return dchi, depth


def offset_grid(period: float) -> np.ndarray:
    half = min(OFFSET_MAX, 0.25 * period)
    return np.arange(-half, half + 1e-12, OFFSET_STEP)


def measure_offset(t, f, tcs, duration, sigma, keep=None, offsets=None) -> OffsetFit | None:
    """Best common offset of ``tcs`` with an inflated Delta-chi2 = 1 error."""
    offsets = offset_grid(np.inf) if offsets is None else offsets
    dchi, depth = scan_offsets(t, f, tcs, duration, sigma, offsets, keep)
    if not np.isfinite(dchi).any():
        return None
    best = int(np.nanargmax(dchi))
    r = joint_fit(t, f, [tc + offsets[best] for tc in tcs], duration, keep)
    cadence = float(np.median(np.diff(np.sort(t))))
    red = np.sum(r[2] ** 2) / sigma**2 / max(r[3], 1)
    infl = np.sqrt(max(red, 1.0)) * _beta(r[2], duration, cadence)
    thresh = dchi[best] - infl**2
    lo = hi = best
    while lo > 0 and np.isfinite(dchi[lo - 1]) and dchi[lo - 1] > thresh:
        lo -= 1
    while hi < offsets.size - 1 and np.isfinite(dchi[hi + 1]) and dchi[hi + 1] > thresh:
        hi += 1
    step = offsets[1] - offsets[0] if offsets.size > 1 else OFFSET_STEP
    err = max(0.5 * (offsets[hi] - offsets[lo]), step / 2)
    curve = np.sqrt(np.clip(dchi, 0, None)) / infl
    far = np.abs(offsets - offsets[best]) > 2 * duration
    return OffsetFit(
        offset=float(offsets[best]),
        err=float(err),
        snr=float(curve[best]),
        depth=float(depth[best]),
        beta=float(infl),
        edge=bool(lo == 0 or hi == offsets.size - 1),
        second_snr=float(np.nanmax(curve[far])) if far.any() else np.nan,
        offsets=offsets,
        snr_curve=curve,
    )


# --- pipeline -------------------------------------------------------------------------


def duration_offset_map(t, f, period, t0, epochs, sigma, keep=None):
    """sqrt(Delta-chi2) of all ``epochs`` stacked over (duration, common offset)."""
    half = min(OFFSET_MAX, 0.25 * period)
    offsets = np.arange(-half, half + 1e-12, GRID_STEP)
    durations = DURATION_GRID_H / 24
    tcs = [t0 + n * period for n in epochs]
    grid = np.full((durations.size, offsets.size), np.nan)
    for i, d in enumerate(durations):
        dchi, _ = scan_offsets(t, f, tcs, d, sigma, offsets, keep)
        grid[i] = np.sqrt(np.clip(dchi, 0, None))
    return offsets, durations, grid


def _injection_errors(t, f, period, t0, duration, depth, sigma, avoid, offsets, n_inj, rng):
    """Inject at random phases away from ``avoid`` and return |recovered - injected| arrays."""
    single, stack = [], []
    for _ in range(n_inj * 4):
        if len(stack) >= n_inj:
            break
        t0p = t0 + rng.uniform(0, period)
        epochs = [
            n
            for n in covered_epochs(t, period, t0p)
            if not avoid(t0p + n * period, HALF_WINDOW + offsets[-1])
        ]
        if not epochs:
            continue
        shift = rng.uniform(-INJ_SHIFT_MAX, INJ_SHIFT_MAX)
        tcs = [t0p + n * period for n in epochs]
        fi = f.copy()
        for tc in tcs:
            fi = fi - depth * trapezoid(t, tc + shift, duration)
        m = measure_offset(t, fi, tcs, duration, sigma, offsets=offsets)
        if m is not None and m.usable:
            stack.append(m.offset - shift)
        for tc in tcs:
            mi = measure_offset(t, fi, [tc], duration, sigma, offsets=offsets)
            if mi is not None and mi.usable:
                single.append(mi.offset - shift)
    return np.asarray(single), np.asarray(stack)


def _sigma68(x):
    return float(np.percentile(np.abs(x - np.median(x)), 68.27)) if x.size >= 5 else np.nan


def run_local_ttv(
    time,
    flux,
    period: float,
    t0: float,
    duration: float | None = None,
    source: str = "",
    n_injections: int = N_INJECTIONS,
    seed: int = 0,
) -> LocalTTVResult | None:
    """Duration/offset map, per-visit and per-transit offsets, calibrated errors.

    ``flux`` is the normalised, *un-flattened* light curve. ``duration`` (days) is
    only a fallback: the best duration of the stacked map is used.
    """
    t, f = np.asarray(time, dtype=float), np.asarray(flux, dtype=float)
    ok = np.isfinite(t) & np.isfinite(f)
    order = np.argsort(t[ok])
    t, f = t[ok][order], f[ok][order]
    epochs = covered_epochs(t, period, t0)
    if not epochs:
        return None
    sigma = _robust_sigma(f)
    g_off, g_dur, grid = duration_offset_map(t, f, period, t0, epochs, sigma)
    if np.isfinite(grid).any():
        i, _ = np.unravel_index(np.nanargmax(grid), grid.shape)
        duration = float(g_dur[i])
    elif duration is None:
        return None
    offsets = offset_grid(period)

    def avoid(tc, pad):
        ph = (tc - t0 + 0.5 * period) % period - 0.5 * period
        return abs(ph) < pad + duration / 2

    rng = np.random.default_rng(seed)
    visits, singles, calib = [], [], {}
    for k, ix in enumerate(split_visits(t)):
        tv, fv = t[ix], f[ix]
        ep = covered_epochs(tv, period, t0)
        if not ep:
            continue
        sig_v = _robust_sigma(fv)
        tcs = [t0 + n * period for n in ep]
        fit = measure_offset(tv, fv, tcs, duration, sig_v, offsets=offsets)
        cal = None
        if n_injections > 0 and fit is not None:
            depth = fit.depth if fit.depth > 0 else 1e-3
            s1, s2 = _injection_errors(
                tv, fv, period, t0, duration, depth, sig_v, avoid, offsets, n_injections, rng
            )
            cal = Calibration(
                single_sigma=_sigma68(s1),
                stack_sigma=_sigma68(s2),
                single_catastrophic=float(np.mean(np.abs(s1) > CATASTROPHIC))
                if s1.size
                else np.nan,
                n_single=int(s1.size),
                n_stack=int(s2.size),
            )
            calib[k] = cal
        visits.append(
            VisitTiming(k, float(np.median(tv)), tuple(ep), fit, cal.stack_sigma if cal else np.nan)
        )
        for n, tc in zip(ep, tcs):
            mi = measure_offset(tv, fv, [tc], duration, sig_v, offsets=offsets)
            if mi is not None:
                singles.append(SingleTiming(k, n, mi, cal.single_sigma if cal else np.nan))
    result = LocalTTVResult(
        period, t0, duration, source, g_off, g_dur, grid, visits, singles, calib
    )
    _visit_linear_test(result)
    return result


def _visit_linear_test(result: LocalTTVResult) -> None:
    """chi2 of a weighted line through the visit-mean offsets (needs >= 3 visits)."""
    vs = result.usable_visits()
    if len(vs) < 3:
        return
    ep = np.array([np.mean(v.epochs) for v in vs])
    oc = np.array([v.fit.offset for v in vs])
    er = np.array([v.err for v in vs])
    A = np.vstack([np.ones_like(ep), ep]).T
    w = 1 / er**2
    coef = np.linalg.solve(A.T @ (A * w[:, None]), A.T @ (w * oc))
    chi = float(np.sum(((oc - A @ coef) / er) ** 2))
    result.visit_chi2 = chi
    result.visit_dof = len(vs) - 2
    result.visit_pvalue = float(chi2_dist.sf(chi, result.visit_dof))


__all__ = [
    "Calibration",
    "LocalTTVResult",
    "OffsetFit",
    "SingleTiming",
    "VisitTiming",
    "covered_epochs",
    "joint_fit",
    "measure_offset",
    "run_local_ttv",
    "split_visits",
]
