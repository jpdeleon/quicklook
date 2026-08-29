"""Notch and LOCoR detrending, adapted for quicklook's flatten step.

Algorithm design follows the Notch filter and Locally Optimized Combination
of Rotations (LOCoR) methods described in Rizzuto et al. 2017
(https://arxiv.org/abs/1709.09670) and implemented in
https://github.com/arizzuto/Notch_and_LOCoR (``core.sliding_window`` and
``core.rcomb``). That reference implementation has no license file and pins
its per-window nonlinear fits to ``mpyfit``, a C-backed fitter that is not on
PyPI and requires a local compiler to build. To avoid that fragile
dependency, the sliding-window fits here are solved with
``scipy.optimize.least_squares`` (plus closed-form weighted-polynomial
seeding) instead of ``mpyfit``, and a few upstream implementation quirks
(see inline notes) are fixed rather than reproduced. Results are therefore
not bit-for-bit identical to the reference implementation, but follow the
same algorithmic structure.

Notch preserves transit-like dips by locally fitting a quadratic trend times
a multiplicative "notch box" within a sliding time window, comparing that
model to a plain quadratic via BIC, and dividing by whichever polynomial
component wins -- so a real transit is never divided out even when the
window is centered on it. LOCoR instead phase-folds on a rotation period and
models each rotation cycle as an optimal linear combination of neighboring
cycles, which tracks starspot evolution without needing a global periodic
model.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

__all__ = ["notch_flatten", "locor_flatten"]


# ---------------------------------------------------------------------------
# shared helpers
# ---------------------------------------------------------------------------
def _fit_quad(t: np.ndarray, y: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Weighted quadratic fit, returned as [a0, a1, a2] (low-to-high order)."""
    return np.polyfit(t, y, 2, w=w)[::-1]


def _robust_mean(y: np.ndarray, cut: float) -> tuple[float, np.ndarray]:
    """Iterative sigma-clipped mean; returns (mean, indices kept)."""
    absdev = np.abs(y - np.median(y))
    medabsdev = np.median(absdev) / 0.6745
    if medabsdev < 1e-24:
        medabsdev = np.mean(absdev) / 0.8
    gind = np.where(absdev <= cut * medabsdev)[0]
    sigma = np.std(y[gind]) if len(gind) > 1 else 0.0
    sc = max(cut, 1.0)
    if sc <= 4.5:
        denom = -0.15405 + 0.90723 * sc - 0.23584 * sc**2 + 0.020142 * sc**3
        sigma = sigma / denom if denom != 0 else sigma
    goodind = np.where(absdev <= cut * sigma)[0] if sigma > 0 else gind
    mean = np.mean(y[goodind]) if len(goodind) else np.mean(y)
    return mean, goodind


def _phase(t: np.ndarray, period: float, t0: float) -> np.ndarray:
    return np.mod((t - t0) / period, 1.0)


# ---------------------------------------------------------------------------
# Notch filter
# ---------------------------------------------------------------------------
def _notch_model_fixed_width(params: np.ndarray, t: np.ndarray, inbox: np.ndarray) -> np.ndarray:
    a0, a1, a2, depth = params
    poly = a0 + a1 * t + a2 * t**2
    box = np.where(inbox, 1.0 - depth, 1.0)
    return poly * box


def _notch_residuals_fixed_width(
    params: np.ndarray, t: np.ndarray, y: np.ndarray, sigma: np.ndarray, inbox: np.ndarray
) -> np.ndarray:
    return (y - _notch_model_fixed_width(params, t, inbox)) / sigma


def _notch_jac_fixed_width(
    params: np.ndarray, t: np.ndarray, y: np.ndarray, sigma: np.ndarray, inbox: np.ndarray
) -> np.ndarray:
    a0, a1, a2, depth = params
    poly = a0 + a1 * t + a2 * t**2
    box = np.where(inbox, 1.0 - depth, 1.0)
    d_a0 = box
    d_a1 = t * box
    d_a2 = t**2 * box
    d_depth = np.where(inbox, -poly, 0.0)
    return -(np.column_stack([d_a0, d_a1, d_a2, d_depth]) / sigma[:, None])


def _fit_notch_point(
    t: np.ndarray,
    y: np.ndarray,
    width_grid: np.ndarray,
    efrac: float,
    sigma_clip: float,
    depth_bounds: tuple[float, float],
    max_outlier_iter: int,
) -> dict | None:
    """Fit poly+notch and null (poly-only) models within one window.

    The notch width is chosen from ``width_grid`` by a closed-form seed
    search and then held fixed during the nonlinear refinement: the box
    function is discontinuous in width (points step in/out as it varies),
    which makes a continuously-floated width badly behaved for a smooth
    gradient-based solver. Holding it fixed keeps the refined parameters
    (poly coefficients + depth) smooth, so an analytic Jacobian applies and
    the fit converges in a handful of iterations instead of hundreds.
    """
    n = len(t)
    if n < 6:
        return None

    sigma0 = efrac * np.abs(y)
    sigma0[sigma0 <= 0] = efrac

    best_p0, best_width, best_cost = None, None, np.inf
    for width in width_grid:
        inbox = np.abs(t) < width / 2.0
        outbox = ~inbox
        if outbox.sum() < 3:
            continue
        a_poly = _fit_quad(t[outbox], y[outbox], 1.0 / sigma0[outbox])
        poly_full = a_poly[0] + a_poly[1] * t + a_poly[2] * t**2
        if inbox.sum() > 0 and np.all(poly_full[inbox] != 0):
            depth0 = 1.0 - np.median(y[inbox] / poly_full[inbox])
        else:
            depth0 = 0.0
        depth0 = float(np.clip(depth0, *depth_bounds))
        box = np.where(inbox, 1.0 - depth0, 1.0)
        model = poly_full * box
        cost = np.sum(((y - model) / sigma0) ** 2)
        if cost < best_cost:
            best_cost = cost
            best_width = width
            best_p0 = np.array([a_poly[0], a_poly[1], a_poly[2], depth0])

    if best_p0 is None:
        return None

    inbox_full = np.abs(t) < best_width / 2.0
    lb = [-np.inf, -np.inf, -np.inf, depth_bounds[0]]
    ub = [np.inf, np.inf, np.inf, depth_bounds[1]]

    mask = np.ones(n, dtype=bool)
    params = best_p0
    prev_n_bad = -1
    for _ in range(max_outlier_iter):
        tt, yy, ss, ib = t[mask], y[mask], sigma0[mask], inbox_full[mask]
        if len(tt) < 5:
            break
        result = least_squares(
            _notch_residuals_fixed_width,
            params,
            jac=_notch_jac_fixed_width,
            bounds=(lb, ub),
            args=(tt, yy, ss, ib),
            max_nfev=100,
        )
        params = result.x
        resid = y - _notch_model_fixed_width(params, t, inbox_full)
        rms = np.sqrt(np.mean(resid[mask] ** 2))
        if not np.isfinite(rms) or rms == 0:
            break
        outliers = np.abs(resid / rms) > sigma_clip
        n_bad = int(outliers.sum())
        if n_bad == prev_n_bad:
            mask = ~outliers
            break
        mask = ~outliers
        prev_n_bad = n_bad

    n_used = int(mask.sum())
    if n_used < 5:
        return None

    model_final = _notch_model_fixed_width(params, t, inbox_full)
    chi2_full = float(np.sum(((y[mask] - model_final[mask]) / sigma0[mask]) ** 2))
    null_coeffs = _fit_quad(t[mask], y[mask], 1.0 / sigma0[mask])
    null_model = null_coeffs[0] + null_coeffs[1] * t + null_coeffs[2] * t**2
    chi2_null = float(np.sum(((y[mask] - null_model[mask]) / sigma0[mask]) ** 2))

    return {
        "poly_coeffs": params[:3],
        "null_coeffs": null_coeffs,
        "chi2_full": chi2_full,
        "chi2_null": chi2_null,
        "n_used": n_used,
        "depth": float(params[3]),
    }


def notch_flatten(
    time: np.ndarray,
    flux: np.ndarray,
    window_length: float = 0.5,
    min_delta_bic: float = -1.0,
    efrac: float = 1e-3,
    sigma_clip: float = 3.5,
    resolvable_trans: bool = False,
    depth_bounds: tuple[float, float] = (0.0, 1.0),
    max_outlier_iter: int = 5,
) -> tuple[np.ndarray, np.ndarray]:
    """Flatten a light curve with the Notch filter.

    At each cadence, fits a quadratic trend times a multiplicative "notch
    box" within a sliding window of ``window_length`` days, compares it via
    BIC to a plain quadratic (``min_delta_bic`` sets how much evidence a
    transit-like dip needs to win), and divides by whichever polynomial
    component was selected -- so real transits are preserved rather than
    divided out.

    Parameters
    ----------
    time, flux : array_like
        Light curve arrays (flux should be roughly normalized to 1).
    window_length : float
        Sliding window width in days. Default 0.5, matching the upstream
        Notch default for short-cadence TESS/K2 data.
    min_delta_bic : float
        BIC threshold to prefer the notch model over the null (poly-only)
        model. Set to ``-np.inf`` to always keep the notch model, or
        ``np.inf`` to always use the plain quadratic.
    efrac : float
        Starting fractional flux uncertainty used to weight the fits.
    sigma_clip : float
        Outlier-rejection threshold (in local RMS units) applied inside
        each window fit.
    resolvable_trans : bool
        If True, drop the shortest (45 min) candidate transit duration from
        the search grid.
    depth_bounds : tuple
        Bounds (fraction) for the fitted notch depth.
    max_outlier_iter : int
        Maximum sigma-clip/refit iterations per window.

    Returns
    -------
    flatten_flux, trend_flux : np.ndarray
        Same convention as ``wotan.flatten``: ``flatten_flux = flux / trend_flux``.
    """
    t_arr = np.asarray(time, dtype=float)
    y_arr = np.asarray(flux, dtype=float)
    n = len(t_arr)

    if resolvable_trans:
        width_grid = np.array([1.0, 2.0, 4.0]) / 24.0
    else:
        width_grid = np.array([0.75, 1.0, 2.0, 4.0]) / 24.0

    trend = np.full(n, np.nan)
    depth_history = np.zeros(n)

    for i in range(n):
        t_i = t_arr[i]
        lo = np.searchsorted(t_arr, t_i - window_length / 2.0, side="left")
        hi = np.searchsorted(t_arr, t_i + window_length / 2.0, side="right")
        idx = np.arange(lo, hi)
        tw = t_arr[idx] - t_i
        yw = y_arr[idx]

        if len(idx) >= 3:
            line = np.polyval(np.polyfit(tw, yw, 1), tw)
            rms_line = np.sqrt(np.nanmean((yw - line) ** 2))
            rms_line = rms_line if rms_line > 0 else 1.0
            flare = (np.abs((yw - line) / rms_line) > 8.0) | (yw <= 0)
        else:
            flare = np.zeros(len(idx), dtype=bool)

        # Protect points with a previously-detected large notch depth from
        # being masked as flares/outliers in a later, overlapping window --
        # a simplified stand-in for upstream's cross-window transit guard.
        protect = np.zeros(len(idx), dtype=bool)
        if i >= 10:
            hist_std = np.std(depth_history[:i])
            if hist_std > 0:
                in_first_third = (idx - idx[0]) < (len(idx) / 3.0)
                protect = (depth_history[idx] > 5.0 * hist_std) & in_first_third

        keep = ~flare & ~protect
        tw_g, yw_g = tw[keep], yw[keep]

        fit = (
            _fit_notch_point(
                tw_g, yw_g, width_grid, efrac, sigma_clip, depth_bounds, max_outlier_iter
            )
            if keep.sum() >= 6
            else None
        )

        if fit is None:
            if keep.sum() >= 3:
                sigma = efrac * np.abs(yw_g)
                sigma[sigma <= 0] = efrac
                trend[i] = _fit_quad(tw_g, yw_g, 1.0 / sigma)[0]
            elif len(yw):
                trend[i] = np.nanmedian(yw)
            else:
                trend[i] = np.nanmedian(y_arr)
            continue

        depth_history[i] = fit["depth"]
        n_used = max(fit["n_used"], 2)
        bic_null = fit["chi2_null"] + 3 * np.log(n_used)
        bic_full = fit["chi2_full"] + 5 * np.log(n_used)
        dbic = bic_null - bic_full
        trend[i] = fit["null_coeffs"][0] if dbic < min_delta_bic else fit["poly_coeffs"][0]

    nanmask = np.isnan(trend)
    if nanmask.any() and not nanmask.all():
        trend[nanmask] = np.interp(t_arr[nanmask], t_arr[~nanmask], trend[~nanmask])
    trend[np.isnan(trend)] = np.nanmedian(y_arr)

    flatten_flux = y_arr / trend
    return flatten_flux, trend


# ---------------------------------------------------------------------------
# LOCoR
# ---------------------------------------------------------------------------
def _label_cycles(phase: np.ndarray) -> np.ndarray:
    """Label each rotation cycle by where phase wraps from ~1 back to ~0.

    Upstream's ``np.roll``-based wrap detector never assigns the final
    partial cycle its own label (it silently falls into cycle 0); this uses
    ``np.diff`` so every cycle, including a trailing partial one, gets its
    own id.
    """
    wraps = np.where(np.diff(phase) < 0)[0] + 1
    bounds = np.concatenate(([0], wraps, [len(phase)]))
    labels = np.zeros(len(phase), dtype=int)
    for k in range(len(bounds) - 1):
        labels[bounds[k] : bounds[k + 1]] = k
    return labels


def _lcbin_robust_goodidx(phase: np.ndarray, y: np.ndarray, nbins: int = 20) -> np.ndarray:
    good = np.zeros(len(phase), dtype=bool)
    edges = np.linspace(0.0, 1.0, nbins + 1)
    for b in range(nbins):
        w = np.where((phase >= edges[b]) & (phase < edges[b + 1]))[0]
        if len(w) == 0:
            continue
        _, goodind = _robust_mean(y[w], 3)
        good[w[goodind]] = True
    return good


def locor_flatten(
    time: np.ndarray,
    flux: np.ndarray,
    period: float,
    alias_num: float = 0.01,
    sigma_clip: float = 3.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Flatten a light curve with LOCoR (Locally Optimized Combination of Rotations).

    Phase-folds the light curve on ``period``, then models each rotation
    cycle as an optimal linear combination of its nearest neighboring
    cycles (fit by linear least squares), iteratively rejecting outliers.
    This tracks starspot evolution without assuming a fixed periodic shape.

    Parameters
    ----------
    time, flux : array_like
        Light curve arrays (flux should be roughly normalized to 1).
    period : float
        Stellar rotation period in days (e.g. from a GLS/Lomb-Scargle
        periodogram).
    alias_num : float
        If ``period`` is shorter than this, it is multiplied up to the
        smallest integer multiple that exceeds it, so each cycle has enough
        points for a well-conditioned fit. Default 0.01 effectively
        disables aliasing for any realistic stellar rotation period,
        matching upstream's ``run_locor`` wrapper default.
    sigma_clip : float
        Outlier-rejection threshold (MAD units) in the per-cycle fit.

    Returns
    -------
    flatten_flux, trend_flux : np.ndarray
        Same convention as ``wotan.flatten``: ``flatten_flux = flux / trend_flux``.
    """
    t = np.asarray(time, dtype=float)
    y = np.asarray(flux, dtype=float)
    n = len(t)

    if not np.isfinite(period) or period <= 0:
        raise ValueError("locor_flatten requires a positive, finite rotation `period`.")

    phase0 = _phase(t, period, t[0])
    good_global = _lcbin_robust_goodidx(phase0, y, nbins=20)

    newprot = period
    if period < alias_num:
        pfactor = int(np.ceil(alias_num / period))
        newprot = period * pfactor

    phase = _phase(t, newprot, t[0])
    cycle_id = _label_cycles(phase)
    cycles = np.unique(cycle_id)
    ncyc = len(cycles)

    flag = np.where(good_global, 1, 0).astype(int)
    flag[y <= 0] = 2

    cyc_idx, cyc_phase, cyc_flux, cyc_flag, cyc_median = [], [], [], [], []
    for c in cycles:
        idx = np.where(cycle_id == c)[0]
        med = np.nanmedian(y[idx])
        med = med if (np.isfinite(med) and med != 0) else 1.0
        cyc_idx.append(idx)
        cyc_phase.append(phase[idx])
        cyc_flux.append(y[idx] / med)
        cyc_flag.append(flag[idx].copy())
        cyc_median.append(med)

    trend = np.ones(n)

    for i in range(ncyc):
        idx = cyc_idx[i]
        this_phase = cyc_phase[i]
        this_flux = cyc_flux[i]
        npts = len(this_flux)
        init_flag = cyc_flag[i]

        use = init_flag == 1
        excluded = init_flag == 0
        if excluded.sum() > npts / 2:
            use = use | excluded  # over-aggressive coarse pass: keep them all

        n_use0 = int(use.sum())
        others = [j for j in range(ncyc) if j != i]
        size_lim = n_use0 - 3
        others = sorted(others, key=lambda j: abs(j - i))[: max(size_lim, 0)]

        ref_curves = []
        for j in others:
            keep_ref = cyc_flag[j] == 1
            if (
                keep_ref.sum() > 5
                and len(cyc_flag[j]) > len(init_flag) * 2.0 / 3.0
                and (np.max(cyc_flux[j][keep_ref]) - np.min(cyc_flux[j][keep_ref])) > 1e-14
            ):
                ref_curves.append(
                    np.interp(this_phase, cyc_phase[j][keep_ref], cyc_flux[j][keep_ref])
                )

        model = None
        if size_lim >= 1 and ref_curves and (np.max(this_flux) - np.min(this_flux)) > 1e-14:
            ref = np.array(ref_curves)
            work = use.copy()
            for _ in range(3):
                w = work.astype(float)
                amat = (ref * w) @ ref.T
                bvec = (ref * w) @ this_flux
                try:
                    coeffs = np.linalg.solve(amat, bvec)
                except np.linalg.LinAlgError:
                    coeffs, *_ = np.linalg.lstsq(amat, bvec, rcond=None)
                model = coeffs @ ref
                used_idx = np.where(work)[0]
                if len(used_idx) == 0:
                    break
                mad = np.median(np.abs(this_flux[used_idx] / model[used_idx] - 1.0)) / 0.6745
                mad = mad if mad > 0 else 1e-12
                offset = np.zeros(npts)
                offset[used_idx] = (this_flux[used_idx] / model[used_idx] - 1.0) / mad
                outliers = np.abs(offset) > sigma_clip
                if not outliers.any():
                    break
                work[outliers] = False
                if work.sum() <= ref.shape[0]:
                    break

            nofit = ~work
            if nofit.any() and work.any():
                model[nofit] = np.interp(this_phase[nofit], this_phase[work], model[work])

        if model is None:
            model = np.ones(npts)

        trend[idx] = model * cyc_median[i]

    flatten_flux = y / trend
    return flatten_flux, trend
