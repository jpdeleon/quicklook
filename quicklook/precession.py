"""Search for nodal precession as a drift in transit shape across TESS sectors (#63).

Nodal precession changes the sky-projected impact parameter ``b`` over years.
The transit gets shorter and shallower as ``b`` grows and longer as it shrinks.
With two or more sectors far enough apart, this shows up as a secular trend in
per-sector ``b`` and ``T14`` (Watanabe et al. 2026, TOI-1355 b,
arXiv:2608.17387; Zasche et al. 2023, V907 Sco, arXiv:2302.04613).

Steps:

1. ``extract_sector_transits`` keeps the data around each transit, split by
   sector, using the existing (TLS or custom) ephemeris.
2. ``fit_sector_transits`` does one joint ``batman`` fit: shared ``a/Rs`` and
   limb darkening, with ``Rp/Rs``, ``b``, a mid-time offset and a baseline free
   per sector. Errors come from the Jacobian covariance, scaled up by the reduced
   chi^2. Derived quantities (``T14``, depth) come from draws of that covariance.
3. ``fit_precession_model`` fits static, linear and (with >= 4 sectors)
   quadratic ``b(t)`` and compares them by BIC.
4. ``assess_precession`` applies the candidate test and the systematics flags.

Geometry (circular orbit): ``b = (a/Rs) cos i``, so
``db/dt = -(a/Rs) sin(i) di/dt`` and ``di/dt = -(db/dt) / ((a/Rs) sin i)``.
Only ``|b|`` is measurable, so a trend through ``b = 0`` looks like a V.

Planet branch only. The eclipsing-binary branch (``ellc``/PHOEBE, triple-star
precession period) is not implemented yet.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional, Sequence

import batman
import numpy as np
from scipy.optimize import least_squares

DAYS_PER_YEAR = 365.25
MIN_SECTORS = 2
MIN_BASELINE_YR = 1.0  # below this db/dt is reported but never a candidate
SIGNIFICANCE_MIN = 3.0  # |db/dt| / sigma
DELTA_BIC_MIN = 2.0  # BIC(static) - BIC(linear): positive evidence for a trend
GRAZING_CORR_MAX = 0.95  # |corr(b, Rp/Rs)| above this: b is not constrained on its own
CROWDSAP_SPREAD_MAX = 0.01  # crowding ratio range across sectors worth flagging
CADENCE_MISMATCH_TOL = 0.1  # relative cadence difference counted as "mixed"
T0_OFFSET_SIGMA = 3.0  # per-sector mid-time shift flagged as TTV-like
LONG_CADENCE = 10.0 / 1440  # days; supersample the model above this cadence
SUPERSAMPLE = 7
HALF_WINDOW_DURATIONS = 1.5  # keep |phase| < this * T14 around each transit
# TLS durations can be far too short (TOI-1355 b: 14 min), which would leave no
# baseline or transit wings in the window; never keep less than this (days),
# and never more than a quarter period.
MIN_HALF_WINDOW = 0.25
MIN_DT0_BOUND = 0.02  # days; mid-time offset bound when the duration is tiny
MIN_IN_TRANSIT_POINTS = 5
N_DRAWS = 2000
DEFAULT_LIMB_DARKENING = (0.35, 0.25)  # quadratic; shared, so it cannot fake a trend
B_STARTS = (0.2, 0.5, 0.8)  # b is multi-modal with a/Rs; try a few starts
B_BOUNDARY = 0.05  # b fitted at the b=0 bound: covariance errors are meaningless
AU_PER_RSUN = 215.032
A_RS_PRIOR_FLOOR = 0.1  # relative a/Rs prior width when R*, M* errors are missing
BLOCKING_FLAGS = (
    "short_baseline",
    "rp_absorbs_b",
    "shape_sign_mismatch",
    "b_at_boundary",
)


@dataclass(frozen=True)
class SectorSpan:
    """Time range and dilution keywords of one stitched sector."""

    sector: int
    tmin: float
    tmax: float
    crowdsap: float = np.nan
    flfrcsap: float = np.nan


@dataclass(frozen=True)
class SectorData:
    """Data within a few durations of each transit in one sector."""

    span: SectorSpan
    time: np.ndarray
    flux: np.ndarray
    flux_err: np.ndarray
    cadence: float
    n_transits: int
    t_mid: float  # median time of in-transit points


@dataclass(frozen=True)
class SectorTransitFit:
    sector: int
    t_mid: float
    n_transits: int
    cadence: float
    crowdsap: float
    flfrcsap: float
    b: float
    b_err: float
    rp_rs: float
    rp_rs_err: float
    depth: float  # (Rp/Rs)^2
    depth_err: float
    model_depth: float  # deepest point of the best-fit model
    t14: float  # days
    t14_err: float
    dt0: float  # mid-time offset from the ephemeris, days
    dt0_err: float
    corr_b_rp: float


@dataclass(frozen=True)
class TrendFit:
    degree: int
    coeffs: tuple  # lowest order first, per year since t_ref
    errors: tuple
    chi2: float
    bic: float
    dof: int


@dataclass(frozen=True)
class PrecessionModels:
    t_ref: float  # days (same time system as the light curve)
    static: TrendFit
    linear: TrendFit
    quadratic: Optional[TrendFit]
    delta_bic: float  # static - linear (> 0 favours a trend)
    delta_bic_quadratic: float  # linear - quadratic (nan without a quadratic fit)


@dataclass(frozen=True)
class PrecessionResult:
    sector_fits: tuple
    models: Optional[PrecessionModels]
    a_rs: float = np.nan
    a_rs_err: float = np.nan
    db_dt: float = np.nan  # per year
    db_dt_err: float = np.nan
    significance: float = np.nan
    delta_bic: float = np.nan
    dt14_dt: float = np.nan  # days per year
    dt14_dt_err: float = np.nan
    drp_dt: float = np.nan  # per year
    drp_dt_err: float = np.nan
    baseline_yr: float = 0.0
    di_dt: float = np.nan  # deg per year
    years_to_exit: float = np.nan  # from the last sector until b = 1 + Rp/Rs
    is_candidate: bool = False
    flags: tuple = ()
    reason: str = ""


# --- data -------------------------------------------------------------------


def _meta_float(meta: Any, key: str) -> float:
    try:
        return float(meta.get(key, np.nan))
    except (TypeError, ValueError):
        return np.nan


def sector_span_from_lc(lc: Any) -> Optional[SectorSpan]:
    """``SectorSpan`` of a single-sector lightkurve light curve, else None."""
    sector = lc.meta.get("SECTOR")
    time = np.asarray(lc.time.value, dtype=float)
    time = time[np.isfinite(time)]
    if sector is None or time.size == 0:
        return None
    return SectorSpan(
        sector=int(sector),
        tmin=float(time.min()),
        tmax=float(time.max()),
        crowdsap=_meta_float(lc.meta, "CROWDSAP"),
        flfrcsap=_meta_float(lc.meta, "FLFRCSAP"),
    )


def correct_dilution(flux: np.ndarray, crowdsap: float) -> np.ndarray:
    """Apply the SPOC crowding and flux-fraction correction to normalised SAP flux.

    PDCSAP already has this applied; use it for the SAP cross-check. After
    removing the ``(1 - CROWDSAP)`` contaminating fraction and renormalising,
    ``FLFRCSAP`` (a pure scale) cancels, so only ``CROWDSAP`` changes the
    depth: ``f' = (f / median - (1 - c)) / c``. Returns the normalised flux
    unchanged when the keyword is missing or invalid.
    """
    flux = np.asarray(flux, dtype=float)
    norm = flux / np.nanmedian(flux)
    if not np.isfinite(crowdsap) or not 0 < crowdsap <= 1:
        return norm
    return (norm - (1 - crowdsap)) / crowdsap


def _phase(time: np.ndarray, period: float, t0: float) -> np.ndarray:
    return (time - t0 + 0.5 * period) % period - 0.5 * period


def _robust_std(x: np.ndarray) -> float:
    return float(1.4826 * np.nanmedian(np.abs(x - np.nanmedian(x))))


def extract_sector_transits(
    time,
    flux,
    flux_err,
    spans: Sequence[SectorSpan],
    period: float,
    t0: float,
    duration: float,
    half_window: float = HALF_WINDOW_DURATIONS,
) -> tuple:
    """Split a stitched light curve into per-sector transit windows.

    Sectors with fewer than ``MIN_IN_TRANSIT_POINTS`` in-transit cadences are
    dropped. Missing or non-finite errors are replaced by the robust scatter
    out of transit.
    """
    time = np.asarray(time, dtype=float)
    flux = np.asarray(flux, dtype=float)
    err = np.full_like(flux, np.nan) if flux_err is None else np.asarray(flux_err, dtype=float)
    phase = _phase(time, period, t0)
    half = min(max(half_window * duration, MIN_HALF_WINDOW), 0.25 * period)
    near = np.abs(phase) < half
    in_transit = np.abs(phase) < 0.5 * duration
    finite = np.isfinite(time) & np.isfinite(flux)
    out = []
    for span in spans:
        sel = finite & (time >= span.tmin) & (time <= span.tmax)
        if np.count_nonzero(sel & in_transit) < MIN_IN_TRANSIT_POINTS:
            continue
        t, f, e = time[sel & near], flux[sel & near], err[sel & near]
        oot = ~in_transit[sel & near]
        scatter = _robust_std(f[oot]) if np.count_nonzero(oot) > 2 else _robust_std(f)
        e = np.where(np.isfinite(e) & (e > 0), e, scatter)
        intr = sel & in_transit
        epochs = np.round((time[intr] - t0) / period)
        out.append(
            SectorData(
                span=span,
                time=t,
                flux=f,
                flux_err=e,
                cadence=float(np.median(np.diff(time[sel]))),
                n_transits=int(np.unique(epochs).size),
                t_mid=float(np.median(time[intr])),
            )
        )
    return tuple(out)


# --- per-sector transit fits ---------------------------------------------------


def transit_duration(period, a_rs, rp_rs, b):
    """Total duration T14 in days for a circular orbit (Winn 2010, eq. 14)."""
    a_rs, rp_rs, b = (np.asarray(x, dtype=float) for x in (a_rs, rp_rs, b))
    with np.errstate(invalid="ignore", divide="ignore"):
        sin_i = np.sqrt(np.clip(1 - (b / a_rs) ** 2, 0, 1))
        chord = np.sqrt(np.clip((1 + rp_rs) ** 2 - b**2, 0, None))
        return period / np.pi * np.arcsin(np.clip(chord / (a_rs * sin_i), 0, 1))


def _make_models(sectors, period: float, t0: float, ld) -> list:
    models = []
    for s in sectors:
        params = batman.TransitParams()
        params.t0, params.per, params.rp, params.a = t0, period, 0.1, 10.0
        params.inc, params.ecc, params.w = 89.0, 0.0, 90.0
        params.limb_dark, params.u = "quadratic", list(ld)
        kwargs = (
            {"supersample_factor": SUPERSAMPLE, "exp_time": s.cadence}
            if s.cadence > LONG_CADENCE
            else {}
        )
        models.append((params, batman.TransitModel(params, s.time, **kwargs)))
    return models


def _light_curve(params, model, t0, a_rs, dt0, rp, b) -> np.ndarray:
    params.t0, params.rp, params.a = t0 + dt0, rp, a_rs
    params.inc = float(np.degrees(np.arccos(np.clip(b / a_rs, 0, 1))))
    return model.light_curve(params)


def _residuals(theta, sectors, models, t0, a_rs_prior=None) -> np.ndarray:
    a_rs = theta[0]
    chunks = []
    for k, (s, (params, model)) in enumerate(zip(sectors, models)):
        dt0, rp, b, c0 = theta[1 + 4 * k : 5 + 4 * k]
        lc = c0 * _light_curve(params, model, t0, a_rs, dt0, rp, b)
        chunks.append((s.flux - lc) / s.flux_err)
    if a_rs_prior is not None:
        chunks.append(np.array([(a_rs - a_rs_prior[0]) / a_rs_prior[1]]))
    return np.concatenate(chunks)


def a_rs_from_star(
    period: float, r_star: float, m_star: float, r_err: float = np.nan, m_err: float = np.nan
) -> Optional[tuple]:
    """``(a/Rs, sigma)`` from Kepler's third law and R*, M* in solar units, else None.

    ``a/Rs = (G M* P^2 / 4 pi^2)^(1/3) / R*``; the error adds ``sigma_M/(3M)``
    and ``sigma_R/R`` in quadrature (a 10% floor when either error is missing).
    Breaks the ``a/Rs``-``b`` degeneracy that a free ``a/Rs`` leaves at
    moderate ``b``.
    """
    if not all(np.isfinite(x) and x > 0 for x in (period, r_star, m_star)):
        return None
    a_au = (m_star * (period / DAYS_PER_YEAR) ** 2) ** (1 / 3)
    a_rs = a_au * AU_PER_RSUN / r_star
    if np.isfinite(r_err) and np.isfinite(m_err) and r_err > 0 and m_err > 0:
        rel = np.hypot(m_err / (3 * m_star), r_err / r_star)
    else:
        rel = A_RS_PRIOR_FLOOR
    return float(a_rs), float(a_rs * max(rel, 1e-3))


def _initial_guess(sectors, period, t0, duration, b0) -> tuple:
    x0, lo, hi = [], [], []
    depths = []
    for s in sectors:
        intr = np.abs(_phase(s.time, period, t0)) < 0.25 * duration
        base = np.median(s.flux[~intr]) if np.any(~intr) else 1.0
        depth = 1 - np.median(s.flux[intr]) / base if np.any(intr) else 1e-3
        depths.append(max(depth, 1e-5))
        x0 += [0.0, np.sqrt(depths[-1]), b0, base]
        dt0_bound = max(0.5 * duration, MIN_DT0_BOUND)
        lo += [-dt0_bound, 1e-3, 0.0, 0.5]
        hi += [dt0_bound, 0.9, 1.5, 1.5]
    k = np.sqrt(np.median(depths))
    a_rs = period / (np.pi * duration) * np.sqrt(max((1 + k) ** 2 - b0**2, 1e-3))
    a_rs = float(np.clip(a_rs, 2.0, 500.0))
    return np.array([a_rs] + x0), np.array([1.5] + lo), np.array([1000.0] + hi)


def _best_fit(sectors, models, period, t0, duration, a_rs_prior=None):
    best = None
    for b0 in B_STARTS:
        x0, lo, hi = _initial_guess(sectors, period, t0, duration, b0)
        if a_rs_prior is not None:
            x0[0] = a_rs_prior[0]
        x0 = np.clip(x0, lo + 1e-9, hi - 1e-9)
        res = least_squares(
            _residuals,
            x0,
            bounds=(lo, hi),
            args=(sectors, models, t0, a_rs_prior),
            x_scale="jac",
        )
        if best is None or res.cost < best.cost:
            best = res
    return best


def _covariance(res) -> np.ndarray:
    jac = res.jac
    dof = max(jac.shape[0] - jac.shape[1], 1)
    redchi2 = 2 * res.cost / dof
    return np.linalg.pinv(jac.T @ jac) * max(redchi2, 1.0)


def _spread(draws: np.ndarray) -> float:
    lo, hi = np.nanpercentile(draws, [15.865, 84.135])
    return float(0.5 * (hi - lo))


def fit_sector_transits(
    sectors: Sequence[SectorData],
    period: float,
    t0: float,
    duration: float,
    limb_darkening=DEFAULT_LIMB_DARKENING,
    a_rs_prior: Optional[tuple] = None,
    seed: int = 0,
) -> tuple:
    """Joint transit fit; returns ``(fits, a_rs, a_rs_err)``.

    ``a/Rs`` and limb darkening are shared (they are constant in time);
    ``Rp/Rs``, ``b``, mid-time offset and baseline are free per sector so a
    shape change is not forced into any one parameter. ``a_rs_prior`` is an
    optional Gaussian ``(mean, sigma)`` (see ``a_rs_from_star``).
    """
    models = _make_models(sectors, period, t0, limb_darkening)
    res = _best_fit(sectors, models, period, t0, duration, a_rs_prior)
    cov = _covariance(res)
    err = np.sqrt(np.clip(np.diag(cov), 0, None))
    rng = np.random.default_rng(seed)
    draws = rng.multivariate_normal(res.x, cov, N_DRAWS, check_valid="ignore")
    a_draws = draws[:, 0]
    fits = []
    for k, (s, (params, _model)) in enumerate(zip(sectors, models)):
        i_dt0, i_rp, i_b = 1 + 4 * k, 2 + 4 * k, 3 + 4 * k
        a_rs, dt0, rp, b = res.x[0], res.x[i_dt0], res.x[i_rp], res.x[i_b]
        denom = err[i_b] * err[i_rp]
        corr = float(cov[i_b, i_rp] / denom) if denom > 0 else np.nan
        params.u = list(limb_darkening)
        fine = np.linspace(-0.01, 0.01, 21) + t0 + dt0
        dense = batman.TransitModel(params, fine)
        model_depth = 1 - float(np.min(_light_curve(params, dense, t0, a_rs, dt0, rp, b)))
        fits.append(
            SectorTransitFit(
                sector=s.span.sector,
                t_mid=s.t_mid,
                n_transits=s.n_transits,
                cadence=s.cadence,
                crowdsap=s.span.crowdsap,
                flfrcsap=s.span.flfrcsap,
                b=float(b),
                b_err=float(err[i_b]),
                rp_rs=float(rp),
                rp_rs_err=float(err[i_rp]),
                depth=float(rp**2),
                depth_err=_spread(draws[:, i_rp] ** 2),
                model_depth=model_depth,
                t14=float(transit_duration(period, a_rs, rp, b)),
                t14_err=_spread(
                    transit_duration(period, a_draws, draws[:, i_rp], np.abs(draws[:, i_b]))
                ),
                dt0=float(dt0),
                dt0_err=float(err[i_dt0]),
                corr_b_rp=corr,
            )
        )
    return tuple(fits), float(res.x[0]), float(err[0])


# --- trend models -------------------------------------------------------------


def fit_trend(t_yr, y, yerr, degree: int) -> TrendFit:
    """Weighted polynomial fit; errors scaled up when reduced chi^2 > 1."""
    t_yr, y, yerr = (np.asarray(x, dtype=float) for x in (t_yr, y, yerr))
    design = np.vander(t_yr, degree + 1, increasing=True) / yerr[:, None]
    yw = y / yerr
    coeffs, *_ = np.linalg.lstsq(design, yw, rcond=None)
    cov = np.linalg.pinv(design.T @ design)
    chi2 = float(np.sum((yw - design @ coeffs) ** 2))
    n, k = y.size, degree + 1
    dof = n - k
    if dof > 0 and chi2 / dof > 1:
        cov = cov * (chi2 / dof)
    return TrendFit(
        degree=degree,
        coeffs=tuple(float(c) for c in coeffs),
        errors=tuple(float(e) for e in np.sqrt(np.clip(np.diag(cov), 0, None))),
        chi2=chi2,
        bic=chi2 + k * np.log(n),
        dof=dof,
    )


def fit_precession_model(t_mid, b, b_err, t_ref: Optional[float] = None) -> PrecessionModels:
    """Static, linear and (>= 4 points) quadratic fits of ``b`` against time."""
    t_mid = np.asarray(t_mid, dtype=float)
    t_ref = float(np.mean(t_mid)) if t_ref is None else float(t_ref)
    t_yr = (t_mid - t_ref) / DAYS_PER_YEAR
    static = fit_trend(t_yr, b, b_err, 0)
    linear = fit_trend(t_yr, b, b_err, 1)
    quadratic = fit_trend(t_yr, b, b_err, 2) if t_mid.size >= 4 else None
    return PrecessionModels(
        t_ref=t_ref,
        static=static,
        linear=linear,
        quadratic=quadratic,
        delta_bic=static.bic - linear.bic,
        delta_bic_quadratic=linear.bic - quadratic.bic if quadratic else np.nan,
    )


# --- derived quantities -------------------------------------------------------


def inclination_rate(db_dt: float, b: float, a_rs: float) -> float:
    """``di/dt`` in deg/yr from ``db/dt`` (per yr): ``-(db/dt) / ((a/Rs) sin i)``."""
    sin_i = np.sqrt(max(1 - (b / a_rs) ** 2, 0.0))
    if sin_i == 0 or not np.isfinite(a_rs):
        return np.nan
    return float(np.degrees(-db_dt / (a_rs * sin_i)))


def years_until_no_transit(b_now: float, rp_rs: float, db_dt: float) -> float:
    """Years until ``b`` reaches ``1 + Rp/Rs`` (no transit), nan if ``b`` is not growing."""
    if not (np.isfinite(db_dt) and db_dt > 0):
        return np.nan
    return float(max(1 + rp_rs - b_now, 0.0) / db_dt)


# --- verdict ------------------------------------------------------------------


def _slope(fits, attr: str, t_ref: float) -> tuple:
    t_yr = (np.array([f.t_mid for f in fits]) - t_ref) / DAYS_PER_YEAR
    y = np.array([getattr(f, attr) for f in fits])
    yerr = np.array([getattr(f, attr + "_err") for f in fits])
    if not np.all(np.isfinite(yerr) & (yerr > 0)):
        return np.nan, np.nan
    trend = fit_trend(t_yr, y, yerr, 1)
    return trend.coeffs[1], trend.errors[1]


def systematics_flags(fits: Sequence[SectorTransitFit]) -> tuple:
    """Flags that do not depend on the trend fit."""
    flags = []
    if any(abs(f.corr_b_rp) > GRAZING_CORR_MAX for f in fits if np.isfinite(f.corr_b_rp)):
        flags.append("grazing_degeneracy")
    if any(f.b < B_BOUNDARY for f in fits):
        flags.append("b_at_boundary")
    crowd = np.array([f.crowdsap for f in fits], dtype=float)
    crowd = crowd[np.isfinite(crowd)]
    if crowd.size >= 2 and np.ptp(crowd) > CROWDSAP_SPREAD_MAX:
        flags.append("contamination_varies")
    cadences = np.array([f.cadence for f in fits], dtype=float)
    if np.ptp(cadences) > CADENCE_MISMATCH_TOL * np.min(cadences):
        flags.append("mixed_cadence")
    if any(abs(f.dt0) > T0_OFFSET_SIGMA * f.dt0_err for f in fits if f.dt0_err > 0):
        flags.append("t0_offset")
    return tuple(flags)


def assess_precession(
    fits: Sequence[SectorTransitFit],
    a_rs: float = np.nan,
    a_rs_err: float = np.nan,
    min_baseline_yr: float = MIN_BASELINE_YR,
) -> PrecessionResult:
    """Trend fits, flags and the candidate decision for per-sector fits.

    A candidate needs ``|db/dt| >= 3 sigma``, ``BIC(static) - BIC(linear) >= 2``,
    and none of: a baseline under ``min_baseline_yr``, a ``b``-``Rp/Rs``
    degeneracy together with an ``Rp/Rs`` trend, ``b`` pinned at 0, or a ``T14`` trend with the same sign as ``db/dt`` (``T14``
    must shrink as ``b`` grows).
    """
    fits = tuple(sorted(fits, key=lambda f: f.t_mid))
    if len(fits) < MIN_SECTORS:
        return PrecessionResult(
            sector_fits=fits, models=None, reason=f"needs >= {MIN_SECTORS} sectors with transits"
        )
    t_mid = np.array([f.t_mid for f in fits])
    b = np.array([f.b for f in fits])
    b_err = np.array([f.b_err for f in fits])
    if not np.all(np.isfinite(b_err) & (b_err > 0)):
        return PrecessionResult(sector_fits=fits, models=None, reason="unusable b uncertainties")
    models = fit_precession_model(t_mid, b, b_err)
    db_dt, db_dt_err = models.linear.coeffs[1], models.linear.errors[1]
    significance = abs(db_dt) / db_dt_err if db_dt_err > 0 else np.nan
    dt14_dt, dt14_dt_err = _slope(fits, "t14", models.t_ref)
    drp_dt, drp_dt_err = _slope(fits, "rp_rs", models.t_ref)
    baseline_yr = float(np.ptp(t_mid) / DAYS_PER_YEAR)

    flags = list(systematics_flags(fits))
    if baseline_yr < min_baseline_yr:
        flags.append("short_baseline")
    if np.isfinite(dt14_dt) and dt14_dt != 0 and np.sign(dt14_dt) == np.sign(db_dt):
        flags.append("shape_sign_mismatch")
    if np.isfinite(drp_dt_err) and drp_dt_err > 0 and abs(drp_dt) / drp_dt_err >= SIGNIFICANCE_MIN:
        flags.append("depth_change")
    if models.delta_bic_quadratic >= DELTA_BIC_MIN:
        flags.append("nonlinear")
    # Near b ~ 1, b and Rp/Rs are correlated; that alone is expected at high
    # b. Block only when Rp/Rs also trends, i.e. it may be absorbing the change.
    if "grazing_degeneracy" in flags and "depth_change" in flags:
        flags.append("rp_absorbs_b")

    blocking = [f for f in flags if f in BLOCKING_FLAGS]
    detected = significance >= SIGNIFICANCE_MIN and models.delta_bic >= DELTA_BIC_MIN
    if not detected:
        reason = f"no significant trend ({significance:.1f} sigma, dBIC={models.delta_bic:.1f})"
    elif blocking:
        reason = f"trend at {significance:.1f} sigma but flagged: {', '.join(blocking)}"
    else:
        reason = f"b drifts at {significance:.1f} sigma (dBIC={models.delta_bic:.1f})"
    last = fits[-1]
    b_last = models.linear.coeffs[0] + db_dt * (last.t_mid - models.t_ref) / DAYS_PER_YEAR
    return PrecessionResult(
        sector_fits=fits,
        models=models,
        a_rs=a_rs,
        a_rs_err=a_rs_err,
        db_dt=db_dt,
        db_dt_err=db_dt_err,
        significance=significance,
        delta_bic=models.delta_bic,
        dt14_dt=dt14_dt,
        dt14_dt_err=dt14_dt_err,
        drp_dt=drp_dt,
        drp_dt_err=drp_dt_err,
        baseline_yr=baseline_yr,
        di_dt=inclination_rate(db_dt, float(np.mean(b)), a_rs),
        years_to_exit=years_until_no_transit(b_last, last.rp_rs, db_dt),
        is_candidate=bool(detected and not blocking),
        flags=tuple(flags),
        reason=reason,
    )


def search_precession(
    time,
    flux,
    flux_err,
    spans: Sequence[SectorSpan],
    period: float,
    t0: float,
    duration: float,
    limb_darkening=DEFAULT_LIMB_DARKENING,
    min_baseline_yr: float = MIN_BASELINE_YR,
    a_rs_prior: Optional[tuple] = None,
    seed: int = 0,
) -> PrecessionResult:
    """Top-level search on a stitched, flattened light curve with sector spans.

    Pass ``a_rs_prior`` (from ``a_rs_from_star``) when the stellar density is
    known: at moderate ``b`` a free ``a/Rs`` lets every sector's ``b`` slide
    together and biases ``db/dt``.
    """
    sectors = extract_sector_transits(time, flux, flux_err, spans, period, t0, duration)
    if len(sectors) < MIN_SECTORS:
        return PrecessionResult(
            sector_fits=(),
            models=None,
            reason=f"{len(sectors)} sector(s) cover a transit; needs >= {MIN_SECTORS}",
        )
    fits, a_rs, a_rs_err = fit_sector_transits(
        sectors,
        period,
        t0,
        duration,
        limb_darkening=limb_darkening,
        a_rs_prior=a_rs_prior,
        seed=seed,
    )
    return assess_precession(fits, a_rs, a_rs_err, min_baseline_yr=min_baseline_yr)
