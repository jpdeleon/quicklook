"""Eclipsing-binary branch of the precession search (#63).

In a hierarchical triple the third star makes the inner binary's orbital plane
precess, so the inner inclination ``i`` drifts by up to degrees per year and
the eclipses deepen or fade (Zasche et al. 2023, V907 Sco, arXiv:2302.04613).
This branch fits ``i`` in each TESS sector and tests for a secular trend.

Model: two spherical, limb-darkened stars on a circular orbit. Each eclipse is
a ``batman`` occultation of one star by an opaque disc (exact for any radius
ratio, including an occulter larger than the star), diluted by the other
star's light and an optional third light::

    primary   (star 2 in front): F = [(1 - l3) (f1 + L) / (1 + L) + l3] c
    secondary (star 1 in front): F = [(1 - l3) (1 + L f2) / (1 + L) + l3] c

with ``L = J k^2`` (``J`` the surface-brightness ratio, ``k = R2/R1``). This
leaves out ellipsoidal variation and reflection, which the flattened light
curve has largely removed for detached systems; ``ellc`` (no wheels; needs a
Fortran build) and PHOEBE (seconds per model) were too heavy for a joint
least-squares fit over many sectors.

Shared across sectors: ``(R1+R2)/a``, ``k``, ``J`` and the secondary-eclipse
offset from phase 0.5 (``~(2P/pi) e cos w``). Free per sector: ``i``, a
mid-time offset and a baseline. Because the radii are shared, any change in
eclipse depth or duration has to come from ``i``.

The trend test is the planet branch's (static/linear/quadratic by BIC) applied
to ``i(t)`` with the joint covariance. ``b = cos(i) / (R1/a)`` and
``db/dt = -(sin i / (R1/a)) di/dt`` are reported too, so EB and planet results
share one ``db/dt`` column. A model-free check compares the measured eclipse
depth trend with ``di/dt``: as ``i`` moves toward 90 deg the eclipse deepens.

Precession period (quadrupole, Söderhjelm 1975; Zasche et al. 2023 eqs. 2-3)::

    nu / n2 = -(3 / (4 eta2^3)) (m3 / (m3 + M1)) (n2 / n1) cos J
              * sqrt(1 + gamma^2 + 2 gamma cos J)
    gamma   = (mu1 / (m3 eta2)) ((m3 + M1) / M1 * n2 / n1)^(1/3)

with ``n = 2 pi / P``, ``eta2 = sqrt(1 - e2^2)``, ``M1 = m1a + m1b``,
``mu1 = m1a m1b / M1`` and ``J`` the mutual inclination. The inner normal
precesses on a cone of half-angle ``<= J`` about the total angular momentum, so
``|d cos i / dt| <= nu sin J`` and ``|di/dt| <= nu sin J / sin i``. A measured
``di/dt`` alone therefore gives ``P_nu <= 2 pi / (sin i |di/dt|)``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
from scipy.optimize import least_squares

from quicklook.precession import (
    DAYS_PER_YEAR,
    DEFAULT_LIMB_DARKENING,
    DELTA_BIC_MIN,
    MIN_BASELINE_YR,
    MIN_DT0_BOUND,
    MIN_SECTORS,
    N_DRAWS,
    SIGNIFICANCE_MIN,
    T0_OFFSET_SIGMA,
    PrecessionModels,
    SectorData,
    SectorSpan,
    _phase,
    covariance,
    data_flags,
    extract_sector_transits,
    fit_precession_model,
    fit_trend,
    import_batman,
    spread,
    transit_duration,
    transit_model,
)

EB_B_STARTS = (0.3, 0.8, 1.3)  # EBs are often grazing: also start beyond b = 1
N_SHARED = 4  # r_sum, k, J, secondary offset
N_PER_SECTOR = 3  # mid-time offset, inclination, baseline
INC_MIN_DEG = 30.0
INC_BOUNDARY_DEG = 0.1  # i fitted within this of 90 deg: covariance errors are meaningless
R_SUM_BOUNDS = (0.005, 0.9)
K_BOUNDS = (0.1, 3.0)
SB_BOUNDS = (0.0, 5.0)
CORE_FRACTION = 0.1  # empirical depth: mean of points within this * duration of mid-eclipse
MIN_CORE_POINTS = 3
DEPTH_MISMATCH_SIGMA = 2.0  # opposite-sign depth trend this significant blocks a candidate
DILUTION_DEPTH_RATIO = 0.5  # CROWDSAP spread this fraction of the depth change blocks
EB_BLOCKING_FLAGS = (
    "short_baseline",
    "i_at_boundary",
    "depth_sign_mismatch",
    "dilution_mimics_depth",
)


@dataclass(frozen=True)
class EclipseWindows:
    """Primary-eclipse data of one sector and, when covered, its secondary eclipses."""

    primary: SectorData
    secondary: Optional[SectorData] = None


@dataclass(frozen=True)
class EclipseSectorFit:
    sector: int
    t_mid: float
    n_eclipses: int
    cadence: float
    crowdsap: float
    flfrcsap: float
    inc: float  # deg
    inc_err: float
    b: float  # cos(i) / (R1/a), in primary radii
    b_err: float
    depth: float  # primary-eclipse depth of the best-fit model (normalised flux)
    depth_err: float
    t14: float  # days
    t14_err: float
    dt0: float  # days
    dt0_err: float
    depth_data: float  # model-free primary depth near mid-eclipse (nan if too few points)
    depth_data_err: float
    has_secondary: bool


@dataclass(frozen=True)
class BinaryParams:
    r_sum: float  # (R1 + R2) / a
    r_sum_err: float
    k: float  # R2 / R1
    k_err: float
    sb_ratio: float  # J, surface brightness of star 2 over star 1
    sb_ratio_err: float
    dt_sec: float  # secondary-eclipse offset from phase 0.5, days
    dt_sec_err: float
    third_light: float


@dataclass(frozen=True)
class EBPrecessionResult:
    sector_fits: tuple
    models: Optional[PrecessionModels]
    binary: Optional[BinaryParams] = None
    di_dt: float = np.nan  # deg per year
    di_dt_err: float = np.nan
    db_dt: float = np.nan  # per year, b in primary radii
    db_dt_err: float = np.nan
    significance: float = np.nan
    delta_bic: float = np.nan
    ddepth_dt: float = np.nan  # model-free primary depth, per year
    ddepth_dt_err: float = np.nan
    baseline_yr: float = 0.0
    max_nodal_period_yr: float = np.nan  # 2 pi / (sin i |di/dt|), only for a detected trend
    is_candidate: bool = False
    flags: tuple = ()
    reason: str = ""


# --- triple-star geometry -----------------------------------------------------


def nodal_precession_period(
    p_inner: float,
    p_outer: float,
    m1a: float,
    m1b: float,
    m_third: float,
    mutual_inc_deg: float,
    e_outer: float = 0.0,
) -> float:
    """Nodal precession period (same unit as the periods) from Zasche et al. 2023 eq. 2.

    Quadrupole approximation for a circular inner orbit; valid for mutual
    inclinations below about 40 deg or above 140 deg (outside the Kozai-Lidov
    regime). Masses in any common unit.
    """
    if not 0 <= e_outer < 1:
        raise ValueError(f"e_outer must be in [0, 1), got {e_outer}")
    if min(p_inner, p_outer, m1a, m1b, m_third) <= 0:
        raise ValueError("periods and masses must be positive")
    m_in = m1a + m1b
    mu1 = m1a * m1b / m_in
    eta2 = np.sqrt(1 - e_outer**2)
    ratio = p_inner / p_outer  # n2 / n1
    cos_j = np.cos(np.radians(mutual_inc_deg))
    gamma = mu1 / (m_third * eta2) * ((m_third + m_in) / m_in * ratio) ** (1 / 3)
    nu_over_n2 = (-3 / (4 * eta2**3) * m_third / (m_third + m_in) * ratio * cos_j) * np.sqrt(
        1 + gamma**2 + 2 * gamma * cos_j
    )
    if abs(cos_j) < 1e-12:  # J = 90 deg: no quadrupole nodal precession
        return np.inf
    return float(p_outer / abs(nu_over_n2))


def max_inclination_rate(period_nodal: float, mutual_inc_deg: float, inc_deg: float) -> float:
    """Largest ``|di/dt|`` (deg per period unit) a given precession can produce.

    The inner orbit normal moves on a cone of half-angle at most ``J``, so
    ``|di/dt| <= (2 pi / P_nu) sin J / sin i``.
    """
    sin_i = np.sin(np.radians(inc_deg))
    if sin_i <= 0 or not np.isfinite(period_nodal) or period_nodal <= 0:
        return np.nan
    rate = 2 * np.pi / period_nodal * np.sin(np.radians(mutual_inc_deg)) / sin_i
    return float(np.degrees(rate))


def max_nodal_period(di_dt_deg: float, inc_deg: float) -> float:
    """Upper limit on the precession period from a measured ``di/dt`` (any J)."""
    rate = abs(np.radians(di_dt_deg)) * np.sin(np.radians(inc_deg))
    if not np.isfinite(rate) or rate == 0:
        return np.nan
    return float(2 * np.pi / rate)


# --- data -----------------------------------------------------------------------


def extract_sector_eclipses(
    time,
    flux,
    flux_err,
    spans: Sequence[SectorSpan],
    period: float,
    t0: float,
    duration: float,
) -> tuple:
    """Per-sector primary-eclipse windows, each with the secondary window when covered.

    Sectors without enough primary-eclipse cadences are dropped; a missing
    secondary window only leaves the surface-brightness ratio less constrained.
    """
    args = (time, flux, flux_err, spans, period)
    primaries = extract_sector_transits(*args, t0, duration)
    secondaries = {
        s.span.sector: s for s in extract_sector_transits(*args, t0 + 0.5 * period, duration)
    }
    return tuple(EclipseWindows(p, secondaries.get(p.span.sector)) for p in primaries)


def empirical_depth(sector: SectorData, period, t0, duration, dt0=0.0) -> tuple:
    """``(depth, error)`` from the data alone: out-of-eclipse median minus core mean."""
    phase = _phase(sector.time, period, t0 + dt0)
    core = np.abs(phase) < max(CORE_FRACTION * duration, sector.cadence)
    out = np.abs(phase) > 0.5 * duration
    if np.count_nonzero(core) < MIN_CORE_POINTS or np.count_nonzero(out) < MIN_CORE_POINTS:
        return np.nan, np.nan
    base = np.median(sector.flux[out])
    depth = (base - np.mean(sector.flux[core])) / base
    err = np.std(sector.flux[core] - np.mean(sector.flux[core]), ddof=1)
    err = max(err, np.median(sector.flux_err[core])) / np.sqrt(np.count_nonzero(core)) / base
    return float(depth), float(err)


# --- joint eclipse fit -------------------------------------------------------------


def _radii(r_sum, k):
    r1 = r_sum / (1 + k)
    return r1, k * r1


def _eclipse_flux(params, model, t0, rp, a_r, inc) -> np.ndarray:
    params.t0, params.rp, params.a, params.inc = t0, rp, a_r, inc
    return model.light_curve(params)


def _sector_model(theta, k_sector, models, period, t0, third_light):
    r_sum, k, sb, dt_sec = theta[:N_SHARED]
    base = N_SHARED + N_PER_SECTOR * k_sector
    dt0, inc, c0 = theta[base : base + N_PER_SECTOR]
    r1, r2 = _radii(r_sum, k)
    light = sb * k**2
    (pp, pm), sec = models
    f1 = _eclipse_flux(pp, pm, t0 + dt0, k, 1 / r1, inc)
    prim = c0 * ((1 - third_light) * (f1 + light) / (1 + light) + third_light)
    if sec is None:
        return prim, None
    sp, sm = sec
    f2 = _eclipse_flux(sp, sm, t0 + 0.5 * period + dt0 + dt_sec, 1 / k, 1 / r2, inc)
    second = c0 * ((1 - third_light) * (1 + light * f2) / (1 + light) + third_light)
    return prim, second


def _residuals(theta, windows, models, period, t0, third_light) -> np.ndarray:
    chunks = []
    for k_sector, (w, m) in enumerate(zip(windows, models)):
        prim, second = _sector_model(theta, k_sector, m, period, t0, third_light)
        chunks.append((w.primary.flux - prim) / w.primary.flux_err)
        if second is not None:
            chunks.append((w.secondary.flux - second) / w.secondary.flux_err)
    return np.concatenate(chunks)


def _make_eb_models(windows, period, t0, ld, batman) -> list:
    models = []
    for w in windows:
        p = transit_model(batman, w.primary.time, w.primary.cadence, period, t0, ld)
        s = None
        if w.secondary is not None:
            s = transit_model(batman, w.secondary.time, w.secondary.cadence, period, t0, ld)
        models.append((p, s))
    return models


def _initial_guess(windows, period, t0, duration, b0) -> tuple:
    def depth(s, centre):
        intr = np.abs(_phase(s.time, period, centre)) < 0.25 * duration
        base = np.median(s.flux[~intr]) if np.any(~intr) else 1.0
        return max(1 - np.median(s.flux[intr]) / base, 1e-5) if np.any(intr) else 1e-5

    r_sum = float(np.clip(np.pi * duration / period, R_SUM_BOUNDS[0] * 2, 0.6))
    d1 = np.median([depth(w.primary, t0) for w in windows])
    secs = [depth(w.secondary, t0 + 0.5 * period) for w in windows if w.secondary is not None]
    sb = float(np.clip(np.median(secs) / d1, 0.01, 3.0)) if secs else 0.5
    r1 = r_sum / 2
    inc0 = float(np.degrees(np.arccos(np.clip(b0 * r1, 0, 1))))
    dt_bound = max(0.5 * duration, MIN_DT0_BOUND)
    x0 = [r_sum, 1.0, sb, 0.0]
    lo = [R_SUM_BOUNDS[0], K_BOUNDS[0], SB_BOUNDS[0], -dt_bound]
    hi = [R_SUM_BOUNDS[1], K_BOUNDS[1], SB_BOUNDS[1], dt_bound]
    for w in windows:
        out = np.abs(_phase(w.primary.time, period, t0)) > 0.5 * duration
        c0 = float(np.median(w.primary.flux[out])) if np.any(out) else 1.0
        x0 += [0.0, min(inc0, 90.0 - 1e-3), c0]
        lo += [-dt_bound, INC_MIN_DEG, 0.5]
        hi += [dt_bound, 90.0, 1.5]
    return np.array(x0), np.array(lo), np.array(hi)


def _best_fit(windows, models, period, t0, duration, third_light):
    best = None
    for b0 in EB_B_STARTS:
        x0, lo, hi = _initial_guess(windows, period, t0, duration, b0)
        x0 = np.clip(x0, lo + 1e-9, hi - 1e-9)
        res = least_squares(
            _residuals,
            x0,
            bounds=(lo, hi),
            args=(windows, models, period, t0, third_light),
            x_scale="jac",
        )
        if best is None or res.cost < best.cost:
            best = res
    return best


def _impact(inc_deg, r1):
    return np.abs(np.cos(np.radians(inc_deg))) / r1


def _model_depth(batman, ld, period, t0, theta, third_light) -> float:
    r_sum, k, sb = theta[:3]
    r1, _ = _radii(r_sum, k)
    light = sb * k**2
    params, model = _single_point_model(batman, ld, period, t0)
    f1 = _eclipse_flux(params, model, t0, k, 1 / r1, theta[1 + N_SHARED])
    return float((1 - third_light) * (1 - f1[0]) / (1 + light))


def _single_point_model(batman, ld, period, t0):
    return transit_model(batman, np.array([t0]), 0.0, period, t0, ld)


def _sector_fit(k_sector, w, x, err, draws, period, t0, duration, batman, ld, third_light):
    base = N_SHARED + N_PER_SECTOR * k_sector
    i_dt0, i_inc = base, base + 1
    r1, _ = _radii(x[0], x[1])
    d_r1, _ = _radii(draws[:, 0], draws[:, 1])
    d_b = _impact(draws[:, i_inc], d_r1)
    d_t14 = transit_duration(period, 1 / d_r1, draws[:, 1], d_b)
    theta_k = np.concatenate([x[:N_SHARED], x[base : base + N_PER_SECTOR]])
    depth = _model_depth(batman, ld, period, t0, theta_k, third_light)
    depth_draws = [
        _model_depth(
            batman,
            ld,
            period,
            t0,
            np.concatenate([d[:N_SHARED], d[base : base + N_PER_SECTOR]]),
            third_light,
        )
        for d in draws[:200]
    ]
    b = float(_impact(x[i_inc], r1))
    data_depth, data_depth_err = empirical_depth(w.primary, period, t0, duration, x[i_dt0])
    p = w.primary
    return EclipseSectorFit(
        sector=p.span.sector,
        t_mid=p.t_mid,
        n_eclipses=p.n_transits,
        cadence=p.cadence,
        crowdsap=p.span.crowdsap,
        flfrcsap=p.span.flfrcsap,
        inc=float(x[i_inc]),
        inc_err=float(err[i_inc]),
        b=b,
        b_err=spread(d_b),
        depth=depth,
        depth_err=spread(np.array(depth_draws)),
        t14=float(transit_duration(period, 1 / r1, x[1], b)),
        t14_err=spread(d_t14),
        dt0=float(x[i_dt0]),
        dt0_err=float(err[i_dt0]),
        depth_data=data_depth,
        depth_data_err=data_depth_err,
        has_secondary=w.secondary is not None,
    )


def fit_sector_eclipses(
    windows: Sequence[EclipseWindows],
    period: float,
    t0: float,
    duration: float,
    limb_darkening=DEFAULT_LIMB_DARKENING,
    third_light: float = 0.0,
    seed: int = 0,
) -> tuple:
    """Joint eclipse fit; returns ``(fits, binary, inc_covariance)``.

    ``third_light`` is a fixed fraction of the total flux (0 when the light
    curve is already corrected for contamination and there is no known third
    star). Both stars share ``limb_darkening``.
    """
    if not 0 <= third_light < 1:
        raise ValueError(f"third_light must be in [0, 1), got {third_light}")
    batman = import_batman()
    models = _make_eb_models(windows, period, t0, limb_darkening, batman)
    res = _best_fit(windows, models, period, t0, duration, third_light)
    cov = covariance(res)
    err = np.sqrt(np.clip(np.diag(cov), 0, None))
    draws = np.random.default_rng(seed).multivariate_normal(
        res.x, cov, N_DRAWS, check_valid="ignore"
    )
    draws[:, 0] = np.clip(draws[:, 0], *R_SUM_BOUNDS)
    draws[:, 1] = np.clip(draws[:, 1], *K_BOUNDS)
    fits = tuple(
        _sector_fit(
            k, w, res.x, err, draws, period, t0, duration, batman, limb_darkening, third_light
        )
        for k, w in enumerate(windows)
    )
    binary = BinaryParams(
        r_sum=float(res.x[0]),
        r_sum_err=float(err[0]),
        k=float(res.x[1]),
        k_err=float(err[1]),
        sb_ratio=float(res.x[2]),
        sb_ratio_err=float(err[2]),
        dt_sec=float(res.x[3]),
        dt_sec_err=float(err[3]),
        third_light=float(third_light),
    )
    inc_idx = N_SHARED + 1 + N_PER_SECTOR * np.arange(len(windows))
    return fits, binary, cov[np.ix_(inc_idx, inc_idx)]


# --- verdict ----------------------------------------------------------------------


def _depth_slope(fits, t_ref) -> tuple:
    usable = [f for f in fits if np.isfinite(f.depth_data_err) and f.depth_data_err > 0]
    if len(usable) < MIN_SECTORS:
        return np.nan, np.nan
    t_yr = (np.array([f.t_mid for f in usable]) - t_ref) / DAYS_PER_YEAR
    trend = fit_trend(t_yr, [f.depth_data for f in usable], [f.depth_data_err for f in usable], 1)
    return trend.coeffs[1], trend.errors[1]


def _dilution_mimics_depth(fits) -> bool:
    crowd = np.array([f.crowdsap for f in fits], dtype=float)
    depth = np.array([f.depth_data for f in fits], dtype=float)
    ok = np.isfinite(crowd) & np.isfinite(depth) & (crowd > 0)
    if np.count_nonzero(ok) < MIN_SECTORS or np.median(depth[ok]) <= 0:
        return False
    rel_crowd = np.ptp(crowd[ok]) / np.median(crowd[ok])
    rel_depth = np.ptp(depth[ok]) / np.median(depth[ok])
    return bool(rel_crowd >= DILUTION_DEPTH_RATIO * rel_depth)


def eb_flags(fits, binary: Optional[BinaryParams]) -> tuple:
    """Flags that do not depend on the trend fit."""
    flags = list(data_flags(fits))
    if any(f.inc > 90.0 - INC_BOUNDARY_DEG for f in fits):
        flags.append("i_at_boundary")
    if "contamination_varies" in flags and _dilution_mimics_depth(fits):
        flags.append("dilution_mimics_depth")
    if not any(f.has_secondary for f in fits):
        flags.append("no_secondary")
    if (
        binary is not None
        and binary.dt_sec_err > 0
        and abs(binary.dt_sec) > T0_OFFSET_SIGMA * binary.dt_sec_err
    ):
        flags.append("eccentric")
    return tuple(flags)


def assess_eb_precession(
    fits: Sequence[EclipseSectorFit],
    binary: Optional[BinaryParams] = None,
    inc_covariance: Optional[np.ndarray] = None,
    min_baseline_yr: float = MIN_BASELINE_YR,
) -> EBPrecessionResult:
    """Trend fits, flags and the candidate decision for per-sector inclinations.

    A candidate needs ``|di/dt| >= 3 sigma``, ``BIC(static) - BIC(linear) >= 2``
    and none of: a baseline under ``min_baseline_yr``, ``i`` pinned at 90 deg,
    a measured depth trend opposite to ``di/dt`` at >= 2 sigma, or a CROWDSAP
    spread large enough to explain the depth change.
    """
    order = np.argsort([f.t_mid for f in fits])
    fits = tuple(fits[i] for i in order)
    if len(fits) < MIN_SECTORS:
        return EBPrecessionResult(
            fits, None, binary, reason=f"needs >= {MIN_SECTORS} sectors with eclipses"
        )
    t_mid = np.array([f.t_mid for f in fits])
    inc = np.array([f.inc for f in fits])
    inc_err = np.array([f.inc_err for f in fits])
    if not np.all(np.isfinite(inc_err) & (inc_err > 0)):
        return EBPrecessionResult(fits, None, binary, reason="unusable inclination uncertainties")
    uncertainty = inc_err if inc_covariance is None else inc_covariance[np.ix_(order, order)]
    try:
        models = fit_precession_model(t_mid, inc, uncertainty)
    except np.linalg.LinAlgError:
        return EBPrecessionResult(fits, None, binary, reason="unusable inclination covariance")
    di_dt, di_dt_err = models.linear.coeffs[1], models.linear.errors[1]
    significance = abs(di_dt) / di_dt_err if di_dt_err > 0 else np.nan
    ddepth_dt, ddepth_dt_err = _depth_slope(fits, models.t_ref)
    baseline_yr = float(np.ptp(t_mid) / DAYS_PER_YEAR)
    mean_inc = float(np.mean(inc))

    flags = list(eb_flags(fits, binary))
    if baseline_yr < min_baseline_yr:
        flags.append("short_baseline")
    if (
        np.isfinite(ddepth_dt_err)
        and ddepth_dt_err > 0
        and np.sign(ddepth_dt) == -np.sign(di_dt)
        and abs(ddepth_dt) / ddepth_dt_err >= DEPTH_MISMATCH_SIGMA
    ):
        flags.append("depth_sign_mismatch")
    if models.delta_bic_quadratic >= DELTA_BIC_MIN:
        flags.append("nonlinear")

    blocking = [f for f in flags if f in EB_BLOCKING_FLAGS]
    detected = significance >= SIGNIFICANCE_MIN and models.delta_bic >= DELTA_BIC_MIN
    if not detected:
        reason = f"no significant trend ({significance:.1f} sigma, dBIC={models.delta_bic:.1f})"
    elif blocking:
        reason = f"trend at {significance:.1f} sigma but flagged: {', '.join(blocking)}"
    else:
        reason = f"i drifts at {significance:.1f} sigma (dBIC={models.delta_bic:.1f})"

    r1 = _radii(binary.r_sum, binary.k)[0] if binary is not None else np.nan
    db_scale = np.sin(np.radians(mean_inc)) / r1 * np.radians(1.0)  # db per degree of i
    return EBPrecessionResult(
        sector_fits=fits,
        models=models,
        binary=binary,
        di_dt=di_dt,
        di_dt_err=di_dt_err,
        db_dt=float(-db_scale * di_dt),
        db_dt_err=float(db_scale * di_dt_err),
        significance=significance,
        delta_bic=models.delta_bic,
        ddepth_dt=ddepth_dt,
        ddepth_dt_err=ddepth_dt_err,
        baseline_yr=baseline_yr,
        max_nodal_period_yr=max_nodal_period(di_dt, mean_inc) if detected else np.nan,
        is_candidate=bool(detected and not blocking),
        flags=tuple(flags),
        reason=reason,
    )


def search_eb_precession(
    time,
    flux,
    flux_err,
    spans: Sequence[SectorSpan],
    period: float,
    t0: float,
    duration: float,
    limb_darkening=DEFAULT_LIMB_DARKENING,
    third_light: float = 0.0,
    min_baseline_yr: float = MIN_BASELINE_YR,
    seed: int = 0,
) -> EBPrecessionResult:
    """Top-level EB search on a stitched, flattened light curve with sector spans.

    ``period`` must be the binary period with ``t0`` at the deeper (primary)
    eclipse; a TLS period of ``P/2`` folds both eclipses together.
    """
    windows = extract_sector_eclipses(time, flux, flux_err, spans, period, t0, duration)
    if len(windows) < MIN_SECTORS:
        return EBPrecessionResult(
            (),
            None,
            reason=f"{len(windows)} sector(s) cover an eclipse; needs >= {MIN_SECTORS}",
        )
    fits, binary, inc_cov = fit_sector_eclipses(
        windows, period, t0, duration, limb_darkening, third_light, seed
    )
    return assess_eb_precession(fits, binary, inc_cov, min_baseline_yr)
