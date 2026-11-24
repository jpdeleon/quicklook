"""Eclipsing-binary precession branch (#63): triple-star geometry, eclipse fits, verdict.

Offline: synthetic detached-EB light curves (two limb-darkened spheres, built
here directly with batman) in four sectors one year apart.
"""

import numpy as np
import pytest

from quicklook.precession import DAYS_PER_YEAR, SectorSpan
from quicklook.precession_eb import (
    BinaryParams,
    EclipseSectorFit,
    assess_eb_precession,
    empirical_depth,
    extract_sector_eclipses,
    max_inclination_rate,
    max_nodal_period,
    nodal_precession_period,
    search_eb_precession,
)

batman = pytest.importorskip("batman")

PERIOD, T0 = 4.0, 1.0
R1, R2 = 0.12, 0.10  # radii in units of a
SB = 0.6  # surface-brightness ratio star 2 / star 1
INC_REF = 82.0  # b = cos(i)/R1 ~ 1.16: partial, as precessing EBs usually are
CADENCE = 10 / 1440
SECTOR_DAYS = 27.0
NOISE = 5e-4
SECTOR_STARTS = np.array([0.0, 1.0, 2.0, 3.0]) * DAYS_PER_YEAR + 2.0
LD = [0.35, 0.25]


def _duration():
    return PERIOD / np.pi * np.arcsin(np.sqrt((R1 + R2) ** 2 - np.cos(np.radians(INC_REF)) ** 2))


def _occultation(t, t0, rp, a_r, inc):
    params = batman.TransitParams()
    params.t0, params.per, params.rp, params.a, params.inc = t0, PERIOD, rp, a_r, inc
    params.ecc, params.w, params.limb_dark, params.u = 0.0, 90.0, "quadratic", LD
    return batman.TransitModel(params, t).light_curve(params)


def _eb_flux(t, inc, third_light=0.0, sb=SB):
    k = R2 / R1
    light = sb * k**2
    f1 = _occultation(t, T0, k, 1 / R1, inc)
    f2 = _occultation(t, T0 + PERIOD / 2, 1 / k, 1 / R2, inc)
    total = (f1 + light * f2) / (1 + light)
    return (1 - third_light) * total + third_light


def _stitched(di_dt, seed=0, starts=SECTOR_STARTS, crowdsap=None, sb=SB, inc_ref=INC_REF):
    t_ref = np.mean(starts + SECTOR_DAYS / 2)
    rng = np.random.default_rng(seed)
    times, fluxes, spans = [], [], []
    for k, start in enumerate(starts):
        inc = inc_ref + di_dt * (start + SECTOR_DAYS / 2 - t_ref) / DAYS_PER_YEAR
        t = np.arange(start, start + SECTOR_DAYS, CADENCE)
        times.append(t)
        fluxes.append(_eb_flux(t, inc, sb=sb) + rng.normal(0, NOISE, t.size))
        crowd = np.nan if crowdsap is None else crowdsap[k]
        spans.append(SectorSpan(sector=20 + k, tmin=t[0], tmax=t[-1], crowdsap=crowd))
    flux = np.concatenate(fluxes)
    return np.concatenate(times), flux, np.full_like(flux, NOISE), spans


def _search(di_dt, seed=0, **kwargs):
    time, flux, err, spans = _stitched(di_dt, seed=seed, **kwargs)
    return search_eb_precession(time, flux, err, spans, PERIOD, T0, _duration())


def _fit(sector, t_mid, inc, inc_err=0.02, depth_data=0.05, depth_data_err=5e-4, **kwargs):
    values = dict(
        sector=sector,
        t_mid=t_mid,
        n_eclipses=6,
        cadence=CADENCE,
        crowdsap=np.nan,
        flfrcsap=np.nan,
        inc=inc,
        inc_err=inc_err,
        b=1.0,
        b_err=0.01,
        depth=depth_data,
        depth_err=depth_data_err,
        t14=0.2,
        t14_err=0.001,
        dt0=0.0,
        dt0_err=1e-3,
        depth_data=depth_data,
        depth_data_err=depth_data_err,
        has_secondary=True,
    )
    values.update(kwargs)
    return EclipseSectorFit(**values)


BINARY = BinaryParams(
    r_sum=R1 + R2,
    r_sum_err=0.001,
    k=R2 / R1,
    k_err=0.01,
    sb_ratio=SB,
    sb_ratio_err=0.01,
    dt_sec=0.0,
    dt_sec_err=1e-3,
    third_light=0.0,
)


def _trend_fits(rate=1.0, depth_rate=0.004, years=(0, 1, 2, 3), **kwargs):
    return [
        _fit(
            20 + n,
            y * DAYS_PER_YEAR,
            INC_REF + rate * (y - 1.5),
            depth_data=0.05 + depth_rate * (y - 1.5),
            **kwargs,
        )
        for n, y in enumerate(years)
    ]


# --- triple-star geometry ---------------------------------------------------------


def test_nodal_period_reproduces_v907_sco():
    # Zasche et al. 2023: P1 = 3.7762 d, P2 = 142.01 d, e2 = 0.467, masses 2.56 +
    # 2.74 Msun, m3 = 1.06 Msun, J = 26.2 deg -> P_nu = 63.5 (+3.3/-2.6) yr.
    period_days = nodal_precession_period(3.7762049, 142.01, 2.56, 2.74, 1.06, 26.2, 0.467)
    assert period_days / DAYS_PER_YEAR == pytest.approx(63.5, abs=1.0)


def test_nodal_period_grows_with_outer_period_and_lighter_third_star():
    base = nodal_precession_period(3.0, 100.0, 1.0, 1.0, 1.0, 20.0)
    assert nodal_precession_period(3.0, 200.0, 1.0, 1.0, 1.0, 20.0) > base
    assert nodal_precession_period(3.0, 100.0, 1.0, 1.0, 0.5, 20.0) > base


def test_nodal_period_rejects_unphysical_input():
    with pytest.raises(ValueError):
        nodal_precession_period(3.0, 100.0, 1.0, 1.0, 1.0, 20.0, e_outer=1.0)
    with pytest.raises(ValueError):
        nodal_precession_period(3.0, 100.0, 1.0, 1.0, 0.0, 20.0)


def test_nodal_period_is_infinite_at_ninety_degrees():
    assert np.isinf(nodal_precession_period(3.0, 100.0, 1.0, 1.0, 1.0, 90.0))


def test_v907_sco_rate_is_within_its_precession_bound():
    # Table 2: i = 73.37 (S13, HJD 58675) -> 76.59 (S39, HJD ~59375): ~1.7 deg/yr.
    rate = (76.59 - 73.37) / ((59375 - 58675) / DAYS_PER_YEAR)
    bound = max_inclination_rate(63.5, 26.2, 75.0)
    assert rate < bound
    assert max_nodal_period(rate, 75.0) > 63.5


def test_max_nodal_period_needs_a_rate():
    assert np.isnan(max_nodal_period(0.0, 80.0))
    assert np.isnan(max_inclination_rate(0.0, 20.0, 80.0))


# --- data ----------------------------------------------------------------------------


def test_extract_sector_eclipses_pairs_primary_and_secondary():
    time, flux, err, spans = _stitched(0.0)
    windows = extract_sector_eclipses(time, flux, err, spans, PERIOD, T0, _duration())
    assert [w.primary.span.sector for w in windows] == [20, 21, 22, 23]
    assert all(w.secondary is not None for w in windows)
    phase = (windows[0].secondary.time - T0) % PERIOD
    assert np.all(np.abs(phase - PERIOD / 2) < PERIOD / 4)


def test_empirical_depth_measures_primary_eclipse():
    time, flux, err, spans = _stitched(0.0)
    windows = extract_sector_eclipses(time, flux, err, spans, PERIOD, T0, _duration())
    depth, depth_err = empirical_depth(windows[0].primary, PERIOD, T0, _duration())
    t_mid = np.array([T0])
    truth = 1 - _eb_flux(t_mid, INC_REF)[0]
    assert depth == pytest.approx(truth, rel=0.1)
    assert 0 < depth_err < 0.1 * depth


# --- verdict --------------------------------------------------------------------------


def test_assess_candidate_for_clean_inclination_trend():
    result = assess_eb_precession(_trend_fits(), BINARY)
    assert result.is_candidate, result.reason
    assert result.di_dt == pytest.approx(1.0, rel=1e-6)
    # b = cos(i)/R1 shrinks as i grows toward 90 deg.
    assert result.db_dt < 0
    assert np.isfinite(result.max_nodal_period_yr)


def test_assess_blocks_depth_trend_opposite_to_inclination():
    result = assess_eb_precession(_trend_fits(depth_rate=-0.004), BINARY)
    assert "depth_sign_mismatch" in result.flags
    assert not result.is_candidate


def test_assess_blocks_inclination_at_ninety_degrees():
    fits = _trend_fits()
    fits[-1] = _fit(23, 3 * DAYS_PER_YEAR, 89.95, depth_data=0.056)
    result = assess_eb_precession(fits, BINARY)
    assert "i_at_boundary" in result.flags
    assert not result.is_candidate


def test_assess_blocks_when_crowding_spread_explains_depth_change():
    crowds = (0.80, 0.86, 0.92, 0.98)
    fits = [
        _fit(f.sector, f.t_mid, f.inc, depth_data=f.depth_data, crowdsap=c)
        for f, c in zip(_trend_fits(depth_rate=0.001), crowds)
    ]
    result = assess_eb_precession(fits, BINARY)
    assert {"contamination_varies", "dilution_mimics_depth"} <= set(result.flags)
    assert not result.is_candidate


def test_assess_blocks_short_baseline_and_flags_missing_secondary():
    fits = _trend_fits(rate=5.0, years=(0, 0.2, 0.4), has_secondary=False)
    result = assess_eb_precession(fits, BINARY)
    assert {"short_baseline", "no_secondary"} <= set(result.flags)
    assert not result.is_candidate


def test_assess_flags_eccentric_secondary_offset():
    binary = BinaryParams(**{**BINARY.__dict__, "dt_sec": 0.01, "dt_sec_err": 1e-3})
    assert "eccentric" in assess_eb_precession(_trend_fits(), binary).flags


def test_assess_needs_two_sectors():
    result = assess_eb_precession(_trend_fits(years=(0,)), BINARY)
    assert result.models is None and not result.is_candidate


# --- end to end -------------------------------------------------------------------------


def test_search_recovers_injected_inclination_drift():
    result = _search(1.0)
    assert result.is_candidate, result.reason
    assert result.di_dt == pytest.approx(1.0, abs=3 * result.di_dt_err + 0.05)
    assert result.binary.k == pytest.approx(R2 / R1, rel=0.15)
    assert result.binary.sb_ratio == pytest.approx(SB, rel=0.2)
    assert result.ddepth_dt > 0  # deeper as i -> 90 deg


def test_search_recovers_decreasing_inclination():
    result = _search(-0.8, seed=3)
    assert result.is_candidate, result.reason
    assert result.di_dt == pytest.approx(-0.8, abs=3 * result.di_dt_err + 0.05)
    assert result.db_dt > 0


def test_search_static_binary_is_not_a_candidate():
    result = _search(0.0, seed=1)
    assert not result.is_candidate
    assert result.significance < 3


def test_search_null_false_positive_rate_below_five_percent():
    hits = sum(_search(0.0, seed=100 + s).is_candidate for s in range(10))
    assert hits / 10 <= 0.05


def test_search_without_secondary_light_flags_nothing_spurious():
    result = _search(1.0, sb=0.0)
    assert result.is_candidate, result.reason
    assert result.binary.sb_ratio < 0.05


def test_search_with_one_sector_returns_reason():
    time, flux, err, spans = _stitched(0.0, starts=SECTOR_STARTS[:1])
    result = search_eb_precession(time, flux, err, spans, PERIOD, T0, _duration())
    assert result.models is None
    assert "needs >= 2" in result.reason


def test_search_rejects_bad_third_light():
    time, flux, err, spans = _stitched(0.0)
    with pytest.raises(ValueError):
        search_eb_precession(time, flux, err, spans, PERIOD, T0, _duration(), third_light=1.0)


def test_unconverged_eclipse_fits_do_not_produce_verdict(monkeypatch):
    from functools import partial

    from scipy.optimize import least_squares

    monkeypatch.setattr("quicklook.precession_eb.least_squares", partial(least_squares, max_nfev=1))
    result = _search(1.0)
    assert not result.is_candidate
    assert result.models is None and result.binary is None
    assert "did not converge" in result.reason


def test_primary_only_fit_leaves_secondary_timing_unavailable():
    time, flux, err, spans = _stitched(1.0)
    primary_half = np.abs((time - T0 + PERIOD / 2) % PERIOD - PERIOD / 2) < PERIOD / 4
    result = search_eb_precession(
        time[primary_half], flux[primary_half], err[primary_half], spans, PERIOD, T0, _duration()
    )
    assert "no_secondary" in result.flags
    assert np.isnan(result.binary.dt_sec)
    assert np.isnan(result.binary.dt_sec_err)


def test_opposite_crowding_trend_keeps_candidate():
    from dataclasses import replace

    fits = [replace(f, crowdsap=c) for f, c in zip(_trend_fits(), (0.98, 0.94, 0.89, 0.84))]
    result = assess_eb_precession(fits, BINARY)
    assert "contamination_varies" in result.flags
    assert result.is_candidate
