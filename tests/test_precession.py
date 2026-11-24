"""Precession search (#63): sector extraction, per-sector fits, trend tests, flags.

Offline: synthetic batman light curves in four sectors spread over a few years.
"""

import numpy as np
import pytest


from quicklook.precession import (
    DAYS_PER_YEAR,
    SectorSpan,
    SectorTransitFit,
    a_rs_from_star,
    assess_precession,
    correct_dilution,
    extract_sector_transits,
    fit_precession_model,
    inclination_rate,
    search_precession,
    systematics_flags,
    transit_duration,
    years_until_no_transit,
)

batman = pytest.importorskip("batman")

PERIOD, T0, RP, A_RS = 3.0, 1.0, 0.1, 10.0
CADENCE = 10 / 1440
SECTOR_DAYS = 27.0
NOISE = 4e-4
SECTOR_STARTS = np.array([0.0, 1.0, 2.0, 3.0]) * DAYS_PER_YEAR + 2.0
B_REF = 0.8  # recent detections sit at high b, where b is well constrained


def _duration(b=B_REF):
    return float(transit_duration(PERIOD, A_RS, RP, b))


def _sector_lc(start, b, seed, rp=RP, cadence=CADENCE):
    t = np.arange(start, start + SECTOR_DAYS, cadence)
    params = batman.TransitParams()
    params.t0, params.per, params.rp, params.a = T0, PERIOD, rp, A_RS
    params.inc = float(np.degrees(np.arccos(b / A_RS)))
    params.ecc, params.w, params.limb_dark, params.u = 0.0, 90.0, "quadratic", [0.35, 0.25]
    flux = batman.TransitModel(params, t).light_curve(params)
    rng = np.random.default_rng(seed)
    return t, flux + rng.normal(0, NOISE, t.size)


def _stitched(db_dt, seed=0, starts=SECTOR_STARTS, rp=RP, crowdsap=None, b_ref=B_REF):
    t_ref = np.mean(starts + SECTOR_DAYS / 2)
    times, fluxes, spans = [], [], []
    for k, start in enumerate(starts):
        b = b_ref + db_dt * (start + SECTOR_DAYS / 2 - t_ref) / DAYS_PER_YEAR
        t, f = _sector_lc(start, b, seed * 100 + k, rp=rp)
        times.append(t)
        fluxes.append(f)
        crowd = np.nan if crowdsap is None else crowdsap[k]
        spans.append(SectorSpan(sector=10 + k, tmin=t[0], tmax=t[-1], crowdsap=crowd))
    time, flux = np.concatenate(times), np.concatenate(fluxes)
    return time, flux, np.full_like(flux, NOISE), spans


def _search(db_dt, seed=0, a_rs_prior=None, **kwargs):
    time, flux, err, spans = _stitched(db_dt, seed=seed, **kwargs)
    return search_precession(time, flux, err, spans, PERIOD, T0, _duration(), a_rs_prior=a_rs_prior)


def _fit(sector, t_mid, b, b_err=0.01, t14=0.12, t14_err=0.001, **kwargs):
    values = dict(
        sector=sector,
        t_mid=t_mid,
        n_transits=8,
        cadence=CADENCE,
        crowdsap=np.nan,
        flfrcsap=np.nan,
        b=b,
        b_err=b_err,
        rp_rs=RP,
        rp_rs_err=0.002,
        depth=RP**2,
        depth_err=4e-4,
        model_depth=RP**2,
        t14=t14,
        t14_err=t14_err,
        dt0=0.0,
        dt0_err=1e-3,
        corr_b_rp=0.3,
    )
    values.update(kwargs)
    return SectorTransitFit(**values)


# --- data helpers ---------------------------------------------------------------


def test_extract_sector_transits_splits_by_span_and_counts_transits():
    time, flux, err, spans = _stitched(0.0)
    sectors = extract_sector_transits(time, flux, err, spans, PERIOD, T0, _duration())
    assert [s.span.sector for s in sectors] == [10, 11, 12, 13]
    for s, span in zip(sectors, spans):
        assert span.tmin <= s.time.min() and s.time.max() <= span.tmax
        assert s.n_transits in (8, 9, 10)
        assert s.cadence == pytest.approx(CADENCE, rel=1e-3)


def test_extract_sector_transits_drops_sector_without_transit():
    time, flux, err, spans = _stitched(0.0)
    empty = SectorSpan(sector=99, tmin=time.max() + 10, tmax=time.max() + 20)
    sectors = extract_sector_transits(time, flux, err, [*spans, empty], PERIOD, T0, _duration())
    assert 99 not in [s.span.sector for s in sectors]


def test_extract_sector_transits_fills_missing_errors_with_scatter():
    time, flux, _, spans = _stitched(0.0)
    sectors = extract_sector_transits(time, flux, None, spans, PERIOD, T0, _duration())
    assert all(np.allclose(s.flux_err, s.flux_err[0]) for s in sectors)
    assert sectors[0].flux_err[0] == pytest.approx(NOISE, rel=0.2)


def test_correct_dilution_deepens_transit_by_crowding():
    flux = np.array([1.0, 1.0, 0.99, 1.0])
    corrected = correct_dilution(flux, crowdsap=0.5)
    assert 1 - corrected[2] == pytest.approx(0.02)


def test_correct_dilution_without_keyword_only_normalises():
    flux = np.array([2.0, 2.0, 1.98, 2.0])
    np.testing.assert_allclose(correct_dilution(flux, np.nan), flux / 2.0)


def test_transit_duration_shrinks_with_impact_parameter():
    assert _duration(0.0) > _duration(0.5) > _duration(0.9) > 0
    assert transit_duration(PERIOD, A_RS, RP, 1.2) == 0


# --- trend models ---------------------------------------------------------------


def test_fit_precession_model_recovers_linear_slope():
    t = np.array([0.0, 1.0, 2.0, 3.0, 4.0]) * DAYS_PER_YEAR
    b = 0.6 + 0.02 * (t - t.mean()) / DAYS_PER_YEAR
    models = fit_precession_model(t, b, np.full(5, 0.002))
    assert models.linear.coeffs[1] == pytest.approx(0.02, abs=1e-9)
    assert models.delta_bic > 10
    assert models.quadratic is not None


def test_fit_precession_model_prefers_static_for_constant_b():
    t = np.array([0.0, 1.0, 2.0]) * DAYS_PER_YEAR
    models = fit_precession_model(t, [0.5, 0.501, 0.499], [0.01, 0.01, 0.01])
    assert models.delta_bic < 0
    assert models.quadratic is None
    assert np.isnan(models.delta_bic_quadratic)


# --- derived quantities -----------------------------------------------------------


def test_inclination_rate_sign_and_size():
    rate = inclination_rate(0.01, b=0.5, a_rs=10.0)
    expected = -np.degrees(0.01 / (10 * np.sqrt(1 - 0.05**2)))
    assert rate == pytest.approx(expected)


def test_a_rs_from_star_earth_sun():
    a_rs, sigma = a_rs_from_star(DAYS_PER_YEAR, 1.0, 1.0, 0.01, 0.03)
    assert a_rs == pytest.approx(215.0, rel=1e-3)
    assert sigma / a_rs == pytest.approx(np.hypot(0.01, 0.01))


def test_a_rs_from_star_floor_and_invalid():
    assert a_rs_from_star(3.0, 1.0, 1.0)[1] / a_rs_from_star(3.0, 1.0, 1.0)[0] == 0.1
    assert a_rs_from_star(3.0, np.nan, 1.0) is None


def test_years_until_no_transit_matches_toi1355_style_numbers():
    # b=0.9275 in 2024 at 0.0163/yr with k~0.07 leaves the disc around 2033
    assert years_until_no_transit(0.9275, 0.07, 0.0163) == pytest.approx(8.74, abs=0.05)
    assert np.isnan(years_until_no_transit(0.5, 0.1, -0.01))


# --- flags and verdict ------------------------------------------------------------


def _trend_fits(db_dt, years=(0, 1, 2, 3), **kwargs):
    fits = []
    for k, yr in enumerate(years):
        b = 0.5 + db_dt * yr
        fits.append(_fit(k, yr * DAYS_PER_YEAR, b, t14=_duration(b), **kwargs))
    return fits


def test_assess_flags_candidate_for_clean_trend():
    result = assess_precession(_trend_fits(0.05), a_rs=A_RS)
    assert result.is_candidate
    assert result.db_dt == pytest.approx(0.05)
    assert result.baseline_yr == pytest.approx(3.0)
    assert result.dt14_dt < 0
    assert result.di_dt < 0


def test_assess_blocks_short_baseline():
    fits = _trend_fits(1.0, years=(0, 0.1, 0.2, 0.3))
    result = assess_precession(fits, a_rs=A_RS)
    assert "short_baseline" in result.flags
    assert not result.is_candidate


def test_assess_blocks_when_t14_moves_with_b():
    fits = [
        _fit(k, yr * DAYS_PER_YEAR, 0.5 + 0.05 * yr, t14=0.12 + 0.005 * yr)
        for k, yr in enumerate((0, 1, 2, 3))
    ]
    result = assess_precession(fits, a_rs=A_RS)
    assert "shape_sign_mismatch" in result.flags
    assert not result.is_candidate


def test_assess_grazing_correlation_alone_does_not_block():
    result = assess_precession(_trend_fits(0.05, corr_b_rp=0.99), a_rs=A_RS)
    assert "grazing_degeneracy" in result.flags
    assert "rp_absorbs_b" not in result.flags
    assert result.is_candidate


def test_assess_blocks_when_rp_rs_trends_with_degenerate_b():
    fits = [
        _fit(
            k,
            yr * DAYS_PER_YEAR,
            0.5 + 0.05 * yr,
            t14=_duration(0.5 + 0.05 * yr),
            corr_b_rp=0.99,
            rp_rs=RP - 0.01 * yr,
        )
        for k, yr in enumerate((0, 1, 2, 3))
    ]
    result = assess_precession(fits, a_rs=A_RS)
    assert {"grazing_degeneracy", "depth_change", "rp_absorbs_b"} <= set(result.flags)
    assert not result.is_candidate


def test_assess_blocks_b_at_zero_bound():
    fits = _trend_fits(0.05)
    fits[0] = _fit(0, 0.0, 0.01, t14=_duration(0.01))
    result = assess_precession(fits, a_rs=A_RS)
    assert "b_at_boundary" in result.flags
    assert not result.is_candidate


def test_assess_needs_two_sectors():
    result = assess_precession([_fit(1, 0.0, 0.5)])
    assert result.models is None and not result.is_candidate


def test_systematics_flags_contamination_cadence_and_t0():
    fits = [
        _fit(1, 0.0, 0.5, crowdsap=0.95),
        _fit(2, 400.0, 0.5, crowdsap=0.99, cadence=30 / 1440, dt0=0.01, dt0_err=1e-3),
    ]
    flags = systematics_flags(fits)
    assert {"contamination_varies", "mixed_cadence", "t0_offset"} <= set(flags)


# TOI-1355 b against its TLS ephemeris: a period error makes dt0 grow linearly.
TOI1355_T_MID = (2868.78, 2894.81, 3404.82, 3443.88, 3597.97)
TOI1355_DT0 = (-0.00680, -0.00602, 0.00727, 0.00826, 0.01211)


def test_t0_offset_ignores_linear_ephemeris_drift():
    fits = [
        _fit(k, t, 0.9, dt0=d, dt0_err=3e-4)
        for k, (t, d) in enumerate(zip(TOI1355_T_MID, TOI1355_DT0))
    ]
    assert "t0_offset" not in systematics_flags(fits)


def test_t0_offset_flags_a_sector_off_the_refitted_ephemeris():
    dt0 = list(TOI1355_DT0)
    dt0[2] += 0.005  # 7 min, ~17 sigma off the drift line
    fits = [
        _fit(k, t, 0.9, dt0=d, dt0_err=3e-4) for k, (t, d) in enumerate(zip(TOI1355_T_MID, dt0))
    ]
    assert "t0_offset" in systematics_flags(fits)


# --- end to end ---------------------------------------------------------------------


def test_search_recovers_injected_db_dt():
    rate = 0.04
    result = _search(rate)
    assert len(result.sector_fits) == 4
    assert abs(result.db_dt - rate) < 3 * result.db_dt_err
    assert result.a_rs == pytest.approx(A_RS, rel=0.15)
    assert result.is_candidate, result.reason


def test_search_moderate_b_recovers_with_a_rs_prior():
    rate = 0.04
    result = _search(rate, b_ref=0.5, a_rs_prior=(A_RS, 0.2))
    assert abs(result.db_dt - rate) < 3 * result.db_dt_err
    assert result.a_rs == pytest.approx(A_RS, rel=0.05)


def test_search_static_transit_is_not_a_candidate():
    result = _search(0.0)
    assert not result.is_candidate
    assert result.significance < 3


def test_search_records_crowding_and_flags_its_change():
    result = _search(0.0, crowdsap=[0.90, 0.95, 0.99, 0.99])
    assert [f.crowdsap for f in result.sector_fits] == [0.90, 0.95, 0.99, 0.99]
    assert "contamination_varies" in result.flags


def test_search_with_one_sector_returns_reason():
    time, flux, err, spans = _stitched(0.0, starts=SECTOR_STARTS[:1])
    result = search_precession(time, flux, err, spans, PERIOD, T0, _duration())
    assert result.models is None and "needs" in result.reason


@pytest.mark.parametrize("rate", [0.03, -0.03])
def test_injection_recovery_grid(rate):
    recovered = [_search(rate, seed=s) for s in range(3)]
    assert all(r.is_candidate for r in recovered)
    assert all(np.sign(r.db_dt) == np.sign(rate) for r in recovered)


def test_null_false_positive_rate_below_five_percent():
    trials = [_search(0.0, seed=s) for s in range(20)]
    false_positives = sum(r.is_candidate for r in trials)
    assert false_positives / len(trials) <= 0.05


def test_extract_keeps_minimum_window_when_duration_is_underestimated():
    time, flux, err, spans = _stitched(0.0)
    short = 0.01  # a TLS duration far below the true T14
    sectors = extract_sector_transits(time, flux, err, spans, PERIOD, T0, short)
    phase = (sectors[0].time - T0 + PERIOD / 2) % PERIOD - PERIOD / 2
    assert np.abs(phase).max() > 0.2


def test_search_recovers_trend_with_underestimated_duration():
    time, flux, err, spans = _stitched(0.04)
    result = search_precession(time, flux, err, spans, PERIOD, T0, 0.3 * _duration())
    assert abs(result.db_dt - 0.04) < 3 * result.db_dt_err


@pytest.mark.parametrize(
    ("r_err", "m_err", "expected"), [(0.3, np.nan, 0.3), (np.nan, 0.9, 0.3), (0.02, np.nan, 0.1)]
)
def test_partial_stellar_errors_respect_known_uncertainty(r_err, m_err, expected):
    a_rs, sigma = a_rs_from_star(3, 1, 1, r_err, m_err)
    assert sigma / a_rs == pytest.approx(expected)


def test_search_skips_sectors_with_no_measurable_noise():
    time, flux, _, spans = _stitched(0)
    result = search_precession(time, np.ones_like(flux), None, spans, PERIOD, T0, _duration())
    assert not result.is_candidate
    assert not result.sector_fits


def test_shared_offset_uncertainty_cancels_in_two_sector_drift():
    # A common offset cancels in the one-year difference. Only the independent
    # 0.01 errors add in quadrature, even with a much larger shared error.
    fits = _trend_fits(0.05, years=(0, 1))
    covariance = np.full((2, 2), 0.1**2) + np.eye(2) * 0.01**2
    result = assess_precession(fits, b_covariance=covariance)
    assert result.db_dt_err == pytest.approx(np.sqrt(2) * 0.01)
    assert result.is_candidate


def test_joint_fit_retains_sensitivity_to_moderate_b_drift():
    data = _stitched(0.012, b_ref=0.55)
    result = search_precession(*data, PERIOD, T0, _duration(0.55))
    assert result.significance >= 3


def test_plot_sector_folds_has_one_panel_per_sector_plus_overlay():
    from quicklook.plot_precession import plot_sector_folds

    time, flux, err, spans = _stitched(0.04)
    result = search_precession(time, flux, err, spans, PERIOD, T0, _duration())
    fig = plot_sector_folds(result, time, flux, err, spans, PERIOD, T0, _duration(), "t")
    visible = [ax for ax in fig.axes if ax.get_visible()]
    assert len(visible) == len(result.sector_fits) + 1
    assert "b = " in visible[0].get_title()


def test_sector_fold_uses_fitted_baseline_without_out_of_transit_points():
    import matplotlib.pyplot as pl

    from quicklook.plot_precession import plot_sector_folds

    time, flux, err, spans = _stitched(0.04)
    first = time <= spans[0].tmax
    phase = (time - T0 + PERIOD / 2) % PERIOD - PERIOD / 2
    keep = ~first | (np.abs(phase) < 0.095)
    flux[first] *= 1.1
    err[first] *= 1.1
    data = (time[keep], flux[keep], err[keep], spans, PERIOD, T0, 0.22)
    result = search_precession(*data, a_rs_prior=(A_RS, 0.1))
    assert result.is_candidate
    baseline = result.sector_fits[0].baseline
    assert baseline == pytest.approx(1.1, abs=0.001)
    fig = plot_sector_folds(result, *data)
    raw, binned = fig.axes[0].lines[:2]
    np.testing.assert_allclose(raw.get_ydata(), flux[first & keep] / baseline)
    assert np.isfinite(binned.get_ydata()).all()
    pl.close(fig)
