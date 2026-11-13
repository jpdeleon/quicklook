"""Dedicated detrending for the single-event search. Offline, synthetic light curves."""

import numpy as np
import pytest

from quicklook.single_event_detrend import (
    MASK_FACTOR,
    choose_window,
    detrend,
    event_mask,
    is_fast_rotator,
    window_floor,
)
from quicklook.single_event import SingleEvent

CADENCE = 120 / 86400


def _event(t0=13.0, duration=22 / 24):
    return SingleEvent(t0, duration, 0.004, 40.0, 500.0, 400.0, 3.0, 600)


def _raw(variability_period=6.0, amp=2e-3, transit=True, seed=0):
    rng = np.random.default_rng(seed)
    t = np.arange(0.0, 27.0, CADENCE)
    f = 1 + amp * np.sin(2 * np.pi * t / variability_period) + rng.normal(0, 3e-4, t.size)
    if transit:
        f[np.abs(t - 13.0) < 11 / 24] -= 0.004
    return t, f


# --- window choice -------------------------------------------------------------


def test_window_floor_preserves_the_longest_trial_duration():
    from quicklook.single_event import duration_grid

    grid = duration_grid(span=27.0, cadence=CADENCE, rho_star=1.0)
    assert window_floor(1.0) == pytest.approx(3 * grid[-1], rel=0.01)
    assert 2.5 < window_floor(1.0) < 3.5  # ~2.9 d for a Sun-like star
    assert window_floor(0.24) > window_floor(1.0)  # lower density -> longer transits
    assert window_floor(0.01) == pytest.approx(4.5)  # capped via MAX_DURATION


def test_default_window_is_raised_to_the_floor():
    window, reason = choose_window(0.5, rho_star=1.0, user_set=False, fast_rotator=False)
    assert window == pytest.approx(window_floor(1.0))
    assert reason == "floor"


def test_wider_pipeline_window_is_kept():
    window, reason = choose_window(3.0, rho_star=1.0, user_set=False, fast_rotator=False)
    assert window == 3.0
    assert reason == "pipeline"


def test_user_window_is_respected():
    window, reason = choose_window(0.6, rho_star=1.0, user_set=True, fast_rotator=False)
    assert window == 0.6
    assert reason == "user"


def test_fast_rotator_keeps_the_pipeline_window():
    window, reason = choose_window(0.5, rho_star=1.0, user_set=False, fast_rotator=True)
    assert window == 0.5
    assert reason == "fast rotator"


# --- rotation check ---------------------------------------------------------------


def test_fast_rotator_detected():
    t, f = _raw(variability_period=1.0, amp=5e-3, transit=False)
    assert is_fast_rotator(t, f)


def test_slow_rotator_is_not_fast():
    t, f = _raw(variability_period=6.0, amp=5e-3, transit=False)
    assert not is_fast_rotator(t, f)


def test_quiet_star_is_not_fast():
    t, f = _raw(amp=0.0, transit=False)
    assert not is_fast_rotator(t, f)


# --- masks and detrending ------------------------------------------------------------


def test_event_mask_covers_the_event_plus_margin():
    t = np.arange(0.0, 27.0, CADENCE)
    ev = _event()
    mask = event_mask(t, [ev])
    half = ev.duration / 2 + MASK_FACTOR * ev.duration / 2
    assert mask[np.argmin(np.abs(t - 13.0))]
    assert not mask[np.argmin(np.abs(t - (13.0 + half + 0.05)))]
    assert mask.sum() == pytest.approx(2 * half / CADENCE, rel=0.02)


def test_event_mask_empty_without_events():
    assert not event_mask(np.arange(0.0, 5.0, CADENCE), []).any()


def _in_transit_depth(t, flat):
    return 1 - np.nanmedian(flat[np.abs(t - 13.0) < 8 / 24])


@pytest.mark.parametrize("window", [1.6, 3.0])
def test_masked_detrending_preserves_a_long_transit_on_slow_variability(window):
    t, f = _raw(variability_period=20.0)
    plain, _ = detrend(t, f, window=window)
    masked, trend = detrend(t, f, window=window, mask=event_mask(t, [_event()]))
    assert _in_transit_depth(t, masked) == pytest.approx(0.004, rel=0.1)
    if window < 3 * _event().duration:
        # Window under ~3x the duration: unmasked, the trend sags into the dip.
        assert _in_transit_depth(t, plain) < 0.9 * _in_transit_depth(t, masked)
    assert np.all(np.isfinite(trend[np.isfinite(masked)]))


def test_masking_helps_even_when_variability_is_faster():
    # P = 6 d: neither method is exact (see detrend's docstring), masking is closer.
    t, f = _raw(variability_period=6.0)
    plain, _ = detrend(t, f, window=1.6)
    masked, _ = detrend(t, f, window=1.6, mask=event_mask(t, [_event()]))
    assert abs(_in_transit_depth(t, masked) - 0.004) < abs(_in_transit_depth(t, plain) - 0.004)


def test_masked_detrending_removes_variability_outside_the_mask():
    t, f = _raw()
    masked, _ = detrend(t, f, window=1.6, mask=event_mask(t, [_event()]))
    out = np.abs(t - 13.0) > 1.5
    assert np.nanstd(masked[out]) < 5e-4  # sine (2 ppt) removed, white noise remains


def test_detrend_handles_non_finite_flux():
    t, f = _raw()
    f[::97] = np.nan
    flat, trend = detrend(t, f, window=1.6)
    assert flat.shape == trend.shape == t.shape
    assert np.all(np.isnan(flat[::97]))  # bad input stays bad
    good = np.ones(t.size, dtype=bool)
    good[::97] = False
    assert np.isfinite(flat[good]).mean() > 0.99


# --- pipeline wiring -------------------------------------------------------------------


def _pipeline_ql(window_length=0.5, user_set=False, variability_period=20.0):
    from lightkurve import LightCurve

    from quicklook.tql import TessQuickLook, _TLSResult

    t, f = _raw(variability_period=variability_period)
    raw = LightCurve(time=t, flux=f, flux_err=np.full_like(f, 3e-4))
    ql = TessQuickLook.__new__(TessQuickLook)
    ql.verbose = False
    ql.single_event_search = True
    ql.raw_lc = raw
    # TLS's flattening: a 0.5 d window that eats most of the 22 h transit.
    tls_flat, _ = detrend(t, f, window=0.5)
    ql.flat_lc = LightCurve(time=t, flux=tls_flat, flux_err=np.full_like(f, 3e-4))
    ql.window_length = window_length
    ql._use_rotation_aware_window = not user_set
    ql.exofop_data = {}
    ql.tls_results = _TLSResult(
        {"SDE": 4.0, "period": 3.3, "T0": 1.0, "duration": 0.1, "distinct_transit_count": 9}
    )
    return ql


def test_pipeline_uses_its_own_wider_window_and_masks_events():
    ql = _pipeline_ql()
    ql.run_single_event_search()
    r = ql.tls_results
    assert r["single_event_window"] == pytest.approx(window_floor(1.0))
    assert r["single_event_window_reason"] == "floor"
    assert r["single_event_redetrended"] is True
    assert r["single_event_max_duration"] == pytest.approx(window_floor(1.0) / 3)
    assert r["event_class"] == "mono"
    assert r["single_event_depth"] == pytest.approx(0.004, rel=0.15)
    # TLS's own flattened light curve is untouched.
    assert ql.se_flux is not None and ql.se_flux is not ql.flat_lc.flux.value


def test_pipeline_respects_a_user_window_and_caps_durations():
    from loguru import logger

    messages = []
    sink = logger.add(messages.append, level="WARNING")
    ql = _pipeline_ql(window_length=0.6, user_set=True)
    try:
        ql.run_single_event_search()
    finally:
        logger.remove(sink)
    r = ql.tls_results
    assert r["single_event_window"] == pytest.approx(0.6)
    assert r["single_event_window_reason"] == "user"
    assert r["single_event_max_duration"] == pytest.approx(0.2)
    # A 0.6 d window removes most of a 22 h transit before the search sees it;
    # the window is still respected, and the log says which durations are lost.
    assert any("skips longer trial durations" in str(m) for m in messages)


def test_pipeline_treats_locor_zero_window_as_unset():
    # LOCoR leaves window_length at 0 d with the rotation-aware window off.
    ql = _pipeline_ql(window_length=0.0, user_set=True)
    ql.run_single_event_search()
    r = ql.tls_results
    assert r["single_event_window_reason"] == "floor"
    assert r["single_event_window"] == pytest.approx(window_floor(1.0))


def test_pipeline_without_events_does_not_redetrend():
    ql = _pipeline_ql()
    t = ql.raw_lc.time.value
    flux = ql.raw_lc.flux.value.copy()
    flux[np.abs(t - 13.0) < 11 / 24] += 0.004  # remove the transit
    ql.raw_lc.flux = flux
    ql.run_single_event_search()
    assert ql.tls_results["single_event_redetrended"] is False
    assert ql.tls_results["event_class"] == "none"


def test_mono_panel_one_shows_the_search_trend():
    import matplotlib.pyplot as pl

    ql = _pipeline_ql()
    ql.run_single_event_search()
    fig, axes = pl.subplots(3, 3)
    ql._plot_single_event_panels(axes.flatten())
    labels = [line.get_label() for line in axes.flatten()[0].lines]
    assert any(label.startswith("single-event trend") for label in labels)
    pl.close(fig)


def test_tls_summary_has_single_event_window_column():
    from pathlib import Path

    from quicklook.app.tls_summary import extract_fields

    data = {"single_event_window": 2.9, "single_event_window_reason": "floor"}
    fields = extract_fields(data, "x_s01_pdcsap_sc_tls")
    assert fields["se_window"] == pytest.approx(2.9)
    assert fields["se_window_reason"] == "floor"
    html = (Path(__file__).parents[1] / "quicklook/app/templates/tls_summary.html").read_text()
    assert 'data-key="se_window"' in html
    assert html.count("<th data-key=") == html.count('<th><input type="text" data-col=')
