"""Unit tests for quicklook.ttv (per-transit timing, O-C and event chaining)."""

from __future__ import annotations

import numpy as np
import pytest

from quicklook.single_event import SingleEvent
from quicklook.ttv import (
    TTV_PVALUE_MAX,
    TransitTime,
    chain_events,
    cycle_count_solutions,
    fit_transit_time,
    fit_transit_times,
    linear_ephemeris,
    trapezoid,
    ttv_mask,
    ttv_statistics,
)

CADENCE = 10 / 1440
PERIOD, T0, DUR, DEPTH, NOISE = 3.0, 1.0, 0.12, 3e-3, 5e-4


def _light_curve(shifts=None, n_transits=15, gap=None, seed=1, durations=None):
    """Flat light curve with trapezoid transits at ``T0 + k*P + shift_k``."""
    rng = np.random.default_rng(seed)
    t = np.arange(0.0, T0 + n_transits * PERIOD, CADENCE)
    if gap is not None:
        t = t[(t < gap[0]) | (t > gap[1])]
    f = 1 + rng.normal(0, NOISE, t.size)
    shifts = np.zeros(n_transits) if shifts is None else np.asarray(shifts)
    durations = np.full(n_transits, DUR) if durations is None else np.asarray(durations)
    for k in range(n_transits):
        f -= DEPTH * trapezoid(t, T0 + k * PERIOD + shifts[k], durations[k])
    return t, f


def _sinusoid(n=15, amplitude=0.02, period_epochs=7.0):
    return amplitude * np.sin(2 * np.pi * np.arange(n) / period_epochs)


# --- model --------------------------------------------------------------------------


def test_trapezoid_is_one_at_centre_and_zero_outside():
    t = np.array([0.0, 0.2, -0.2, 0.059])
    shape = trapezoid(t, 0.0, 0.12)
    assert shape[0] == pytest.approx(1.0)
    assert shape[1] == shape[2] == 0.0
    assert 0 < shape[3] < 1  # in ingress


# --- single transit fits ------------------------------------------------------------


def test_fit_transit_time_recovers_a_shifted_transit():
    t, f = _light_curve(shifts=np.r_[0.03, np.zeros(14)])
    tt = fit_transit_time(t, f, T0, DUR, delta=0.1, epoch=0)
    assert tt is not None
    assert tt.tc == pytest.approx(T0 + 0.03, abs=3 * tt.tc_err + 1e-3)
    assert tt.tc_err < 0.005
    assert tt.depth == pytest.approx(DEPTH, rel=0.2)
    assert tt.duration == pytest.approx(DUR, rel=0.2)


def test_fit_transit_time_returns_none_in_a_data_gap():
    t, f = _light_curve(gap=(T0 + 3 * PERIOD - 0.5, T0 + 3 * PERIOD + 0.5))
    assert fit_transit_time(t, f, T0 + 3 * PERIOD, DUR, delta=0.1, epoch=3) is None


def test_fit_transit_time_returns_none_without_a_dip():
    rng = np.random.default_rng(3)
    t = np.arange(0, 5, CADENCE)
    f = 1 + rng.normal(0, NOISE, t.size)
    assert fit_transit_time(t, f, 2.5, DUR, delta=0.1, epoch=0) is None


def test_fit_transit_times_skips_gaps_and_keeps_epochs():
    t, f = _light_curve(gap=(T0 + 5 * PERIOD - 0.5, T0 + 5 * PERIOD + 0.5))
    times = fit_transit_times(t, f, PERIOD, T0, DUR)
    epochs = [tt.epoch for tt in times]
    assert 5 not in epochs
    assert len(times) == 14


def test_fit_transit_times_measures_each_duration():
    durations = np.where(np.arange(15) < 8, 0.10, 0.16)
    t, f = _light_curve(durations=durations)
    times = fit_transit_times(t, f, PERIOD, T0, 0.13)
    early = np.median([tt.duration for tt in times if tt.epoch < 8])
    late = np.median([tt.duration for tt in times if tt.epoch >= 8])
    assert early == pytest.approx(0.10, rel=0.15)
    assert late == pytest.approx(0.16, rel=0.15)


# --- ephemeris and O-C --------------------------------------------------------------


def _times(epochs, tcs, err=0.001):
    return tuple(
        TransitTime(int(e), float(tc), err, DUR, DEPTH, 10.0) for e, tc in zip(epochs, tcs)
    )


def test_linear_ephemeris_recovers_period_and_t0():
    epochs = np.arange(10)
    eph = linear_ephemeris(_times(epochs, T0 + PERIOD * epochs))
    assert eph.period == pytest.approx(PERIOD)
    assert eph.t0 == pytest.approx(T0)
    assert eph.chi2 == pytest.approx(0.0, abs=1e-9)
    assert eph.dof == 8


def test_linear_ephemeris_needs_two_epochs():
    with pytest.raises(ValueError):
        linear_ephemeris(_times([0, 0], [T0, T0]))


def test_no_ttv_gives_a_high_pvalue():
    rng = np.random.default_rng(5)
    epochs = np.arange(12)
    tcs = T0 + PERIOD * epochs + rng.normal(0, 0.001, epochs.size)
    stats = ttv_statistics(_times(epochs, tcs))
    assert stats.pvalue > TTV_PVALUE_MAX
    assert not stats.detected


def test_sinusoidal_ttv_is_detected_with_its_amplitude():
    epochs = np.arange(15)
    tcs = T0 + PERIOD * epochs + _sinusoid()
    stats = ttv_statistics(_times(epochs, tcs))
    assert stats.detected
    assert stats.pvalue < 1e-6
    assert stats.amplitude == pytest.approx(0.02, rel=0.15)
    assert stats.oc.size == 15


def test_ttv_statistics_needs_three_transits():
    stats = ttv_statistics(_times([0, 1], [T0, T0 + PERIOD]))
    assert not stats.detected
    assert np.isnan(stats.pvalue)


def test_end_to_end_ttv_from_light_curve():
    t, f = _light_curve(shifts=_sinusoid())
    times = fit_transit_times(t, f, PERIOD, T0, DUR)
    stats = ttv_statistics(times)
    assert len(times) == 15
    assert stats.detected
    assert stats.amplitude == pytest.approx(0.02, rel=0.3)


def test_end_to_end_no_ttv_from_light_curve():
    t, f = _light_curve()
    stats = ttv_statistics(fit_transit_times(t, f, PERIOD, T0, DUR))
    assert not stats.detected


def test_fit_follows_ttvs_larger_than_half_a_duration():
    # 0.08 d shifts exceed half the 0.12 d duration: a linear-ephemeris mask misses them.
    shifts = 0.08 * np.where(np.arange(15) % 2, 1, -1)
    t, f = _light_curve(shifts=shifts)
    times = fit_transit_times(t, f, PERIOD, T0, DUR)
    fitted = np.array([tt.tc for tt in times])
    expected = (
        T0 + PERIOD * np.array([tt.epoch for tt in times]) + shifts[[tt.epoch for tt in times]]
    )
    assert np.allclose(fitted, expected, atol=0.01)


# --- masking ------------------------------------------------------------------------


def test_ttv_mask_follows_fitted_times():
    t = np.arange(0, 10, CADENCE)
    times = _times([0, 1], [1.08, 4.0])
    mask = ttv_mask(t, times, margin=1.0)
    assert mask[np.argmin(np.abs(t - 1.08))]
    assert not mask[np.argmin(np.abs(t - 1.0 + 0.07))]  # linear time 1.0 - margin is outside
    assert mask[np.argmin(np.abs(t - 4.0))]


# --- cycle-count ambiguity ----------------------------------------------------------


def test_cycle_count_is_unique_with_dense_coverage():
    epochs = np.arange(10)
    sols = cycle_count_solutions(_times(epochs, T0 + PERIOD * epochs))
    assert len(sols) == 1
    assert sols[0].period == pytest.approx(PERIOD)


def test_cycle_count_ambiguity_across_a_long_gap():
    # Two short segments 600 d apart with loose timings: N or N+-1 cycles all fit.
    p = 10.0
    epochs = np.r_[0, 1, 60, 61]
    tcs = T0 + p * epochs
    sols = cycle_count_solutions(_times(epochs, tcs, err=0.1), gap_min=50.0)
    periods = sorted(s.period for s in sols)
    assert len(sols) >= 3
    assert any(np.isclose(per, p) for per in periods)
    assert sols[0].shifts == (0,)


# --- chaining single events ---------------------------------------------------------


def _event(t0, depth=DEPTH, duration=DUR, snr=12.0):
    return SingleEvent(t0, duration, depth, snr, 50.0, 30.0, 1.0, 20)


def test_chain_events_finds_the_period_with_ttvs():
    events = [_event(T0 + 20 * k + s) for k, s in zip([0, 1, 3, 4], [0.0, 0.03, -0.04, 0.02])]
    chain = chain_events(events)
    assert chain is not None
    period, t0, members = chain
    assert period == pytest.approx(20.0, abs=0.05)
    assert len(members) == 4


def test_chain_events_needs_three_consistent_events():
    assert chain_events([_event(1.0), _event(21.0)]) is None
    # Third event twice as deep: not the same planet.
    assert chain_events([_event(1.0), _event(21.0), _event(41.0, depth=2 * DEPTH)]) is None


def test_chain_events_prefers_the_longest_period_that_fits():
    events = [_event(T0 + 30 * k) for k in range(3)]
    period, _, _ = chain_events(events)
    assert period == pytest.approx(30.0, abs=0.01)


def test_a_single_outlier_does_not_make_a_ttv_detection():
    rng = np.random.default_rng(7)
    epochs = np.arange(12)
    tcs = T0 + PERIOD * epochs + rng.normal(0, 0.001, epochs.size)
    tcs[6] += 0.02  # one 20-sigma timing, e.g. a flare in ingress
    stats = ttv_statistics(_times(epochs, tcs))
    assert stats.pvalue < TTV_PVALUE_MAX  # the full chi-square is fooled ...
    assert stats.pvalue_robust > TTV_PVALUE_MAX  # ... dropping the worst transit is not
    assert not stats.detected


def test_consecutive_long_period_transits_are_not_ambiguous():
    # P = 34.5 d: consecutive transits are > 30 d apart but no epoch is missing,
    # so shifting each by a cycle (giving P/2) must not count as an alternative.
    p = 34.5
    epochs = np.arange(5)
    sols = cycle_count_solutions(_times(epochs, T0 + p * epochs, err=0.002))
    assert len(sols) == 1
    assert sols[0].period == pytest.approx(p)
