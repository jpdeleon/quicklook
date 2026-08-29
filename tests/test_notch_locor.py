"""Tests for quicklook.notch_locor (Notch filter and LOCoR detrending)."""

import numpy as np
import pytest

from quicklook.notch_locor import _label_cycles, locor_flatten, notch_flatten


# --- notch_flatten -----------------------------------------------------------


@pytest.fixture
def quiet_lc():
    rng = np.random.default_rng(0)
    t = np.arange(0, 1.0, 2 / 60 / 24)  # 1 day, 2-min cadence
    trend = 1.0 + 0.005 * np.sin(2 * np.pi * t / 0.7)
    flux = trend * (1 + rng.normal(0, 0.0006, len(t)))
    return t, flux


def test_notch_flatten_output_shapes(quiet_lc):
    t, flux = quiet_lc
    flat, trend = notch_flatten(t, flux, window_length=0.3)
    assert flat.shape == t.shape
    assert trend.shape == t.shape
    assert np.isfinite(flat).all()
    assert np.isfinite(trend).all()


def test_notch_flatten_reduces_scatter_from_slow_trend(quiet_lc):
    t, flux = quiet_lc
    flat, _ = notch_flatten(t, flux, window_length=0.3)
    assert np.std(flat) < np.std(flux)
    assert np.median(flat) == pytest.approx(1.0, abs=0.01)


def test_notch_flatten_preserves_injected_transit_depth(quiet_lc):
    t, flux = quiet_lc
    t0, dur, depth = 0.5, 0.08, 0.015
    in_transit = np.abs(t - t0) < dur / 2
    injected = flux.copy()
    injected[in_transit] *= 1 - depth

    flat, _ = notch_flatten(t, injected, window_length=0.3)
    measured_depth = 1 - np.median(flat[in_transit])
    # a windowed local detrender attenuates depth somewhat; require it to
    # survive in the right ballpark rather than exact recovery.
    assert 0.3 * depth < measured_depth < 3 * depth


def test_notch_flatten_handles_gaps_and_short_arrays():
    rng = np.random.default_rng(1)
    t = np.concatenate([np.arange(0, 0.3, 2 / 60 / 24), np.arange(0.6, 0.9, 2 / 60 / 24)])
    flux = 1 + rng.normal(0, 0.001, len(t))
    flat, trend = notch_flatten(t, flux, window_length=0.3)
    assert np.isfinite(flat).all() and np.isfinite(trend).all()

    t_short = t[:10]
    flux_short = flux[:10]
    flat_short, trend_short = notch_flatten(t_short, flux_short, window_length=0.3)
    assert np.isfinite(flat_short).all() and np.isfinite(trend_short).all()


def test_notch_flatten_min_delta_bic_inf_disables_notch(quiet_lc):
    """With min_delta_bic=inf, the notch model can never win over the null."""
    t, flux = quiet_lc
    t0, dur, depth = 0.5, 0.08, 0.03
    in_transit = np.abs(t - t0) < dur / 2
    injected = flux.copy()
    injected[in_transit] *= 1 - depth

    flat, _ = notch_flatten(t, injected, window_length=0.3, min_delta_bic=np.inf)
    # transit gets divided out (attenuated) by the plain quadratic trend
    measured_depth = 1 - np.median(flat[in_transit])
    assert measured_depth < depth


# --- LOCoR helpers -----------------------------------------------------------


def test_label_cycles_assigns_trailing_partial_cycle_its_own_id():
    """Regression check for the upstream off-by-one: the last (possibly
    partial) rotation cycle must get its own label, not silently merge into
    cycle 0."""
    phase = np.concatenate([np.linspace(0, 0.9, 5), np.linspace(0, 0.9, 5), np.linspace(0, 0.4, 3)])
    labels = _label_cycles(phase)
    assert labels[:5].tolist() == [0] * 5
    assert labels[5:10].tolist() == [1] * 5
    assert labels[10:].tolist() == [2] * 3


# --- locor_flatten -------------------------------------------------------------


@pytest.fixture
def spotted_lc():
    rng = np.random.default_rng(2)
    prot = 0.6
    t = np.arange(0, 6.0, 2 / 60 / 24)  # 6 days, 2-min cadence -> 10 cycles
    trend = 1.0 + 0.02 * np.sin(2 * np.pi * t / prot) + 0.004 * np.sin(2 * np.pi * t / prot * 2)
    flux = trend * (1 + rng.normal(0, 0.0005, len(t)))
    return t, flux, prot


def test_locor_flatten_output_shapes(spotted_lc):
    t, flux, prot = spotted_lc
    flat, trend = locor_flatten(t, flux, period=prot)
    assert flat.shape == t.shape
    assert trend.shape == t.shape
    assert np.isfinite(flat).all()
    assert np.isfinite(trend).all()


def test_locor_flatten_reduces_rotation_scatter(spotted_lc):
    t, flux, prot = spotted_lc
    flat, _ = locor_flatten(t, flux, period=prot)
    assert np.std(flat) < np.std(flux)


def test_locor_flatten_preserves_injected_transit_depth(spotted_lc):
    t, flux, prot = spotted_lc
    t0, dur, depth = 3.1, 0.1, 0.01
    in_transit = np.abs(t - t0) < dur / 2
    injected = flux.copy()
    injected[in_transit] *= 1 - depth

    flat, _ = locor_flatten(t, injected, period=prot)
    measured_depth = 1 - np.median(flat[in_transit])
    assert 0.3 * depth < measured_depth < 3 * depth


def test_locor_flatten_aliases_short_period(spotted_lc):
    """A period shorter than alias_num gets multiplied up so each cycle has
    enough points for a well-conditioned fit, instead of crashing."""
    t, flux, _ = spotted_lc
    short_period = 0.05  # shorter than the default alias_num window
    flat, trend = locor_flatten(t, flux, period=short_period, alias_num=0.2)
    assert np.isfinite(flat).all() and np.isfinite(trend).all()


def test_locor_flatten_rejects_nonpositive_period(spotted_lc):
    t, flux, _ = spotted_lc
    with pytest.raises(ValueError):
        locor_flatten(t, flux, period=0.0)
    with pytest.raises(ValueError):
        locor_flatten(t, flux, period=-1.0)
