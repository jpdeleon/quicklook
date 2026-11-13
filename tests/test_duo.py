"""Duotransit pairing and period aliases. Offline, synthetic light curves."""

import numpy as np
import pytest

from quicklook.duo import (
    DuoPair,
    duo_aliases,
    duration_likelihood,
    events_consistent,
    find_duo,
)
from quicklook.single_event import SingleEvent, classify_events, single_event_search

CADENCE = 120 / 86400
SECTORS = ((0.0, 27.0), (700.0, 727.0))  # e.g. TESS year 1 and year 3


def _time(sectors=SECTORS):
    return np.concatenate([np.arange(lo, hi, CADENCE) for lo, hi in sectors])


def _event(t0, depth=0.004, duration=0.3, snr=30.0):
    return SingleEvent(
        t0=t0,
        duration=duration,
        depth=depth,
        snr=snr,
        delta_bic_null=200.0,
        delta_bic_systematics=150.0,
        dist_to_gap=5.0,
        n_in=200,
    )


# --- consistency -----------------------------------------------------------------


def test_matching_events_are_consistent():
    assert events_consistent(_event(10.0), _event(710.0, depth=0.0042, duration=0.32))


def test_depth_mismatch_is_inconsistent():
    assert not events_consistent(_event(10.0, depth=0.004), _event(710.0, depth=0.008))


def test_duration_mismatch_is_inconsistent():
    assert not events_consistent(_event(10.0, duration=0.2), _event(710.0, duration=0.45))


def test_small_depth_difference_at_high_snr_is_tolerated():
    # Box fits bias depths by a few percent, which exceeds depth/SNR at SNR ~ 100.
    a, b = _event(10.0, depth=0.00396, snr=95.0), _event(710.0, depth=0.00372, snr=92.0)
    assert events_consistent(a, b)


def test_low_snr_depth_difference_is_tolerated():
    a, b = _event(10.0, depth=0.004, snr=8.0), _event(710.0, depth=0.0055, snr=8.0)
    assert events_consistent(a, b)


# --- aliases -------------------------------------------------------------------------


def test_aliases_are_delta_t_over_n_and_respect_coverage():
    aliases = duo_aliases(_time(), _event(10.0), _event(710.0))
    periods = np.array([a.period for a in aliases])
    ns = np.array([a.n for a in aliases])
    assert np.allclose(periods, 700.0 / ns)
    assert 700.0 in periods  # n=1
    # Shorter than ~17 d puts another transit inside sector 1 (10 + P < 27).
    assert periods.min() >= 17.0 - 0.5
    assert len(aliases) == len(set(ns))


def test_alias_with_extra_transit_in_data_is_excluded():
    # A third sector covering 360 d excludes n=2 (P=350 -> transit at 360).
    time = _time(SECTORS + ((350.0, 377.0),))
    ns = {a.n for a in duo_aliases(time, _event(10.0), _event(710.0))}
    assert 2 not in ns
    assert 1 in ns


def test_alias_weights_are_normalised_and_disfavour_too_long_durations():
    aliases = duo_aliases(_time(), _event(10.0, duration=0.3), _event(710.0, duration=0.3))
    weights = {a.n: a.weight for a in aliases}
    assert sum(weights.values()) == pytest.approx(1.0)
    # P=700/41~17 d: a 7 h transit is far longer than a central transit there.
    shortest = max(weights)
    assert weights[shortest] < max(weights.values()) / 100


def test_duration_likelihood_peaks_below_central_duration_and_has_a_tail():
    assert duration_likelihood(0.2, 0.3) > duration_likelihood(0.05, 0.3)
    assert duration_likelihood(0.6, 0.3) < duration_likelihood(0.29, 0.3) / 100
    assert duration_likelihood(0.6, 0.3) > 0


# --- pairing ---------------------------------------------------------------------------


def test_find_duo_pairs_consistent_events():
    pair = find_duo(_time(), (_event(10.0), _event(710.0)))
    assert isinstance(pair, DuoPair)
    assert pair.delta_t == pytest.approx(700.0)
    assert pair.first.t0 < pair.second.t0
    assert pair.aliases


def test_find_duo_rejects_inconsistent_events():
    assert find_duo(_time(), (_event(10.0, depth=0.002), _event(710.0, depth=0.01))) is None


def test_find_duo_needs_two_events():
    assert find_duo(_time(), (_event(10.0),)) is None


def test_find_duo_skips_pairs_with_no_allowed_alias():
    # 5 d apart inside one sector: every P = 5/n puts a third transit in the data.
    time = _time(((0.0, 27.0),))
    assert find_duo(time, (_event(10.0), _event(15.0))) is None


def test_classify_duo():
    assert classify_events(n_events=2, tls_periodic=False, n_unexplained=2, has_duo=True) == "duo"
    assert classify_events(n_events=2, tls_periodic=True, n_unexplained=2, has_duo=True) == (
        "periodic"
    )


# --- end to end on a two-sector light curve ----------------------------------------------


def test_detector_and_pairing_on_two_sector_light_curve():
    rng = np.random.default_rng(2)
    t = _time()
    f = 1 + rng.normal(0, 5e-4, t.size)
    for t0 in (10.0, 710.0):
        f[np.abs(t - t0) < 0.15] -= 0.004
    events = single_event_search(t, f).events
    assert len(events) == 2
    pair = find_duo(t, events)
    assert pair is not None
    assert pair.delta_t == pytest.approx(700.0, abs=0.1)
