"""Single-event (mono-transit) search: detector, period lower limit, classification.

Offline: all light curves are synthetic.
"""

from pathlib import Path

import numpy as np
import pytest

from quicklook.single_event import (
    SingleEventResult,
    classify_events,
    duration_grid,
    min_period,
    single_event_search,
)

CADENCE = 120 / 86400  # 2-min


def _lc(
    span=27.0,
    events=((14.0, 0.25, 0.005),),
    noise=5e-4,
    gaps=((13.0 - 0.5, 13.0 + 0.5),),
    step=None,
    seed=1,
):
    """Synthetic flattened light curve with box dips, optional gaps and a step."""
    rng = np.random.default_rng(seed)
    t = np.arange(0.0, span, CADENCE)
    for lo, hi in gaps:
        t = t[(t < lo) | (t > hi)]
    f = 1.0 + rng.normal(0, noise, t.size)
    for t0, dur, depth in events:
        f[np.abs(t - t0) < dur / 2] -= depth
    if step is not None:
        t_step, height = step
        f[t > t_step] -= height
    return t, f


# --- duration grid -----------------------------------------------------------


def test_duration_grid_is_log_spaced_and_bounded():
    grid = duration_grid(span=27.0, cadence=CADENCE)
    assert np.all(np.diff(grid) > 0)
    assert grid[0] >= 3 * CADENCE
    assert grid[-1] <= 1.5  # days
    ratios = grid[1:] / grid[:-1]
    assert np.allclose(ratios, ratios[0])


def test_duration_grid_widens_for_low_density_star():
    sun = duration_grid(span=27.0, cadence=CADENCE, rho_star=1.0)
    giant = duration_grid(span=27.0, cadence=CADENCE, rho_star=0.01)
    assert giant[-1] > sun[-1]


# --- detector ----------------------------------------------------------------


def test_recovers_injected_single_transit():
    t, f = _lc()
    res = single_event_search(t, f)
    assert isinstance(res, SingleEventResult)
    assert len(res.events) == 1
    ev = res.events[0]
    assert ev.t0 == pytest.approx(14.0, abs=0.05)
    assert ev.depth == pytest.approx(0.005, rel=0.25)
    assert ev.duration == pytest.approx(0.25, rel=0.5)
    assert ev.snr > 20
    assert ev.delta_bic_systematics > 0


def test_real_toi2180b_transit_at_sector_start_is_recovered():
    # Regression: TOI-2180 b (P=260 d, ~3 ppt) in TESS S48. In both the raw and
    # flattened flux the dip runs 2610.75-2611.68 BTJD (centre ~2611.22, ~22 h),
    # starting 0.25 d after the data begin. A 2x-duration baseline requirement
    # and an unrefined transit model let a step model win. (The ExoFOP T0,
    # 2611.45, is ~5 h later than the dip seen here.)
    data = np.load(Path(__file__).parent / "data" / "toi2180_s48_flat_cutout.npz")
    res = single_event_search(data["time"], data["flux"].astype(float), rho_star=0.3)
    assert len(res.events) >= 1
    ev = res.events[0]
    assert ev.t0 == pytest.approx(2611.22, abs=0.08)
    assert 0.75 < ev.duration < 1.1
    assert ev.depth == pytest.approx(0.0035, rel=0.3)
    assert ev.delta_bic_systematics > 10


def test_long_transit_near_data_start_is_found_with_its_full_duration():
    # Regression (TOI-2180 b, TESS S48): a ~24 h, 3 ppt transit whose ingress is
    # 0.44 d after the data start, with a +2 ppt detrending edge ramp before it.
    t, f = _lc(span=25.0, events=((0.94 + 0.5, 1.0, 0.003),), gaps=((11.0, 15.0),))
    f = f + 0.002 * np.exp(-t / 0.3)
    res = single_event_search(t, f, rho_star=0.3)
    assert len(res.events) == 1
    ev = res.events[0]
    assert ev.t0 == pytest.approx(1.44, abs=0.1)
    assert ev.duration == pytest.approx(1.0, rel=0.3)
    assert ev.delta_bic_systematics > 10


def test_noise_only_has_no_events():
    t, f = _lc(events=())
    res = single_event_search(t, f)
    assert res.events == ()


@pytest.mark.parametrize("amp", [3e-4, 1e-3])
def test_red_noise_alone_has_no_events(amp):
    rng = np.random.default_rng(3)
    t = np.arange(0.0, 27.0, CADENCE)
    red = np.convolve(rng.normal(0, 1, t.size), np.ones(90) / 90, "same")
    f = 1 + rng.normal(0, 5e-4, t.size) + amp * red / red.std()
    assert single_event_search(t, f).events == ()


def test_transit_still_found_in_red_noise():
    rng = np.random.default_rng(4)
    t = np.arange(0.0, 27.0, CADENCE)
    red = np.convolve(rng.normal(0, 1, t.size), np.ones(90) / 90, "same")
    f = 1 + rng.normal(0, 5e-4, t.size) + 3e-4 * red / red.std()
    f[np.abs(t - 9.0) < 0.15] -= 0.004
    events = single_event_search(t, f).events
    assert len(events) >= 1
    assert events[0].t0 == pytest.approx(9.0, abs=0.05)


def test_step_discontinuity_is_rejected_as_systematic():
    t, f = _lc(events=(), step=(20.0, 0.004))
    res = single_event_search(t, f)
    assert all(abs(ev.t0 - 20.0) > 0.5 for ev in res.events)


def test_brightening_is_not_an_event():
    t, f = _lc(events=((14.0, 0.25, -0.005),))
    assert single_event_search(t, f).events == ()


def test_event_at_gap_edge_is_excluded():
    # Dip straddling the start of the mid-sector gap: only one side has baseline.
    t, f = _lc(events=((12.45, 0.25, 0.005),))
    res = single_event_search(t, f)
    assert all(abs(ev.t0 - 12.45) > 0.2 for ev in res.events)


def test_two_separate_events_are_both_reported_strongest_first():
    t, f = _lc(events=((5.0, 0.2, 0.004), (20.0, 0.2, 0.008)))
    res = single_event_search(t, f)
    assert len(res.events) == 2
    assert res.events[0].t0 == pytest.approx(20.0, abs=0.05)
    assert res.events[1].t0 == pytest.approx(5.0, abs=0.05)


def test_windowed_search_only_looks_inside_windows():
    t, f = _lc(events=((5.0, 0.2, 0.004), (20.0, 0.2, 0.008)))
    res = single_event_search(t, f, windows=[(4.5, 5.5)])
    assert len(res.events) == 1
    assert res.events[0].t0 == pytest.approx(5.0, abs=0.05)


def test_statistic_time_series_is_returned_for_plotting():
    t, f = _lc()
    res = single_event_search(t, f)
    assert res.grid_t0.shape == res.grid_snr.shape
    assert np.nanargmax(res.grid_snr) == np.argmin(np.abs(res.grid_t0 - res.events[0].t0))


def test_ignores_non_finite_points():
    t, f = _lc()
    f[::50] = np.nan
    res = single_event_search(t, f)
    assert len(res.events) == 1


def test_reports_distance_to_nearest_gap():
    t, f = _lc(events=((14.0, 0.25, 0.005),))
    ev = single_event_search(t, f).events[0]
    assert ev.dist_to_gap == pytest.approx(14.0 - 13.5, abs=0.05)


# --- period lower limit --------------------------------------------------------


def test_min_period_without_gaps_is_distance_to_far_edge():
    t, _ = _lc(gaps=())
    p = min_period(t, t0=10.0, duration=0.2)
    assert p == pytest.approx(27.0 - 10.0, abs=0.1)


def test_min_period_lets_second_transit_hide_in_gap():
    # Data cover 2.0-12.5 and 13.5-14.0. For t0=7.5, P <= 5.5 puts t0-P in the data;
    # just above 5.5, t0+P ~ 13 falls in the gap, so that is the lower limit even
    # though the far data edge (14.0) is 6.5 d away.
    t = np.arange(0.0, 27.0, CADENCE)
    t = t[(t > 2.0) & ((t < 12.5) | (t > 13.5)) & (t < 14.0)]
    p = min_period(t, t0=7.5, duration=0.2)
    assert p == pytest.approx(5.5, abs=0.15)


def test_min_period_never_below_covered_span_to_one_side():
    t, _ = _lc()
    p = min_period(t, t0=14.0, duration=0.25)
    assert p >= 27.0 - 14.0 - 0.2


# --- classification -------------------------------------------------------------


def test_classify_none_without_events_or_tls():
    assert classify_events(n_events=0, tls_periodic=False, n_unexplained=0) == "none"


def test_classify_mono_for_single_unexplained_event():
    assert classify_events(n_events=1, tls_periodic=False, n_unexplained=1) == "mono"


def test_classify_periodic_when_tls_explains_events():
    assert classify_events(n_events=3, tls_periodic=True, n_unexplained=0) == "periodic"


def test_classify_periodic_wins_even_with_extra_event():
    assert classify_events(n_events=4, tls_periodic=True, n_unexplained=1) == "periodic"


def test_short_transit_found_in_multi_year_gappy_light_curve():
    # Duration grid must follow the longest contiguous segment, not the total
    # span: a 2000-day stitched light curve still contains 2 h transits.
    rng = np.random.default_rng(5)
    t = np.concatenate([np.arange(s, s + 27.0, CADENCE) for s in range(0, 2000, 350)])
    f = 1 + rng.normal(0, 5e-4, t.size)
    f[np.abs(t - 1060.0) < 1 / 24] -= 0.004
    res = single_event_search(t, f)
    assert res.durations[0] < 2 / 24
    assert len(res.events) == 1
    assert res.events[0].t0 == pytest.approx(1060.0, abs=0.02)
