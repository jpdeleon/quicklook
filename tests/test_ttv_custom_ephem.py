"""--custom-ephem drives the TTV timing, whatever TLS locked onto.

Young, spotted stars make TLS lock onto the rotation period, so a TTV search
from the TLS ephemeris times starspots. A user-supplied ephemeris must win.
Offline: the pipeline object is built with ``__new__`` and a synthetic light curve.
"""

import os

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
from lightkurve import LightCurve

from quicklook.tql import TessQuickLook, _TLSResult
from quicklook.ttv import trapezoid

CADENCE = 10 / 1440
PERIOD, T0, DUR, DEPTH, NOISE = 3.0, 1.0, 0.12, 3e-3, 5e-4
N_TRANSITS = 12
WRONG_PERIOD = 2.2  # what TLS "found" (e.g. the rotation period)


def _lc(seed=1):
    rng = np.random.default_rng(seed)
    t = np.arange(0.0, T0 + N_TRANSITS * PERIOD, CADENCE)
    f = 1 + rng.normal(0, NOISE, t.size)
    for k in range(N_TRANSITS):
        f -= DEPTH * trapezoid(t, T0 + k * PERIOD, DUR)
    return LightCurve(time=t, flux=f, flux_err=np.full_like(f, NOISE))


def _ql(ephem_source="custom", sde=20.0, mask_ephem=False):
    ql = TessQuickLook.__new__(TessQuickLook)
    ql.verbose = False
    ql.ttv_search = True
    ql.single_event_search = False
    ql.mask_ephem = mask_ephem
    ql.flat_lc = ql.raw_lc = _lc()
    ql.exofop_data = {}
    ql.ephem_source = ephem_source
    ql.toi_period = (PERIOD, 1e-5)
    ql.toi_epoch = (T0, 1e-3)
    ql.toi_dur = (DUR, 0.01)
    ql.tls_results = _TLSResult(
        {
            "SDE": sde,
            "period": WRONG_PERIOD,
            "T0": 0.5,
            "duration": 0.25,
            "distinct_transit_count": 20,
            "depth_variance_ratio": 1.0,
        }
    )
    return ql


def test_custom_ephem_overrides_a_periodic_tls_signal():
    ql = _ql()
    stats = ql.run_ttv_search()
    assert ql.ttv_source == "custom"
    assert stats.n == N_TRANSITS
    assert abs(stats.ephemeris.period - PERIOD) < 1e-3
    assert not stats.detected


def test_custom_ephem_times_transits_when_tls_is_not_periodic():
    ql = _ql(sde=4.0)
    stats = ql.run_ttv_search()
    assert ql.ttv_source == "custom"
    assert stats is not None and stats.n == N_TRANSITS


def test_catalogue_ephem_does_not_override_tls():
    ql = _ql(ephem_source="TFOP")
    ql.run_ttv_search()
    assert ql.ttv_source == "tls"


def test_masked_custom_ephem_is_not_timed():
    # --mask-ephem cuts those transits from the light curve, so TLS keeps the lead.
    ql = _ql(mask_ephem=True)
    ql.run_ttv_search()
    assert ql.ttv_source == "tls"


def test_diagnostics_use_the_custom_ephemeris():
    ql = _ql()
    ql.run_ttv_search()
    period, t0, duration, source = ql._diagnostic_ephemeris()
    assert source == "custom"
    assert abs(period - PERIOD) < 1e-3


def test_diagnostics_fall_back_to_custom_when_too_few_transits_are_timed():
    ql = _ql()
    ql.flat_lc = ql.raw_lc = _lc()[: int((T0 + 1.5 * PERIOD) / CADENCE)]  # two transits only
    ql.run_ttv_search()
    period, t0, duration, source = ql._diagnostic_ephemeris()
    assert source == "custom"
    assert period == PERIOD and t0 == T0


def test_summary_row_names_the_custom_source():
    ql = _ql()
    ql.run_ttv_search()
    rows = dict(ql._ttv_candidate_rows())
    assert rows["TTV source"] == "custom ephemeris"
