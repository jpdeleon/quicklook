"""TTV vetting from the young-planet validation run (issue #42, step 6).

- Red-noise beta inflates white-noise timing errors (hot Jupiters flagged with
  O-C chi2/dof 11-41 in the TGLC run).
- A detection counts only on a catalogue period and away from the rotation
  period (HIP 67522, Kepler-616 and TIC 434398831 timed stellar variability).
"""

import numpy as np
import pytest

from quicklook.ttv import (
    BETA_MAX,
    TransitTime,
    inflate_errors,
    red_noise_beta,
)

PERIOD, T0, DUR = 3.0, 1.0, 0.12
CADENCE = 10 / 1440


def _time(span=27.0):
    return np.arange(0.0, span, CADENCE)


@pytest.mark.unit
def test_beta_is_one_for_white_noise():
    t = _time()
    f = 1 + np.random.default_rng(0).normal(0, 1e-3, t.size)

    assert red_noise_beta(t, f, PERIOD, T0, DUR) == pytest.approx(1.0, abs=0.15)


@pytest.mark.unit
def test_beta_grows_with_correlated_noise():
    rng = np.random.default_rng(1)
    t = _time()
    white = rng.normal(0, 1e-3, t.size)
    # Red noise on ~1 h timescales: smoothed random walk.
    red = np.convolve(rng.normal(0, 1e-3, t.size), np.ones(6) / np.sqrt(6), mode="same")

    beta = red_noise_beta(t, 1 + white + red, PERIOD, T0, DUR)

    assert 1.5 < beta <= BETA_MAX


@pytest.mark.unit
def test_beta_ignores_in_transit_cadences():
    from quicklook.ttv import trapezoid

    t = _time()
    f = 1 + np.random.default_rng(2).normal(0, 1e-3, t.size)
    for n in range(9):
        f -= 0.01 * trapezoid(t, T0 + n * PERIOD, DUR)  # deep transits, 10x the noise

    assert red_noise_beta(t, f, PERIOD, T0, DUR) == pytest.approx(1.0, abs=0.15)


@pytest.mark.unit
@pytest.mark.parametrize(
    "time,flux,period",
    [
        (np.arange(5.0), np.ones(5), PERIOD),  # too few cadences
        (_time(), np.ones(_time().size), PERIOD),  # zero scatter
        (_time(), np.ones(_time().size), np.nan),  # bad ephemeris
    ],
    ids=["short", "flat", "nan-period"],
)
def test_beta_defaults_to_one(time, flux, period):
    assert red_noise_beta(time, flux, period, T0, DUR) == 1.0


@pytest.mark.unit
def test_inflate_errors_scales_only_tc_err():
    times = (TransitTime(0, 1.0, 0.001, 0.1, 0.01, 10.0),)

    (out,) = inflate_errors(times, 2.5)

    assert out.tc_err == pytest.approx(0.0025)
    assert (out.epoch, out.tc, out.duration, out.depth, out.snr) == (0, 1.0, 0.1, 0.01, 10.0)
    assert times[0].tc_err == 0.001  # input untouched


# --- vetting ---------------------------------------------------------------------------


def _vetted_ql(period=PERIOD, detected=True, toi_period=None, planets=(), prot=None):
    from types import SimpleNamespace

    from quicklook.tql import TessQuickLook, _TLSResult

    ql = TessQuickLook.__new__(TessQuickLook)
    ql.verbose = False
    eph = SimpleNamespace(period=period)
    ql.ttv_stats = SimpleNamespace(ephemeris=eph, detected=detected, n=6)
    ql.toi_period = None if toi_period is None else (toi_period, 0.0)
    ql.exofop_data = {"planet_parameters": [{"per": str(p)} for p in planets]}
    ql.tls_results = _TLSResult({})
    if prot is not None:
        ql.Prot_ls = prot
    return ql


@pytest.mark.unit
def test_vetting_keeps_detection_on_catalogue_period():
    ql = _vetted_ql(period=2.7444, toi_period=2.7443, prot=9.0)

    assert ql.vet_ttv_detection() is True
    r = ql.tls_results
    assert (r["ttv_detected"], r["ttv_detected_raw"]) == (True, True)
    assert (r["ttv_period_match"], r["ttv_rotation_alias"]) == (True, False)


@pytest.mark.unit
def test_vetting_rejects_period_off_catalogue():
    # HIP 67522: timed 1.415 d, planets at 6.96 and 14.33 d.
    ql = _vetted_ql(period=1.4153, toi_period=6.9595, planets=(6.9595, 14.335))

    assert ql.vet_ttv_detection() is False
    assert ql.tls_results["ttv_detected_raw"] is True
    assert ql.tls_results["ttv_period_match"] is False
    assert ql._ttv_verdict_label() == "rejected (vetting)"


@pytest.mark.unit
def test_vetting_matches_any_exofop_planet():
    ql = _vetted_ql(period=14.334, toi_period=6.9595, planets=(6.9595, 14.335))

    assert ql.vet_ttv_detection() is True


@pytest.mark.unit
@pytest.mark.parametrize("prot", [1.4153, 2 * 1.4153, 1.4153 / 2])
def test_vetting_rejects_rotation_alias(prot):
    ql = _vetted_ql(period=1.4153, prot=prot)  # no catalogue period known

    assert ql.vet_ttv_detection() is False
    assert ql.tls_results["ttv_period_match"] is None
    assert ql.tls_results["ttv_rotation_alias"] is True


@pytest.mark.unit
def test_vetting_without_catalogue_or_rotation_keeps_raw_verdict():
    ql = _vetted_ql(period=4.2)

    assert ql.vet_ttv_detection() is True
    assert ql.tls_results["ttv_period_match"] is None
    assert ql.tls_results["ttv_rotation_alias"] is None


@pytest.mark.unit
def test_vetting_ignores_placeholder_periods():
    # TFOP fills missing periods with 0.1 d.
    ql = _vetted_ql(period=4.2, toi_period=0.1, planets=("", 0.1))

    assert ql.vet_ttv_detection() is True
    assert ql.tls_results["ttv_period_match"] is None


@pytest.mark.unit
def test_vetting_of_non_detection_stays_negative():
    ql = _vetted_ql(period=2.7444, detected=False, toi_period=2.7443)

    assert ql.vet_ttv_detection() is False
    assert ql._ttv_verdict_label() == "not significant"


@pytest.mark.unit
def test_vetting_prefers_gls_rotation_period():
    ql = _vetted_ql(period=1.4153)
    ql.tls_results["Prot_gls"] = (1.42, 0.01)
    ql.Prot_ls = 9.0

    assert ql.vet_ttv_detection() is False
    assert ql.tls_results["ttv_rotation_alias"] is True
