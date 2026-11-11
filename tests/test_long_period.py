"""Long-period support: extended TLS period ceiling and the catalogue-epoch check.

Offline: the TLS engine is replaced by a scripted stand-in, light curves are synthetic.
"""

import os

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import pytest

from quicklook import tql as tql_mod
from quicklook.long_period import (
    EPOCH_SNR_MIN,
    default_period_ceiling,
    epoch_check,
    needs_standard_fallback,
)
from quicklook.tql import TessQuickLook, _TLSResult


def _lc(period=15.0, t0=5.0, dur=0.2, depth=0.01, noise=0.001, span=27.0, cad=600 / 86400, seed=0):
    rng = np.random.default_rng(seed)
    t = np.arange(0, span, cad)
    f = 1 + rng.normal(0, noise, len(t))
    if depth:
        phase = ((t - t0 + period / 2) % period) - period / 2
        f[np.abs(phase) < dur / 2] -= depth
    return t, f


# --- ceiling ----------------------------------------------------------------


def test_default_ceiling_is_half_baseline_when_not_extended():
    assert default_period_ceiling(27.0, extend=False) == pytest.approx(13.5)


def test_default_ceiling_extends_to_nearly_full_baseline():
    assert default_period_ceiling(27.0, extend=True) == pytest.approx(0.98 * 27.0)


@pytest.mark.parametrize(
    "ntransits, period_max, span, expected",
    [
        (1, 26.0, 27.0, True),  # extended range, single-transit winner -> fall back
        (2, 26.0, 27.0, False),  # extended range but >= 2 transits -> keep
        (1, 13.0, 27.0, False),  # standard range: TLS already enforces >= 2
    ],
)
def test_needs_standard_fallback(ntransits, period_max, span, expected):
    assert needs_standard_fallback(ntransits, period_max, span) is expected


# --- epoch check ------------------------------------------------------------


def test_epoch_check_detects_dip_at_catalogue_epoch():
    t, f = _lc()
    res = epoch_check(t, f, period=15.0, epoch=5.0, duration=0.2)
    assert res["epoch_n_covered"] == 2
    assert res["epoch_snr"] > 10
    assert res["epoch_detected"] is True
    assert res["epoch_control_snr_p95"] < EPOCH_SNR_MIN


def test_epoch_check_rejects_noise_only():
    t, f = _lc(depth=0.0)
    res = epoch_check(t, f, period=15.0, epoch=5.0, duration=0.2)
    assert res["epoch_n_covered"] == 2
    assert res["epoch_detected"] is False


def test_epoch_check_wrong_epoch_not_detected():
    t, f = _lc()  # dips at 5 and 20; claim epoch 12
    res = epoch_check(t, f, period=15.0, epoch=12.0, duration=0.2)
    assert res["epoch_detected"] is False


def test_epoch_check_no_covered_transit():
    t, f = _lc(span=3.0)  # data 0-3 d; transits at 5, 20 -> nothing covered
    res = epoch_check(t, f, period=15.0, epoch=5.0, duration=0.2)
    assert res["epoch_n_covered"] == 0
    assert np.isnan(res["epoch_snr"])
    assert res["epoch_detected"] is False


def test_epoch_check_single_covered_transit():
    t, f = _lc(span=10.0)  # only the transit at t=5 is covered
    res = epoch_check(t, f, period=15.0, epoch=5.0, duration=0.2)
    assert res["epoch_n_covered"] == 1
    assert res["epoch_detected"] is True


def test_epoch_check_ignores_transit_in_gap():
    t, f = _lc()
    keep = np.abs(t - 5.0) > 1.0  # blank out the first transit
    res = epoch_check(t[keep], f[keep], period=15.0, epoch=5.0, duration=0.2)
    assert res["epoch_n_covered"] == 1


@pytest.mark.parametrize(
    "bad", [dict(period=0), dict(period=np.nan), dict(duration=0), dict(epoch=np.nan)]
)
def test_epoch_check_invalid_ephemeris_returns_empty(bad):
    t, f = _lc()
    kwargs = dict(period=15.0, epoch=5.0, duration=0.2) | bad
    res = epoch_check(t, f, **kwargs)
    assert res["epoch_n_covered"] == 0
    assert res["epoch_detected"] is False


# --- TLS search wiring ------------------------------------------------------


class _RecordingTLS:
    """Scripted stand-in for transitleastsquares; records power() kwargs."""

    queue = []
    calls = []

    def __init__(self, *args, **kwargs):
        pass

    def power(self, **kwargs):
        _RecordingTLS.calls.append(kwargs)
        return _RecordingTLS.queue.pop(0)


def _result(period, ntransits):
    return _TLSResult({"period": period, "SDE": 10.0, "distinct_transit_count": ntransits})


def _ql(monkeypatch, flat_lc, porb_max, extend=True):
    monkeypatch.setattr(tql_mod, "tls", _RecordingTLS)
    _RecordingTLS.calls = []
    ql = TessQuickLook.__new__(TessQuickLook)
    ql.verbose = False
    ql.tls_use_threads = None
    ql.use_star_priors = False
    ql.Porb_min = 0.5
    ql.Porb_max = porb_max
    ql.extend_period_ceiling = extend
    ql.flat_lc = flat_lc
    return ql


def test_extended_range_requests_single_transit_grid(monkeypatch, mock_light_curve):
    span = float(np.ptp(mock_light_curve.time.value))
    ql = _ql(monkeypatch, mock_light_curve, porb_max=0.98 * span)
    _RecordingTLS.queue = [_result(span * 0.6, 2)]
    res = ql._tls_search(mock_light_curve)
    assert res.distinct_transit_count == 2
    assert len(_RecordingTLS.calls) == 1
    assert _RecordingTLS.calls[0]["n_transits_min"] == 1
    assert _RecordingTLS.calls[0]["period_max"] == pytest.approx(0.98 * span)


def test_single_transit_winner_falls_back_to_standard_range(monkeypatch, mock_light_curve):
    span = float(np.ptp(mock_light_curve.time.value))
    ql = _ql(monkeypatch, mock_light_curve, porb_max=0.98 * span)
    _RecordingTLS.queue = [_result(span * 0.9, 1), _result(3.0, 8)]
    res = ql._tls_search(mock_light_curve)
    assert res.period == pytest.approx(3.0)
    assert len(_RecordingTLS.calls) == 2
    assert _RecordingTLS.calls[1]["period_max"] == pytest.approx(span / 2)
    assert "n_transits_min" not in _RecordingTLS.calls[1]


def test_standard_range_is_unchanged(monkeypatch, mock_light_curve):
    span = float(np.ptp(mock_light_curve.time.value))
    ql = _ql(monkeypatch, mock_light_curve, porb_max=span / 2, extend=False)
    _RecordingTLS.queue = [_result(3.0, 8)]
    ql._tls_search(mock_light_curve)
    assert len(_RecordingTLS.calls) == 1
    assert "n_transits_min" not in _RecordingTLS.calls[0]


def test_user_period_limits_below_half_baseline_are_respected(monkeypatch, mock_light_curve):
    span = float(np.ptp(mock_light_curve.time.value))
    ql = _ql(monkeypatch, mock_light_curve, porb_max=5.0)
    _RecordingTLS.queue = [_result(3.0, 8)]
    ql._tls_search(mock_light_curve)
    assert _RecordingTLS.calls[0]["period_max"] == pytest.approx(5.0)
    assert "n_transits_min" not in _RecordingTLS.calls[0]
    assert 5.0 < span / 2


# --- epoch check on the pipeline object / CLI wiring --------------------------


def _epoch_ql(period=15.0, epoch=5.0, dur=0.2, with_ephem=True):
    import astropy.units as u
    from lightkurve import LightCurve

    t, f = _lc(period=period, t0=epoch, dur=dur)
    ql = TessQuickLook.__new__(TessQuickLook)
    ql.verbose = False
    ql.ephem_source = "TFOP"
    ql.flat_lc = LightCurve(time=t, flux=f, flux_err=np.full_like(f, 1e-3))
    ql.tls_results = _TLSResult({"period": period})
    ql.toi_epoch = np.array((epoch, 0.01)) if with_ephem else None
    ql.toi_period = np.array((period, 0.001)) if with_ephem else None
    ql.toi_dur = np.array((dur, 0.01)) if with_ephem else None
    return ql


def test_run_epoch_check_stores_results_in_tls_results():
    ql = _epoch_ql()
    res = ql.run_epoch_check()
    assert res["epoch_detected"] is True
    for key in ("epoch_snr", "epoch_n_covered", "epoch_control_snr_p95", "epoch_detected"):
        assert key in ql.tls_results


def test_run_epoch_check_skipped_without_ephemeris():
    ql = _epoch_ql(with_ephem=False)
    assert ql.run_epoch_check() is None
    assert "epoch_snr" not in ql.tls_results


@pytest.mark.parametrize("flag, expected", [([], True), (["--no-extend-period-ceiling"], False)])
def test_cli_passes_extend_period_ceiling(flag, expected):
    from unittest.mock import patch
    from typer.testing import CliRunner
    from quicklook.cli.app import app

    with patch("quicklook.tql.TessQuickLook") as ql:
        CliRunner().invoke(app, ["run", "--name", "TOI-1234", "--sector", "5", "--save"] + flag)
    assert ql.call_args.kwargs["extend_period_ceiling"] is expected
