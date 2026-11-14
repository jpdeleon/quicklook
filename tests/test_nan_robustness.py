"""Robustness to NaN inputs seen in batch runs on TGLC HLSP light curves.

- Some TGLC HLSP files have all-NaN PSF flux (HATS-67 s10); quicklook must
  fall back to aperture photometry instead of crashing on an empty light curve.
- A degenerate TTV fit can yield a NaN period (CoRoT-18 s6); the local-baseline
  TTV diagnostics must skip it rather than raise on ``int(nan)``.
"""

from types import SimpleNamespace

import lightkurve as lk
import numpy as np
import pytest

from quicklook.tql import _select_tglc_flux
from quicklook.ttv_local import covered_epochs, run_local_ttv
from quicklook.ttv_pipeline import TTVMixin


def _tglc_lc(psf_finite: bool):
    n = 50
    psf = np.ones(n) if psf_finite else np.full(n, np.nan)
    lc = lk.LightCurve(time=np.arange(n, dtype=float), flux=psf)
    lc["cal_psf_flux"] = psf
    lc["cal_aper_flux"] = np.full(n, 2.0)
    lc.meta["AUTHOR"] = "TGLC"
    return lc


@pytest.mark.unit
def test_select_tglc_flux_uses_psf_when_finite():
    lc, used = _select_tglc_flux(_tglc_lc(psf_finite=True), "psf")

    assert used == "psf"
    assert np.allclose(lc.flux.value, 1.0)


@pytest.mark.unit
def test_select_tglc_flux_falls_back_to_aperture_when_psf_all_nan():
    lc, used = _select_tglc_flux(_tglc_lc(psf_finite=False), "psf")

    assert used == "aperture"
    assert np.allclose(lc.flux.value, 2.0)


@pytest.mark.unit
def test_select_tglc_flux_auto_falls_back_when_default_flux_all_nan():
    lc, used = _select_tglc_flux(_tglc_lc(psf_finite=False), "auto")

    assert used == "aperture"
    assert np.isfinite(lc.flux.value).all()


@pytest.mark.unit
@pytest.mark.parametrize("period,t0", [(np.nan, 1.0), (5.0, np.nan), (0.0, 1.0)])
def test_covered_epochs_rejects_invalid_ephemeris(period, t0):
    assert covered_epochs(np.linspace(0, 27, 1000), period, t0) == []


@pytest.mark.unit
def test_run_local_ttv_returns_none_for_nan_period():
    t = np.linspace(0, 27, 1000)

    assert run_local_ttv(t, np.ones_like(t), np.nan, 1.0, 0.1) is None


@pytest.mark.unit
def test_diagnostic_ephemeris_falls_back_to_toi_when_timed_period_nan():
    ql = SimpleNamespace(
        ttv_stats=SimpleNamespace(ephemeris=SimpleNamespace(period=np.nan, t0=1.0), n=6),
        ttv_times=[SimpleNamespace(duration=0.1)] * 6,
        ttv_source="tls",
        toi_period=(3.0, 0.0),
        toi_epoch=(2.0, 0.0),
        toi_dur=(0.1, 0.0),
        mask_ephem=False,
        ephem_source="toi",
    )

    period, t0, _, source = TTVMixin._diagnostic_ephemeris(ql)

    assert (period, t0, source) == (3.0, 2.0, "toi")


def _bare_ql(pipeline="tglc"):
    from quicklook.tql import TessQuickLook

    ql = object.__new__(TessQuickLook)
    ql.pipeline, ql.verbose, ql.sector = pipeline, False, 5
    return ql


@pytest.mark.unit
def test_load_tpf_returns_false_when_tesscut_and_spoc_fail(monkeypatch):
    from lightkurve.utils import LightkurveError

    from quicklook.tql import TessQuickLook

    def boom(*a, **k):
        raise LightkurveError("no data")

    monkeypatch.setattr(TessQuickLook, "get_tpf_tesscut", boom)
    monkeypatch.setattr(TessQuickLook, "get_tpf", boom)

    assert _bare_ql()._load_tpf() is False


@pytest.mark.unit
def test_load_tpf_falls_back_to_spoc(monkeypatch):
    from quicklook.tql import TessQuickLook

    def boom(*a, **k):
        raise TimeoutError("mast")

    monkeypatch.setattr(TessQuickLook, "get_tpf_tesscut", boom)
    monkeypatch.setattr(TessQuickLook, "get_tpf", lambda self, **k: "spoc-tpf")
    ql = _bare_ql()

    assert ql._load_tpf() is True
    assert (ql.tpf, ql.sap_mask) == ("spoc-tpf", "pipeline")
