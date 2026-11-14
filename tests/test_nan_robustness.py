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


@pytest.mark.unit
@pytest.mark.parametrize(
    "counts,depths,expected",
    [
        (np.nan, np.nan, 0),  # degenerate TLS fit returns scalars (HD 114082)
        ([10, 0, 8], [0.99, np.nan, 0.98], 2),
        ([10, 9], [0.99, 0.98, 0.97], 2),
    ],
)
def test_count_valid_transits_handles_scalars_and_arrays(counts, depths, expected):
    from quicklook.tql import _count_valid_transits

    assert _count_valid_transits(counts, depths) == expected


@pytest.mark.unit
def test_cycle_count_solutions_skips_singular_fits():
    # Infinite timing errors give zero weights and a singular normal matrix (WASP-93).
    from quicklook.ttv import TransitTime, cycle_count_solutions

    times = [TransitTime(n, 1.0 + 2.0 * n, np.inf, 0.1, 0.01, 10.0) for n in (0, 1, 2)]

    assert cycle_count_solutions(times) == []


@pytest.mark.unit
def test_odd_even_plot_with_no_cadences_near_transit():
    # Degenerate TLS fit: nothing falls inside the transit window (HD 114082).
    import matplotlib

    matplotlib.use("Agg")
    from quicklook.plot import plot_odd_even_transit
    from quicklook.tql import _TLSResult

    t = np.linspace(0, 27, 500)
    lc = lk.LightCurve(time=t, flux=np.ones_like(t))
    fold = lc.fold(period=10.0, epoch_time=5.0)
    fold = fold[np.abs(fold.time.value) > 2.0]  # empty the transit window
    tls = _TLSResult(
        {
            "duration": 0.1,
            "period": 10.0,
            "T0": 5.0,
            "depth": 0.99,
            "model_folded_phase": np.linspace(0, 1, 50),
            "model_folded_model": np.ones(50),
        }
    )

    plot_odd_even_transit(fold, tls)


@pytest.mark.unit
def test_aperture_if_flux_all_nan_swaps_flux_and_error():
    # Local ePSF extraction near the cutout edge (Kepler-1892 s14).
    from quicklook.tql import _aperture_if_flux_all_nan

    lc = _tglc_lc(psf_finite=False)
    lc.flux_err = np.full(len(lc), np.nan)
    lc["cal_aper_flux"] = 2.0 + 0.01 * np.sin(np.arange(len(lc)))
    lc.sector = 14

    out = _aperture_if_flux_all_nan(lc)

    assert np.allclose(out.flux.value, lc["cal_aper_flux"])
    assert np.isfinite(out.flux_err.value).all() and (out.flux_err.value > 0).all()
    assert out.sector == 14
    assert np.isnan(lc.flux.value).all()  # input untouched


@pytest.mark.unit
def test_aperture_if_flux_all_nan_keeps_finite_flux():
    from quicklook.tql import _aperture_if_flux_all_nan

    lc = _tglc_lc(psf_finite=True)

    assert _aperture_if_flux_all_nan(lc) is lc


@pytest.mark.unit
def test_tls_search_falls_back_when_extended_range_raises(monkeypatch):
    # Sparse data over the extended range empties TLS's duration cache (Kepler-1928).
    import quicklook.tql as tql

    calls = []

    class FakeTLS:
        def __init__(self, *a, **k):
            pass

        def power(self, **kw):
            calls.append(dict(kw))
            if kw.get("n_transits_min") == 1:
                raise ValueError("zero-size array to reduction operation minimum")
            return SimpleNamespace(period=1.0, distinct_transit_count=5)

    monkeypatch.setattr(tql, "tls", FakeTLS)
    ql = _bare_ql()
    ql.Porb_min, ql.Porb_max = 0.1, 26.0
    ql.tls_use_threads, ql.use_star_priors = None, False
    t = np.linspace(0, 26, 400)
    lc = lk.LightCurve(time=t, flux=np.ones_like(t), flux_err=np.full_like(t, 1e-3))

    result = ql._tls_search(lc)

    assert result.period == 1.0
    assert len(calls) == 2
    assert calls[1]["period_max"] == pytest.approx(13.0)
    assert "n_transits_min" not in calls[1]


@pytest.mark.unit
def test_tls_search_reraises_outside_extended_range(monkeypatch):
    import quicklook.tql as tql

    class FakeTLS:
        def __init__(self, *a, **k):
            pass

        def power(self, **kw):
            raise ValueError("boom")

    monkeypatch.setattr(tql, "tls", FakeTLS)
    ql = _bare_ql()
    ql.Porb_min, ql.Porb_max = 0.1, 5.0
    ql.tls_use_threads, ql.use_star_priors = None, False
    t = np.linspace(0, 26, 400)
    lc = lk.LightCurve(time=t, flux=np.ones_like(t), flux_err=np.full_like(t, 1e-3))

    with pytest.raises(ValueError, match="boom"):
        ql._tls_search(lc)


@pytest.mark.unit
def test_odd_even_plot_with_scalar_tls_model():
    # Degenerate TLS fit returns scalar NaN model light curves (HD 114082).
    import matplotlib

    matplotlib.use("Agg")
    from quicklook.plot import plot_odd_even_transit
    from quicklook.tql import _TLSResult

    t = np.linspace(0, 27, 500)
    fold = lk.LightCurve(time=t, flux=np.ones_like(t)).fold(period=10.0, epoch_time=5.0)
    tls = _TLSResult(
        {
            "duration": 0.1,
            "period": 10.0,
            "T0": 5.0,
            "depth": 0.99,
            "model_folded_phase": np.linspace(0, 1, 50),
            "model_folded_model": np.ones(50),
            "model_lightcurve_time": np.nan,
            "model_lightcurve_model": np.nan,
        }
    )

    plot_odd_even_transit(fold, tls)


@pytest.mark.unit
@pytest.mark.parametrize(
    "flux",
    [np.full(50, np.nan), np.ones(50), np.array([1.0, 1.01])],
    ids=["all-nan", "constant", "too-few"],
)
def test_require_usable_flux_rejects_degenerate_extraction(flux):
    # CCD-edge stars: PSF all NaN, aperture constant (Kepler-1892 s14).
    from quicklook.tql import _require_usable_flux

    lc = lk.LightCurve(time=np.arange(len(flux), dtype=float), flux=flux)

    with pytest.raises(ValueError, match="no usable photometry"):
        _require_usable_flux(lc)


@pytest.mark.unit
def test_require_usable_flux_accepts_varying_flux():
    from quicklook.tql import _require_usable_flux

    t = np.arange(50, dtype=float)
    _require_usable_flux(lk.LightCurve(time=t, flux=1 + 1e-3 * np.sin(t)))


@pytest.mark.unit
def test_odd_even_plot_when_tls_fit_no_transit():
    # TLS "No transit were fit": period set, duration and T0 NaN (HD 114082).
    import matplotlib

    matplotlib.use("Agg")
    from quicklook.plot import plot_odd_even_transit
    from quicklook.tql import _TLSResult

    t = np.linspace(0, 27, 500)
    fold = lk.LightCurve(time=t, flux=np.ones_like(t)).fold(period=10.0, epoch_time=5.0)
    tls = _TLSResult({"duration": np.nan, "period": 10.0, "T0": np.nan, "depth": np.nan})

    ax = plot_odd_even_transit(fold, tls)

    assert ax.texts and ax.texts[0].get_text() == "No TLS transit fit"


@pytest.mark.unit
def test_secondary_eclipse_plot_when_tls_fit_no_transit():
    import matplotlib

    matplotlib.use("Agg")
    from quicklook.plot import plot_secondary_eclipse
    from quicklook.tql import _TLSResult

    t = np.linspace(0, 27, 500)
    lc = lk.LightCurve(time=t, flux=np.ones_like(t))
    tls = _TLSResult({"duration": np.nan, "period": 10.0, "T0": np.nan, "depth": np.nan})

    ax = plot_secondary_eclipse(lc, tls, np.zeros_like(t, dtype=bool))

    assert ax.texts and ax.texts[0].get_text() == "No TLS transit fit"
