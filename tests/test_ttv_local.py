"""Local-baseline TTV diagnostics on a synthetic spotted star (offline)."""

import os

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as pl
import numpy as np
import pytest
from lightkurve import LightCurve

from quicklook.plot_ttv_local import plot_transit_gallery, plot_ttv_diagnostics
from quicklook.tql import TessQuickLook, _TLSResult
from quicklook.ttv import trapezoid
from quicklook.ttv_local import covered_epochs, measure_offset, run_local_ttv, split_visits

CADENCE = 2 / 1440
P, T0, DUR, DEPTH, NOISE = 6.0, 1.0, 2.5 / 24, 2.5e-3, 1e-3
VISITS = ((0.0, 20.0), (300.0, 320.0), (600.0, 620.0))
SPOT_AMP, P_ROT = 20e-3, 3.1


def _star(visit_offsets_min=(0, 0, 0), seed=3):
    rng = np.random.default_rng(seed)
    t = np.concatenate([np.arange(a, b, CADENCE) for a, b in VISITS])
    f = (
        1
        + SPOT_AMP * np.sin(2 * np.pi * t / P_ROT)
        + 0.4 * SPOT_AMP * np.sin(4 * np.pi * t / P_ROT + 1)
    )
    f += rng.normal(0, NOISE, t.size)
    for (a, b), off in zip(VISITS, visit_offsets_min):
        for n in range(int(np.ceil((a - T0) / P)), int((b - T0) / P) + 1):
            f -= DEPTH * trapezoid(t, T0 + n * P + off / 1440, DUR)
    return t, f


def test_split_visits_and_covered_epochs():
    t, _ = _star()
    visits = split_visits(t)
    assert len(visits) == 3
    assert covered_epochs(t[visits[0]], P, T0) == [0, 1, 2, 3]


def test_single_transit_offset_on_spotted_baseline():
    t, f = _star((45, 45, 45))
    tc = T0 + 2 * P
    m = measure_offset(t, f, [tc], DUR, NOISE)
    assert m.usable
    assert m.offset * 1440 == pytest.approx(45, abs=8)


def test_visit_offsets_recovered_and_flagged_nonlinear():
    offsets = (0, 70, -20)
    t, f = _star(offsets)
    r = run_local_ttv(t, f, P, T0, DUR, source="toi", n_injections=5)
    assert r.duration * 24 == pytest.approx(2.5, abs=0.6)
    got = [v.fit.offset * 1440 for v in r.usable_visits()]
    assert np.allclose(got, offsets, atol=8)
    assert r.visit_pvalue < 1e-3
    assert all(np.isfinite(v.err_calibrated) for v in r.visits)
    assert len(r.usable_singles()) >= 9


def test_linear_ephemeris_not_flagged():
    t, f = _star((0, 0, 0), seed=5)
    r = run_local_ttv(t, f, P, T0, DUR, n_injections=0)
    assert r.visit_pvalue > 1e-3
    assert abs(r.best_common_offset * 1440) < 10


def test_no_covered_transits_returns_none():
    t, f = _star()
    assert run_local_ttv(t[:100], f[:100], P, T0 + 0.5 * P, DUR, n_injections=0) is None


def test_figures_render():
    t, f = _star((0, 40, -20))
    r = run_local_ttv(t, f, P, T0, DUR, source="toi", n_injections=0)
    for fig in (plot_ttv_diagnostics(r, "synthetic"), plot_transit_gallery(r, t, f, "synthetic")):
        assert len(fig.axes) >= 3
        pl.close(fig)


# --- pipeline wiring ------------------------------------------------------------------


def _ql(enabled=True, sde=4.0):
    t, f = _star((0, 60, -30))
    ql = TessQuickLook.__new__(TessQuickLook)
    ql.verbose = False
    ql.ttv_search = enabled
    ql.single_event_search = False
    ql.ttv_injections = 0
    ql.target_name = "synthetic"
    ql.raw_lc = LightCurve(time=t, flux=f, flux_err=np.full_like(f, NOISE))
    ql.flat_lc = ql.raw_lc
    ql.exofop_data = {}
    ql.toi_epoch, ql.toi_period, ql.toi_dur = (T0, 0.001), (P, 1e-5), (DUR, 0.01)
    ql.tls_results = _TLSResult(
        {
            "SDE": sde,
            "period": 18.4,
            "T0": 5.0,
            "duration": 0.1,
            "distinct_transit_count": 2,
            "depth_variance_ratio": 1.0,
        }
    )
    return ql


def test_pipeline_falls_back_to_toi_ephemeris_when_tls_is_not_periodic(tmp_path):
    ql = _ql()
    ql.run_ttv_search()
    r = ql.run_ttv_diagnostics()
    assert r is not None and r.source == "toi"
    res = ql.tls_results
    assert res["ttv_local_source"] == "toi"
    assert np.allclose(res["ttv_local_visit_oc"] * 1440, [0, 60, -30], atol=8)
    assert res["ttv_local_visit_pvalue"] < 1e-3
    paths = ql._save_ttv_figures(tmp_path / "x_s01")
    assert [p.name for p in paths] == ["x_s01_ttv.png", "x_s01_ttv_transits.png"]
    assert all(p.stat().st_size > 0 for p in paths)


def test_pipeline_skips_without_ttv_search_or_with_masked_ephemeris(tmp_path):
    assert _ql(enabled=False).run_ttv_diagnostics() is None
    ql = _ql()
    ql.mask_ephem = True
    ql.run_ttv_search()
    assert ql.run_ttv_diagnostics() is None
    assert ql._save_ttv_figures(tmp_path / "x") == []


def test_ttv_local_fields_survive_hdf5(tmp_path):
    from quicklook import h5io

    ql = _ql()
    ql.run_ttv_search()
    ql.run_ttv_diagnostics()
    h5io.save(tmp_path / "x_tls.h5", dict(ql.tls_results))
    back = h5io.load(tmp_path / "x_tls.h5")
    assert np.allclose(back["ttv_local_visit_oc"], ql.tls_results["ttv_local_visit_oc"])
    assert str(back["ttv_local_source"]) == "toi"


def test_custom_ephem_tuple_from_cli_is_parsed_in_place():
    """Regression: typer passes --custom-ephem as a tuple; JD/hour conversion used to crash."""
    ql = TessQuickLook.__new__(TessQuickLook)
    ql.verbose = False
    ql.custom_ephem = (2458791.0538, 0.001, 5.690476, 1e-5, 1.75, 0.2)
    ql.parse_custom_ephem()
    assert ql.ephem_source == "custom"
    assert ql.toi_epoch[0] == pytest.approx(1791.0538)
    assert ql.toi_period[0] == pytest.approx(5.690476)
    assert ql.toi_dur[0] == pytest.approx(1.75 / 24)
