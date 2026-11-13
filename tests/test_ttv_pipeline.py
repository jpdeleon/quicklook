"""TTV search wired into TessQuickLook, the iterative search, figure, CLI, GUI and summary.

Offline: the pipeline object is built with ``__new__`` and synthetic light curves.
"""

import os

os.environ.setdefault("MPLBACKEND", "Agg")

from pathlib import Path
from unittest.mock import patch

import matplotlib.pyplot as pl
import numpy as np
import pytest
from lightkurve import LightCurve

from quicklook.tql import TessQuickLook, _TLSResult
from quicklook.ttv import trapezoid

REPO = Path(__file__).resolve().parent.parent
CADENCE = 10 / 1440
PERIOD, T0, DUR, DEPTH, NOISE = 3.0, 1.0, 0.12, 3e-3, 5e-4
N_TRANSITS = 15


def _lc(shifts, period=PERIOD, t0=T0, dur=DUR, depth=DEPTH, span=None, seed=1):
    rng = np.random.default_rng(seed)
    t = np.arange(0.0, span or t0 + len(shifts) * period, CADENCE)
    f = 1 + rng.normal(0, NOISE, t.size)
    for k, s in enumerate(shifts):
        f -= depth * trapezoid(t, t0 + k * period + s, dur)
    return LightCurve(time=t, flux=f, flux_err=np.full_like(f, NOISE))


def _ql(shifts=None, sde=20.0, enabled=True):
    shifts = np.zeros(N_TRANSITS) if shifts is None else shifts
    ql = TessQuickLook.__new__(TessQuickLook)
    ql.verbose = False
    ql.ttv_search = enabled
    ql.single_event_search = False
    ql.flat_lc = ql.raw_lc = _lc(shifts)
    ql.exofop_data = {}
    ql.tls_results = _TLSResult(
        {
            "SDE": sde,
            "period": PERIOD,
            "T0": T0,
            "duration": DUR,
            "distinct_transit_count": N_TRANSITS,
            "depth_variance_ratio": 1.0,
        }
    )
    return ql


SINUSOID = 0.02 * np.sin(2 * np.pi * np.arange(N_TRANSITS) / 7)
ALTERNATING = 0.08 * np.where(np.arange(N_TRANSITS) % 2, 1, -1)  # > half a duration


# --- run ------------------------------------------------------------------------------


def test_disabled_search_adds_nothing():
    ql = _ql(enabled=False)
    assert ql.run_ttv_search() is None
    assert "ttv_search" not in ql.tls_results


def test_non_periodic_tls_signal_is_not_timed():
    ql = _ql(sde=4.0)
    assert ql.run_ttv_search() is None
    r = ql.tls_results
    assert r["ttv_search"] is True
    assert r["n_ttv_transits"] == 0
    assert r["ttv_detected"] is False
    assert np.isnan(r["ttv_pvalue"])


def test_periodic_signal_with_ttvs_is_detected():
    ql = _ql(SINUSOID)
    stats = ql.run_ttv_search()
    r = ql.tls_results
    assert stats.detected and r["ttv_detected"] is True
    assert r["ttv_source"] == "tls"
    assert r["n_ttv_transits"] == N_TRANSITS == r["ttv_tcs"].size == r["ttv_oc"].size
    assert r["ttv_amplitude"] == pytest.approx(0.02, rel=0.3)
    assert r["ttv_pvalue"] < 1e-6
    assert r["ttv_period"] == pytest.approx(PERIOD, abs=3e-3)  # the sinusoid tilts the line
    assert r["ttv_n_cycle_solutions"] == 1
    assert r["ttv_alt_periods"].size == 0


def test_periodic_signal_without_ttvs_is_not_flagged():
    ql = _ql()
    ql.run_ttv_search()
    assert ql.tls_results["ttv_detected"] is False
    assert ql.tls_results["n_ttv_transits"] == N_TRANSITS


# --- iterative search masks at fitted times ---------------------------------------


def _residual_depth(lc, tcs):
    t, f = lc.time.value, lc.flux.value
    near = np.zeros(t.size, dtype=bool)
    for tc in tcs:
        near |= np.abs(t - tc) < DUR / 2
    return 1 - np.mean(f[near]) if near.any() else 0.0


def test_primary_mask_follows_fitted_times():
    ql = _ql(ALTERNATING)
    ql.run_ttv_search()
    tcs = T0 + PERIOD * np.arange(N_TRANSITS) + ALTERNATING
    linear = ql._transit_mask_for_lc(ql.flat_lc, ql.tls_results)
    fitted = ql.primary_transit_mask(ql.flat_lc)
    # The linear-ephemeris mask leaves part of every displaced transit behind.
    assert _residual_depth(ql.flat_lc[~linear], tcs) > DEPTH / 4
    assert _residual_depth(ql.flat_lc[~fitted], tcs) == 0.0


def test_primary_mask_is_linear_without_ttv_search():
    ql = _ql(ALTERNATING, enabled=False)
    ql.run_ttv_search()
    linear = ql._transit_mask_for_lc(ql.flat_lc, ql.tls_results)
    assert np.array_equal(ql.primary_transit_mask(ql.flat_lc), linear)


def test_iterative_search_runs_on_the_fitted_time_mask():
    ql = _ql(ALTERNATING)
    ql.min_sde_iterative, ql.max_planets = 5.0, None
    ql.run_ttv_search()
    seen = []

    def fake_search(lc, *args, **kwargs):
        seen.append(lc)
        return _TLSResult({"SDE": 1.0})

    with patch.object(ql, "_tls_search", side_effect=fake_search):
        ql.run_iterative_tls()
    tcs = T0 + PERIOD * np.arange(N_TRANSITS) + ALTERNATING
    assert _residual_depth(seen[0], tcs) == 0.0


# --- chaining single events -----------------------------------------------------------


def test_three_chained_single_events_make_a_periodic_ttv_target():
    shifts = [0.0, 0.06, -0.05, 0.03]
    lc = _lc(shifts, period=20.0, t0=5.0, dur=0.2, depth=4e-3, span=70.0)
    ql = _ql(sde=4.0)
    ql.single_event_search = True
    ql.flat_lc = ql.raw_lc = lc
    ql.window_length = 0.5
    ql.run_ttv_search()
    ql.run_single_event_search()
    assert len(ql.unexplained_events) >= 3
    stats = ql.run_ttv_from_chain()
    r = ql.tls_results
    assert stats is not None and stats.n >= 3
    assert r["ttv_source"] == "chain"
    assert r["event_class"] == "periodic" == ql.event_class
    assert r["ttv_period"] == pytest.approx(20.0, abs=0.05)


def test_chain_is_skipped_when_tls_was_already_timed():
    ql = _ql(SINUSOID)
    ql.single_event_search = True
    ql.single_event_result = object()
    ql.run_ttv_search()
    assert ql.run_ttv_from_chain() is None
    assert ql.tls_results["ttv_source"] == "tls"


# --- figure and summary ---------------------------------------------------------------


def test_ttv_layout_and_panels_render():
    ql = _ql(SINUSOID)
    ql.run_ttv_search()
    assert ql.has_ttv_layout
    fig, axes = pl.subplots(3, 3)
    ql._plot_ttv_panels(axes.flatten())
    assert "aligned" in axes.flatten()[6].get_title()
    assert axes.flatten()[7].containers  # O-C error bars
    assert "TTVs detected" in axes.flatten()[7].get_title()
    pl.close(fig)


def test_no_ttv_layout_for_mono_targets_or_without_search():
    ql = _ql(SINUSOID)
    ql.run_ttv_search()
    ql.event_class = "mono"
    assert not ql.has_ttv_layout
    assert not _ql(SINUSOID, enabled=False).has_ttv_layout


def test_ttv_summary_rows():
    ql = _ql(SINUSOID)
    ql.run_ttv_search()
    rows = dict(ql._ttv_candidate_rows())
    assert rows["TTVs"] == "detected"
    assert rows["Transits timed"] == str(N_TRANSITS)
    assert rows["TTV amplitude"].endswith("min")


# --- CLI / GUI / summary-table wiring -------------------------------------------------


@pytest.mark.parametrize("flag, expected", [([], False), (["--ttv-search"], True)])
def test_cli_passes_ttv_search(flag, expected):
    from typer.testing import CliRunner

    from quicklook.cli.app import app

    with patch("quicklook.tql.TessQuickLook") as ql:
        CliRunner().invoke(app, ["run", "--name", "TOI-1234", "--sector", "5", "--save"] + flag)
    assert ql.call_args.kwargs["ttv_search"] is expected


def test_cli_each_sector_helper_accepts_ttv_search():
    import inspect

    from quicklook.cli.app import _run_ql_for_sector

    assert "ttv_search" in inspect.signature(_run_ql_for_sector).parameters


def test_gui_ttv_search_wired():
    html = (REPO / "quicklook/app/templates/index.html").read_text()
    app_py = (REPO / "quicklook/app/app.py").read_text()
    assert 'name="ttv_search"' in html
    assert "'ttv_search'" in html  # persisted between visits
    assert '"ttv_search": _is_truthy(form.get("ttv_search"))' in app_py
    assert 'ttv_search=kwargs.get("ttv_search"' in app_py


def test_tls_summary_has_ttv_columns():
    from quicklook.app.tls_summary import extract_fields

    data = {"ttv_detected": True, "ttv_amplitude": 0.01, "ttv_pvalue": 1e-8}
    fields = extract_fields(data, "x_s01_pdcsap_sc_tls")
    assert fields["ttv_detected"] is True
    assert fields["ttv_amp_min"] == pytest.approx(14.4)
    assert fields["ttv_pvalue"] == pytest.approx(1e-8)
    html = (REPO / "quicklook/app/templates/tls_summary.html").read_text()
    assert 'data-key="ttv_amp_min"' in html and 'data-key="ttv_pvalue"' in html
    assert html.count("<th data-key=") == html.count('<th><input type="text" data-col=')


def test_ttv_fields_survive_hdf5_round_trip(tmp_path):
    from quicklook import h5io

    ql = _ql(SINUSOID)
    ql.run_ttv_search()
    path = tmp_path / "x_tls.h5"
    h5io.save(path, dict(ql.tls_results))
    back = h5io.load(path)
    assert bool(back["ttv_detected"]) is True
    assert np.allclose(back["ttv_tcs"], ql.tls_results["ttv_tcs"])
