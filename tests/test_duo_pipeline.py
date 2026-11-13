"""Duotransit classification wired into TessQuickLook, the figure and the summary table."""

import os

os.environ.setdefault("MPLBACKEND", "Agg")

from pathlib import Path

import matplotlib.pyplot as pl
import numpy as np
import pytest
from lightkurve import LightCurve

from quicklook.tql import TessQuickLook, _TLSResult

CADENCE = 120 / 86400
REPO = Path(__file__).resolve().parent.parent


def _two_sector_lc(events=(10.0, 710.0), depths=(0.004, 0.004), seed=3):
    rng = np.random.default_rng(seed)
    t = np.concatenate([np.arange(0.0, 27.0, CADENCE), np.arange(700.0, 727.0, CADENCE)])
    f = 1 + rng.normal(0, 5e-4, t.size)
    for t0, depth in zip(events, depths):
        f[np.abs(t - t0) < 0.15] -= depth
    return LightCurve(time=t, flux=f, flux_err=np.full_like(f, 5e-4))


def _ql(lc, tls=None):
    ql = TessQuickLook.__new__(TessQuickLook)
    ql.verbose = False
    ql.single_event_search = True
    ql.flat_lc = ql.raw_lc = lc
    ql.window_length = 0.5
    ql.exofop_data = {}
    ql.tls_results = _TLSResult(
        tls or {"SDE": 4.0, "period": 3.3, "T0": 1.0, "duration": 0.1, "distinct_transit_count": 9}
    )
    return ql


def test_two_consistent_events_are_a_duo():
    ql = _ql(_two_sector_lc())
    ql.run_single_event_search()
    r = ql.tls_results
    assert r["event_class"] == "duo"
    assert r["duo_delta_t"] == pytest.approx(700.0, abs=0.1)
    assert np.allclose(r["duo_t0s"], [10.0, 710.0], atol=0.05)
    assert r["n_aliases"] == len(r["p_aliases"]) == len(r["p_alias_weights"]) > 1
    assert r["p_alias_weights"].sum() == pytest.approx(1.0)
    assert r["p_min"] == pytest.approx(r["p_aliases"].min())
    assert r["p_best"] in r["p_aliases"]
    assert ql.duo is not None


def test_tls_with_two_transits_does_not_block_a_duo():
    # TLS can fit P=dT with exactly 2 transits; that is one alias, not a period.
    tls = {"SDE": 20.0, "period": 700.0, "T0": 10.0, "duration": 0.3, "distinct_transit_count": 2}
    ql = _ql(_two_sector_lc(), tls=tls)
    ql.run_single_event_search()
    assert ql.tls_results["event_class"] == "duo"


def test_inconsistent_depths_stay_mono():
    ql = _ql(_two_sector_lc(depths=(0.002, 0.010)))
    ql.run_single_event_search()
    assert ql.tls_results["event_class"] == "mono"
    assert ql.tls_results["n_aliases"] == 0


def test_duo_summary_rows():
    ql = _ql(_two_sector_lc())
    ql.run_single_event_search()
    labels = [label for label, _ in ql._single_event_candidate_rows(r_star=1.0)]
    for expected in ("Class", "T0 (1st)", "T0 (2nd)", r"$\Delta$T", "Aliases", "P (most likely)"):
        assert expected in labels


def test_duo_panels_render():
    ql = _ql(_two_sector_lc())
    ql.run_single_event_search()
    fig, axes = pl.subplots(3, 3)
    ql._plot_single_event_panels(axes.flatten())
    assert axes.flatten()[6].get_legend() is not None
    assert len(axes.flatten()[6].lines) >= 4  # data + model for each event
    assert axes.flatten()[7].get_xlabel().startswith("Period")
    assert len(axes.flatten()[7].patches) == ql.tls_results["n_aliases"]
    pl.close(fig)


def test_tls_summary_has_alias_count_column():
    from quicklook.app.tls_summary import extract_fields

    fields = extract_fields({"event_class": "duo", "n_aliases": 12}, "x_s01_pdcsap_sc_tls")
    assert fields["n_aliases"] == 12
    html = (REPO / "quicklook/app/templates/tls_summary.html").read_text()
    assert 'data-key="n_aliases"' in html
    assert html.count("<th data-key=") == html.count('<th><input type="text" data-col=')


def test_stellar_density_skips_entries_without_mass():
    # Regression (TOI-2180): the newest ExoFOP entry has a radius but no mass,
    # which silently fell back to solar density and mis-ranked duo aliases.
    ql = _ql(_two_sector_lc())
    ql.exofop_data = {
        "stellar_parameters": [
            {"srad": None, "mass": None},
            {"srad": "1.62762", "mass": "1.025", "sdate": "2019-04-15"},
            {"srad": "1.6276", "mass": "", "sdate": "2026-05-05 12:03:40"},
        ]
    }
    assert ql._stellar_density() == pytest.approx(1.025 / 1.62762**3)


def test_duo_time_panels_are_broken_between_sectors():
    ql = _ql(_two_sector_lc())
    ql.run_single_event_search()
    fig, axes = pl.subplots(3, 3)
    flat = axes.flatten()
    ql._plot_single_event_panels(flat)
    # Panels 4 and 5 (flattened flux, single-event SNR) each become 2 sub-axes.
    for idx in (3, 4):
        subs = ql._broken_panels[idx]
        assert len(subs) == 2
        assert subs[0].get_xlim()[1] < 100 < subs[1].get_xlim()[0]
    assert flat[3] not in fig.axes and flat[4] not in fig.axes
    pl.close(fig)


def test_single_sector_mono_time_panels_are_not_broken():
    from test_single_event_pipeline import _ql as mono_ql

    ql = mono_ql()
    ql.run_single_event_search()
    fig, axes = pl.subplots(3, 3)
    ql._plot_single_event_panels(axes.flatten())
    assert all(len(v) == 1 for v in ql._broken_panels.values())
    pl.close(fig)
