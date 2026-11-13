"""Single-event search wired into TessQuickLook, the figure, CLI, GUI and summary table.

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

CADENCE = 120 / 86400
REPO = Path(__file__).resolve().parent.parent


def _flat_lc(events=((14.0, 0.25, 0.005),), seed=0):
    rng = np.random.default_rng(seed)
    t = np.arange(0.0, 27.0, CADENCE)
    t = t[(t < 12.5) | (t > 13.5)]
    f = 1 + rng.normal(0, 5e-4, t.size)
    for t0, dur, depth in events:
        f[np.abs(t - t0) < dur / 2] -= depth
    return LightCurve(time=t, flux=f, flux_err=np.full_like(f, 5e-4))


def _ql(events=((14.0, 0.25, 0.005),), tls=None, enabled=True):
    ql = TessQuickLook.__new__(TessQuickLook)
    ql.verbose = False
    ql.single_event_search = enabled
    ql.flat_lc = _flat_lc(events)
    ql.exofop_data = {}
    ql.tls_results = _TLSResult(
        tls or {"SDE": 4.0, "period": 3.3, "T0": 1.0, "duration": 0.1, "distinct_transit_count": 7}
    )
    return ql


# --- run + classification -------------------------------------------------------


def test_disabled_search_adds_nothing():
    ql = _ql(enabled=False)
    assert ql.run_single_event_search() is None
    assert "event_class" not in ql.tls_results


def test_single_dip_without_tls_signal_is_mono():
    ql = _ql()
    ql.run_single_event_search()
    r = ql.tls_results
    assert r["event_class"] == "mono"
    assert r["n_events"] == 1
    assert r["single_event_t0"] == pytest.approx(14.0, abs=0.05)
    assert r["single_event_depth"] == pytest.approx(0.005, rel=0.25)
    assert r["p_min"] >= 27.0 - 14.0 - 0.3
    assert r["single_event_search"] is True
    assert isinstance(r["event_t0s"], np.ndarray) and r["event_t0s"].size == 1
    assert ql.event_class == "mono"


def test_no_dip_and_no_tls_signal_is_none():
    ql = _ql(events=())
    ql.run_single_event_search()
    assert ql.tls_results["event_class"] == "none"
    assert ql.tls_results["n_events"] == 0
    assert np.isnan(ql.tls_results["p_min"])


def test_tls_periodic_signal_explaining_the_events_is_periodic():
    events = tuple((2.0 + 4.0 * k, 0.15, 0.006) for k in range(6) if not 12 < 2 + 4 * k < 14)
    tls = {"SDE": 25.0, "period": 4.0, "T0": 2.0, "duration": 0.15, "distinct_transit_count": 6}
    ql = _ql(events=events, tls=tls)
    ql.run_single_event_search()
    assert ql.tls_results["event_class"] == "periodic"
    assert ql.tls_results["n_unexplained_events"] == 0


def test_extra_event_outside_tls_transits_is_kept_as_unexplained():
    events = tuple((2.0 + 4.0 * k, 0.15, 0.006) for k in range(6) if not 12 < 2 + 4 * k < 14)
    events = events + ((15.3, 0.3, 0.008),)
    tls = {"SDE": 25.0, "period": 4.0, "T0": 2.0, "duration": 0.15, "distinct_transit_count": 6}
    ql = _ql(events=events, tls=tls)
    ql.run_single_event_search()
    assert ql.tls_results["event_class"] == "periodic"
    assert ql.tls_results["n_unexplained_events"] == 1
    assert ql.tls_results["single_event_t0"] == pytest.approx(15.3, abs=0.05)


def test_tls_with_one_transit_is_not_periodic():
    tls = {"SDE": 30.0, "period": 20.0, "T0": 14.0, "duration": 0.25, "distinct_transit_count": 1}
    ql = _ql(tls=tls)
    ql.run_single_event_search()
    assert ql.tls_results["event_class"] == "mono"


def test_stellar_density_used_when_known():
    ql = _ql()
    ql.exofop_data = {"stellar_parameters": [{"srad": "2.0", "mass": "1.0", "sdate": "2020-01-01"}]}
    assert ql._stellar_density() == pytest.approx(1.0 / 8.0)


def test_stellar_density_defaults_to_solar():
    assert _ql()._stellar_density() == 1.0


# --- summary + figure -----------------------------------------------------------


def test_mono_summary_rows_describe_the_event():
    ql = _ql()
    ql.run_single_event_search()
    labels = [label for label, _ in ql._single_event_candidate_rows(r_star=1.0)]
    for expected in ("Class", "T0 (event)", "Depth (event)", "Dur. (event)", "SNR (event)"):
        assert expected in labels
    assert any(label.startswith("P min") for label in labels)
    assert any(label.startswith(r"$\Delta$BIC") for label in labels)


def test_mono_panels_render():
    ql = _ql()
    ql.run_single_event_search()
    fig, axes = pl.subplots(3, 3)
    ql._plot_single_event_panels(axes.flatten())
    assert axes.flatten()[4].get_xlabel().startswith("Time")
    assert "Period" in axes.flatten()[7].get_xlabel()
    assert len(axes.flatten()[6].lines) >= 2  # transit + systematics models
    pl.close(fig)


def test_extra_event_figure_written_for_periodic_targets(tmp_path):
    events = tuple((2.0 + 4.0 * k, 0.15, 0.006) for k in range(6) if not 12 < 2 + 4 * k < 14)
    events = events + ((15.3, 0.3, 0.008),)
    tls = {"SDE": 25.0, "period": 4.0, "T0": 2.0, "duration": 0.15, "distinct_transit_count": 6}
    ql = _ql(events=events, tls=tls)
    ql.run_single_event_search()
    out = ql._save_single_event_figure(tmp_path / "target_s01_pdcsap_sc")
    assert out == tmp_path / "target_s01_pdcsap_sc_mono.png"
    assert out.exists()


def test_extra_event_figure_skipped_without_unexplained_events(tmp_path):
    ql = _ql(events=())
    ql.run_single_event_search()
    assert ql._save_single_event_figure(tmp_path / "x") is None


# --- CLI / GUI / summary-table wiring ---------------------------------------------


@pytest.mark.parametrize("flag, expected", [([], False), (["--single-event-search"], True)])
def test_cli_passes_single_event_search(flag, expected):
    from typer.testing import CliRunner

    from quicklook.cli.app import app

    with patch("quicklook.tql.TessQuickLook") as ql:
        CliRunner().invoke(app, ["run", "--name", "TOI-1234", "--sector", "5", "--save"] + flag)
    assert ql.call_args.kwargs["single_event_search"] is expected


def test_cli_each_sector_helper_accepts_single_event_search():
    import inspect

    from quicklook.cli.app import _run_ql_for_sector

    assert "single_event_search" in inspect.signature(_run_ql_for_sector).parameters


def test_gui_single_event_search_wired():
    html = (REPO / "quicklook/app/templates/index.html").read_text()
    app_py = (REPO / "quicklook/app/app.py").read_text()
    assert 'name="single_event_search"' in html
    assert "'single_event_search'" in html  # persisted between visits
    assert '"single_event_search": _is_truthy(form.get("single_event_search"))' in app_py
    assert 'single_event_search=kwargs.get("single_event_search"' in app_py


def test_tls_summary_extracts_event_class_and_count():
    from quicklook.app.tls_summary import extract_fields

    fields = extract_fields(
        {"event_class": "mono", "n_events": 1, "p_min": 13.2}, "x_s01_pdcsap_sc_tls"
    )
    assert fields["event_class"] == "mono"
    assert fields["n_events"] == 1
    assert fields["p_min"] == pytest.approx(13.2)


def test_tls_summary_table_has_event_class_column():
    html = (REPO / "quicklook/app/templates/tls_summary.html").read_text()
    assert 'data-key="event_class"' in html
    assert 'data-key="n_events"' in html
    n_headers = html.count("<th data-key=")
    n_filters = html.count('<th><input type="text" data-col=')
    assert n_headers == n_filters


def test_gui_progress_tracks_single_event_plotting():
    import re

    app_py = (REPO / "quicklook/app/app.py").read_text()
    pattern = re.search(r'\(r"([^"]*Plotting single-event[^"]*)", "Generating plots"\)', app_py)
    assert pattern and re.search(pattern.group(1), "Plotting single-event panels...")
