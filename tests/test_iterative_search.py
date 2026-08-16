"""Offline tests for the iterative multi-planet TLS search.

Covers the masking loop in ``run_iterative_tls`` (stop conditions: SDE
threshold, duplicate period, max planets), the per-planet PNG/HDF5 outputs,
the h5 summary metadata, and the CLI/GUI wiring of the opt-in flags. All
tests bypass ``TessQuickLook.__init__`` so nothing touches the network.
"""

import os

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib

matplotlib.use("Agg")

from pathlib import Path

import numpy as np
import pytest
from unittest.mock import patch

from quicklook import tql as tql_mod
from quicklook.tql import TessQuickLook, _TLSResult


class _ScriptedTLS:
    """Stand-in for ``transitleastsquares.transitleastsquares``.

    Each ``.power()`` call returns the next queued result; an unexpected
    extra search raises so tests fail loudly if the loop mis-terminates.
    """

    results = []

    def __init__(self, *args, **kwargs):
        pass

    def power(self, **kwargs):
        if not _ScriptedTLS.results:
            raise AssertionError("Unexpected extra TLS search")
        return _ScriptedTLS.results.pop(0)


def _candidate(period, T0, sde, dur=0.1, **extra):
    """Build a TLS-compatible results dict for the primary or a candidate."""
    values = {
        "period": period,
        "period_uncertainty": 1e-4,
        "T0": T0,
        "duration": dur,
        "duration_expected": dur,
        "depth": 0.98,
        "rp_rs": 0.05,
        "SDE": sde,
        "FAP": 1e-6,
        "odd_even_mismatch": 0.5,
        "distinct_transit_count": 5,
        "transit_depths": None,
    }
    values.update(extra)
    return _TLSResult(values)


def _build_iterative_ql(monkeypatch, primary, flat_lc):
    """Build a TessQuickLook with just enough state for the iterative search."""
    monkeypatch.setattr(tql_mod, "tls", _ScriptedTLS)
    monkeypatch.setattr(tql_mod, "_get_gpu_tls", lambda: None)

    ql = TessQuickLook.__new__(TessQuickLook)
    ql.flat_lc = flat_lc
    ql.tls_results = primary
    ql.verbose = False
    ql.tls_use_threads = None
    ql.use_star_priors = False
    ql.Porb_min = 0.1
    ql.Porb_max = 13.0
    ql.min_sde_iterative = 5.0
    ql.max_planets = None
    ql.iterative_search = True
    return ql


# --- masking helpers --------------------------------------------------------


def test_period_is_duplicate():
    primary = _candidate(3.0, 5.0, 20.0)
    other = _candidate(2.0, 5.0, 12.0)
    assert TessQuickLook._period_is_duplicate(3.0005, [primary, other]) is True
    assert TessQuickLook._period_is_duplicate(1.999, [primary, other]) is True
    assert TessQuickLook._period_is_duplicate(2.5, [primary, other]) is False


def test_transit_mask_for_lc(mock_light_curve):
    results = _candidate(3.0, 5.0, 20.0)
    mask = TessQuickLook._transit_mask_for_lc(mock_light_curve, results)
    assert mask.dtype == bool
    assert len(mask) == len(mock_light_curve.time)
    assert mask.sum() > 0
    # The cadence nearest the transit epoch must be masked.
    i_nearest = np.abs(mock_light_curve.time.value - 5.0).argmin()
    assert mask[i_nearest]


# --- run_iterative_tls loop -------------------------------------------------


def test_run_iterative_tls_finds_additional_planets(monkeypatch, mock_light_curve):
    primary = _candidate(3.0, 5.0, 21.0)
    _ScriptedTLS.results = [
        _candidate(2.0, 5.5, 12.0, transit_depths=np.array([0.02, 0.018, 0.021])),
        _candidate(1.5, 5.8, 4.0),  # below min_sde -> loop must stop
    ]
    ql = _build_iterative_ql(monkeypatch, primary, mock_light_curve)

    results = ql.run_iterative_tls()

    assert len(results) == 2
    assert results[0] is primary
    assert results[1].period == pytest.approx(2.0)
    assert results[1].SDE == pytest.approx(12.0)
    assert ql.iterative_tls_results == results
    # Vetting metrics were attached to the candidate.
    assert np.isfinite(results[1]["depth_variance_ratio"])
    assert results[1]["duration_consistency_ratio"] == pytest.approx(1.0)
    assert np.isfinite(results[1]["secondary_sde"])


def test_run_iterative_tls_stops_on_duplicate_period(monkeypatch, mock_light_curve):
    primary = _candidate(3.0, 5.0, 21.0)
    _ScriptedTLS.results = [_candidate(3.001, 5.5, 15.0)]  # same period as primary
    ql = _build_iterative_ql(monkeypatch, primary, mock_light_curve)

    results = ql.run_iterative_tls()

    assert len(results) == 1
    assert results[0] is primary


def test_run_iterative_tls_respects_max_planets(monkeypatch, mock_light_curve):
    primary = _candidate(3.0, 5.0, 21.0)
    _ScriptedTLS.results = [
        _candidate(2.0, 5.5, 12.0),
        _candidate(1.5, 5.8, 11.0),  # must never be reached
    ]
    ql = _build_iterative_ql(monkeypatch, primary, mock_light_curve)
    ql.max_planets = 2

    results = ql.run_iterative_tls()

    assert len(results) == 2
    assert results[1].period == pytest.approx(2.0)
    # The scripted TLS would raise on a third call; max_planets broke the loop.


def test_run_iterative_tls_skips_when_primary_below_threshold(monkeypatch, mock_light_curve):
    primary = _candidate(3.0, 5.0, 4.0)
    _ScriptedTLS.results = []  # any search would raise
    ql = _build_iterative_ql(monkeypatch, primary, mock_light_curve)

    results = ql.run_iterative_tls()

    assert len(results) == 1
    assert results[0] is primary
    assert _ScriptedTLS.results == []


def test_run_iterative_tls_without_primary_returns_empty():
    ql = TessQuickLook.__new__(TessQuickLook)
    ql.min_sde_iterative = 5.0
    ql.max_planets = None
    ql.verbose = False

    results = ql.run_iterative_tls()

    assert results == []
    assert ql.iterative_tls_results == []


# --- h5 metadata + per-planet outputs ---------------------------------------


def test_append_tls_results_records_iterative_metadata(monkeypatch):
    primary = _candidate(3.0, 5.0, 21.0)
    cand = _candidate(2.0, 5.5, 12.0)

    ql = TessQuickLook.__new__(TessQuickLook)
    ql.pipeline = "spoc"
    ql.flux_type = "pdcsap"
    ql.exptime = 120
    ql.cadence = "short"
    ql.Porb_min = 0.5
    ql.Porb_max = 10.0
    ql.toi_period = None
    ql.toi_epoch = None
    ql.toi_dur = None
    ql.toi_depth = None
    ql.gaiaid = 1
    ql.ticid = 2
    ql.toiid = None
    ql.sector = 1
    ql.flatten_method = "biweight"
    ql.window_length = 0.5
    ql.gls = None
    ql.simbad_obj_type = None
    ql.iterative_search = True
    ql.iterative_tls_results = [primary, cand]
    ql.tls_results = primary
    ql.use_star_priors = False
    ql.exofop_data = {"stub": True}

    class _Col:
        def __init__(self, a):
            self.value = a

    class _LC:
        time = _Col(np.array([0.0]))
        flux = _Col(np.array([1.0]))
        flux_err = _Col(np.array([0.001]))

    ql.raw_lc = _LC()
    ql.flat_lc = _LC()

    monkeypatch.setattr(TessQuickLook, "get_toi_radius", lambda self: None)

    ql.append_tls_results()

    assert ql.tls_results["iterative_search"] is True
    assert ql.tls_results["n_planets"] == 2
    rows = ql.tls_results["iterative_candidates"]
    assert len(rows) == 2
    assert [r["planet"] for r in rows] == [1, 2]
    assert rows[1]["period"] == pytest.approx(2.0)
    assert rows[1]["SDE"] == pytest.approx(12.0)


def test_save_iterative_candidates_writes_per_planet_files(tmp_path, monkeypatch, mock_light_curve):
    primary = _candidate(3.0, 5.0, 21.0)
    cand2 = _candidate(2.0, 5.5, 12.0)
    cand3 = _candidate(1.5, 5.8, 9.0)

    def fake_plot_tls(tls_results, toi_period=None, period_min=0.1, period_max=None, ax=None):
        ax.plot([0, 1], [0, 1])

    def fake_odd_even(fold_lc, tls_results, bin_mins=20, markersize=3, ax=None, phase_xlim=None):
        ax.scatter([0, 1], [0, 1])

    monkeypatch.setattr(tql_mod, "plot_tls", fake_plot_tls)
    monkeypatch.setattr(tql_mod, "plot_odd_even_transit", fake_odd_even)

    ql = _build_iterative_ql(monkeypatch, primary, mock_light_curve)
    ql.iterative_tls_results = [primary, cand2, cand3]
    ql.outdir = str(tmp_path)
    ql.savefig = True
    ql.savetls = True
    ql.phase_xlim = None

    fp = Path(tmp_path, "TIC1234_s01_pdcsap_1200c")
    ql._save_iterative_candidates(fp)

    base = "TIC1234_s01_pdcsap_1200c"
    assert (tmp_path / f"{base}_p2.png").exists()
    assert (tmp_path / f"{base}_p3.png").exists()
    assert (tmp_path / f"{base}_tls_p2.h5").exists()
    assert (tmp_path / f"{base}_tls_p3.h5").exists()
    # The primary's own outputs are written by plot_tql, not this helper.
    assert not (tmp_path / f"{base}.png").exists()
    assert not (tmp_path / f"{base}_tls.h5").exists()


# --- CLI + GUI wiring -------------------------------------------------------


def test_cli_iterative_flags_wired(monkeypatch):
    from typer.testing import CliRunner

    from quicklook.cli.app import app

    calls = {}

    class FakeQuickLook:
        def __init__(self, **kwargs):
            calls.update(kwargs)

        def plot_tql(self):
            return None

    monkeypatch.setattr(tql_mod, "TessQuickLook", FakeQuickLook)

    result = CliRunner().invoke(
        app,
        [
            "run",
            "--name",
            "TOI-1234",
            "--iterative",
            "--min-sde-iterative",
            "6",
            "--max-planets",
            "3",
            "--save",
        ],
    )

    assert result.exit_code == 0, result.output
    assert calls["iterative_search"] is True
    assert calls["min_sde_iterative"] == 6.0
    assert calls["max_planets"] == 3


def test_cli_iterative_defaults_off_and_sde5():
    from typer.testing import CliRunner

    from quicklook.cli.app import app

    calls = {}

    class FakeQuickLook:
        def __init__(self, **kwargs):
            calls.update(kwargs)

        def plot_tql(self):
            return None

    with patch("quicklook.tql.TessQuickLook", FakeQuickLook):
        result = CliRunner().invoke(app, ["run", "--name", "TOI-1234", "--save"])

    assert result.exit_code == 0, result.output
    assert calls["iterative_search"] is False
    assert calls["min_sde_iterative"] == 5.0
    assert calls["max_planets"] is None


def test_gui_iterative_search_wired():
    """GUI checkbox + form parser + Flask worker must all agree."""
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    html = open(os.path.join(here, "quicklook/app/templates/index.html")).read()
    app_py = open(os.path.join(here, "quicklook/app/app.py")).read()

    assert 'name="iterative_search"' in html
    assert "'iterative_search'" in html  # in PERSIST_CHECKS
    assert '"iterative_search": _is_truthy(form.get("iterative_search"))' in app_py
    assert 'iterative_search=kwargs.get("iterative_search"' in app_py
