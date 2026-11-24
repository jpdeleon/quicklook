"""Precession search wired into TessQuickLook, the saved figure, the h5 summary and the CLI.

Offline: the pipeline object is built with ``__new__`` and a synthetic stitched light curve.
"""

import os

os.environ.setdefault("MPLBACKEND", "Agg")

from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest
import typer
from lightkurve import LightCurve
from typer.testing import CliRunner

from quicklook import h5io
from quicklook.cli.app import app
from quicklook.precession import SectorSpan, sector_span_from_lc
from quicklook.precession_pipeline import PRECESSION_FIGURE_SUFFIX, PRECESSION_FOLDS_SUFFIX
from quicklook.tql import TessQuickLook, _TLSResult
from test_precession import NOISE, PERIOD, T0, _duration, _stitched


def _ql(db_dt=0.04, enabled=True, sde=20.0, spans=None):
    time, flux, err, stitched_spans = _stitched(db_dt)
    ql = TessQuickLook.__new__(TessQuickLook)
    ql.verbose = False
    ql.target_name = "synthetic"
    ql.precession_search = enabled
    ql.single_event_search = False
    ql.exofop_data = {}
    ql._sector_spans = tuple(stitched_spans if spans is None else spans)
    ql.flat_lc = LightCurve(time=time, flux=flux, flux_err=err)
    ql.tls_results = _TLSResult(
        {
            "SDE": sde,
            "period": PERIOD,
            "T0": T0,
            "duration": _duration(),
            "distinct_transit_count": 30,
        }
    )
    return ql


def test_disabled_search_adds_nothing():
    ql = _ql(enabled=False)
    assert ql.run_precession_search() is None
    assert "precession_search" not in ql.tls_results


def test_single_sector_is_skipped_with_reason():
    ql = _ql(spans=[SectorSpan(1, 0.0, 27.0)])
    result = ql.run_precession_search()
    assert result.reason == "single sector"
    assert ql.tls_results["precession_search"] is True
    assert ql.tls_results["precession_candidate"] is False


def test_non_periodic_signal_is_skipped():
    ql = _ql(sde=4.0)
    result = ql.run_precession_search()
    assert result.reason == "no periodic signal"


def test_search_stores_trend_and_per_sector_arrays():
    ql = _ql(db_dt=0.04)
    result = ql.run_precession_search()
    r = ql.tls_results
    assert result.is_candidate and r["precession_candidate"] is True
    assert r["precession_db_dt"] == pytest.approx(0.04, abs=3 * r["precession_db_dt_err"])
    assert r["precession_sectors"].tolist() == [10, 11, 12, 13]
    assert r["precession_b"].shape == (4,)
    assert r["precession_a_rs_prior"] is False


def test_a_rs_prior_comes_from_stellar_parameters():
    ql = _ql()
    star = {"R_star": 1.0, "R_star_max": 1.05, "M_star": 1.0, "M_star_max": 1.1}
    with patch.object(TessQuickLook, "_stellar_prior_kwargs", return_value=star):
        ql.run_precession_search()
    mean, sigma = ql.precession_a_rs_prior
    assert mean == pytest.approx(215.032 * (PERIOD / 365.25) ** (2 / 3), rel=1e-3)
    assert sigma > 0
    assert ql.tls_results["precession_a_rs_prior"] is True


def test_results_round_trip_through_h5(tmp_path):
    ql = _ql()
    ql.run_precession_search()
    path = tmp_path / "out_tls.h5"
    h5io.save(path, dict(ql.tls_results))
    loaded = h5io.load(path)
    np.testing.assert_allclose(loaded["precession_b"], ql.tls_results["precession_b"])
    assert bool(loaded["precession_candidate"]) is True


def test_saves_precession_figure(tmp_path):
    ql = _ql()
    ql.run_precession_search()
    paths = ql._save_precession_figure(tmp_path / "target")
    assert paths and paths[0].name == "target" + PRECESSION_FIGURE_SUFFIX
    assert paths[0].stat().st_size > 0


def test_candidate_saves_per_sector_folds_figure(tmp_path):
    ql = _ql()
    assert ql.run_precession_search().is_candidate
    paths = ql._save_precession_figure(tmp_path / "target")
    assert [p.name for p in paths] == [
        "target" + PRECESSION_FIGURE_SUFFIX,
        "target" + PRECESSION_FOLDS_SUFFIX,
    ]
    assert all(p.stat().st_size > 0 for p in paths)


def test_non_candidate_skips_folds_figure(tmp_path):
    ql = _ql(db_dt=0.0)
    result = ql.run_precession_search()
    assert not result.is_candidate
    names = [p.name for p in ql._save_precession_figure(tmp_path / "target")]
    assert names == ["target" + PRECESSION_FIGURE_SUFFIX]


def test_no_figure_without_fits(tmp_path):
    ql = _ql(sde=4.0)
    ql.run_precession_search()
    assert ql._save_precession_figure(tmp_path / "target") == []


def test_sector_span_from_lc_reads_dilution_keywords():
    lc = LightCurve(time=np.arange(10.0, 20.0), flux=np.ones(10), flux_err=np.full(10, NOISE))
    lc.meta.update({"SECTOR": 42, "CROWDSAP": 0.97, "FLFRCSAP": 0.88})
    span = sector_span_from_lc(lc)
    assert (span.sector, span.tmin, span.tmax) == (42, 10.0, 19.0)
    assert (span.crowdsap, span.flfrcsap) == (0.97, 0.88)


def test_sector_span_from_lc_without_sector_is_none():
    lc = SimpleNamespace(meta={}, time=SimpleNamespace(value=np.arange(3.0)))
    assert sector_span_from_lc(lc) is None


def test_cli_exposes_precession_flag():
    command = typer.main.get_command(app).commands["run"]
    opts = {opt for param in command.params for opt in param.opts}
    assert "--precession-search" in opts


def test_cli_passes_precession_flag_to_quicklook():
    with patch("quicklook.tql.TessQuickLook") as ql_cls:
        ql_cls.return_value.plot_tql.return_value = None
        result = CliRunner().invoke(
            app,
            [
                "run",
                "--name",
                "TOI-1355",
                "--sector",
                "all",
                "--exptime",
                "120",
                "--precession-search",
                "--save",
            ],
        )
    assert ql_cls.call_args is not None, result.output
    assert ql_cls.call_args.kwargs["precession_search"] is True


def test_record_sector_spans_tolerates_non_iterable_collection():
    ql = TessQuickLook.__new__(TessQuickLook)
    ql._record_sector_spans(object())
    assert ql._sector_spans == ()


def test_record_sector_spans_sorts_by_time():
    def lc(sector, start):
        out = LightCurve(time=np.arange(start, start + 5.0), flux=np.ones(5))
        out.meta["SECTOR"] = sector
        return out

    ql = TessQuickLook.__new__(TessQuickLook)
    ql._record_sector_spans([lc(100, 50.0), lc(99, 10.0)])
    assert [s.sector for s in ql._sector_spans] == [99, 100]


# --- web GUI ---------------------------------------------------------------------------


def test_gui_form_parses_precession_checkbox():
    from quicklook.app.app import _job_dedup_signature, _parse_params

    on = _parse_params({"precession_search": "on"})
    off = _parse_params({})
    assert on["precession_search"] is True and off["precession_search"] is False
    assert _job_dedup_signature(on) != _job_dedup_signature(off)


def test_gui_index_has_precession_checkbox():
    from quicklook.app.app import app as flask_app

    html = flask_app.test_client().get("/").get_data(as_text=True)
    assert 'name="precession_search"' in html


def test_gui_summary_extracts_precession_fields():
    from quicklook.app.tls_summary import extract_fields

    ql = _ql()
    ql.run_precession_search()
    fields = extract_fields(dict(ql.tls_results), "synthetic_s10-13_pdcsap_10c")
    assert fields["prec_candidate"] is True
    assert fields["prec_db_dt"] == pytest.approx(ql.precession_result.db_dt)
    assert fields["prec_sigma"] > 3
