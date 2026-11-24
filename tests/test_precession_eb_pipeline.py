"""EB precession branch wired into TessQuickLook, the figure, h5 summary, CLI and GUI.

Offline: the pipeline object is built with ``__new__`` and a synthetic stitched EB.
"""

import os

os.environ.setdefault("MPLBACKEND", "Agg")

from unittest.mock import patch

import numpy as np
import pytest
from lightkurve import LightCurve
from typer.testing import CliRunner

from quicklook import h5io
from quicklook.cli.app import app
from quicklook.precession import PrecessionResult
from quicklook.precession_eb import EBPrecessionResult
from quicklook.precession_pipeline import PRECESSION_FIGURE_SUFFIX, is_eclipsing_binary_otype
from quicklook.tql import TessQuickLook, _TLSResult
from test_precession_eb import PERIOD, T0, _duration, _stitched


def _ql(di_dt=1.0, mode="eb", otype=None):
    time, flux, err, spans = _stitched(di_dt)
    ql = TessQuickLook.__new__(TessQuickLook)
    ql.verbose = False
    ql.target_name = "synthetic EB"
    ql.precession_search = True
    ql.precession_mode = mode
    ql.simbad_obj_type = otype
    ql.single_event_search = False
    ql.exofop_data = {}
    ql._sector_spans = tuple(spans)
    ql.flat_lc = LightCurve(time=time, flux=flux, flux_err=err)
    ql.tls_results = _TLSResult(
        {
            "SDE": 30.0,
            "period": PERIOD,
            "T0": T0,
            "duration": _duration(),
            "distinct_transit_count": 25,
        }
    )
    return ql


@pytest.mark.parametrize(
    "otype, expected",
    [
        ("EclBin", True),
        ("Eclipsing Binary", True),
        ("Eclipsing Binary of Algol type", True),
        ("EB*", True),
        ("Star showing Eclipses by its Planet", False),
        ("Star", False),
        (None, False),
    ],
)
def test_eclipsing_binary_otype(otype, expected):
    assert is_eclipsing_binary_otype(otype) is expected


def test_auto_mode_uses_simbad_type():
    assert _ql(mode="auto", otype="EclBin").resolve_precession_mode() == "eb"
    assert _ql(mode="auto", otype="Star").resolve_precession_mode() == "planet"
    assert _ql(mode="planet", otype="EclBin").resolve_precession_mode() == "planet"


def test_unknown_mode_raises():
    with pytest.raises(ValueError):
        _ql(mode="triple").resolve_precession_mode()


def test_eb_search_stores_inclination_trend_and_binary():
    ql = _ql()
    result = ql.run_precession_search()
    r = ql.tls_results
    assert isinstance(result, EBPrecessionResult)
    assert result.is_candidate and r["precession_candidate"] is True
    assert r["precession_mode"] == "eb"
    assert r["precession_di_dt"] == pytest.approx(1.0, abs=3 * r["precession_di_dt_err"] + 0.05)
    assert r["precession_db_dt"] < 0
    assert r["precession_inc"].shape == (4,)
    assert r["precession_has_secondary"].tolist() == [1.0] * 4
    assert 0 < r["precession_eb_k"] < 3
    assert "precession_a_rs_prior" not in r


def test_planet_mode_keeps_planet_branch():
    ql = _ql(mode="planet")
    assert isinstance(ql.run_precession_search(), PrecessionResult)
    assert ql.tls_results["precession_mode"] == "planet"


def test_eb_results_round_trip_through_h5(tmp_path):
    ql = _ql()
    ql.run_precession_search()
    path = tmp_path / "eb_tls.h5"
    h5io.save(path, dict(ql.tls_results))
    loaded = h5io.load(path)
    np.testing.assert_allclose(loaded["precession_inc"], ql.tls_results["precession_inc"])
    assert str(loaded["precession_mode"]) == "eb"


def test_saves_eb_precession_figure(tmp_path):
    ql = _ql()
    ql.run_precession_search()
    paths = ql._save_precession_figure(tmp_path / "eb")
    assert paths and paths[0].name == "eb" + PRECESSION_FIGURE_SUFFIX
    assert paths[0].stat().st_size > 0


def test_cli_passes_precession_mode():
    with patch("quicklook.tql.TessQuickLook") as ql_cls:
        ql_cls.return_value.plot_tql.return_value = None
        result = CliRunner().invoke(
            app,
            ["run", "--name", "V907 Sco", "--sector", "all", "--exptime", "120",
             "--precession-search", "--precession-mode", "EB", "--save"],
        )  # fmt: skip
    assert ql_cls.call_args is not None, result.output
    assert ql_cls.call_args.kwargs["precession_mode"] == "eb"


def test_cli_rejects_unknown_precession_mode():
    result = CliRunner().invoke(app, ["run", "--name", "X", "--precession-mode", "triple"])
    assert result.exit_code != 0
    assert "precession_mode" in result.output


def test_gui_form_parses_precession_mode():
    from quicklook.app.app import _job_dedup_signature, _parse_params

    eb = _parse_params({"precession_search": "on", "precession_mode": "EB"})
    planet = _parse_params({"precession_search": "on", "precession_mode": "planet"})
    assert eb["precession_mode"] == "eb"
    assert _parse_params({"precession_mode": "bogus"})["precession_mode"] == "auto"
    assert _job_dedup_signature(eb) != _job_dedup_signature(planet)


def test_gui_index_has_precession_mode_select():
    from quicklook.app.app import app as flask_app

    html = flask_app.test_client().get("/").get_data(as_text=True)
    assert 'name="precession_mode"' in html


def test_gui_summary_extracts_eb_fields():
    from quicklook.app.tls_summary import extract_fields

    ql = _ql()
    ql.run_precession_search()
    fields = extract_fields(dict(ql.tls_results), "synthetic_s20-23_pdcsap_10c")
    assert fields["prec_mode"] == "eb"
    assert fields["prec_di_dt"] == pytest.approx(ql.precession_result.di_dt)
    assert fields["prec_candidate"] is True


@pytest.mark.parametrize("reason", ["single sector", "no periodic signal"])
def test_skipped_auto_eb_preserves_mode(reason):
    ql = _ql(mode="auto", otype="EclBin")
    if reason == "single sector":
        ql._sector_spans = ql._sector_spans[:1]
    else:
        ql.tls_results["SDE"] = 0.0
    result = ql.run_precession_search()
    assert isinstance(result, EBPrecessionResult)
    assert result.reason == reason
    assert ql.tls_results["precession_mode"] == "eb"
