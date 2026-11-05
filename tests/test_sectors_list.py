"""`--sectors 66-68,93`: stitch an explicit set of sectors (e.g. consecutive ones)."""

from unittest.mock import patch

import numpy as np
import pytest
from astropy.table import Table
from typer.testing import CliRunner

from quicklook.cli.app import _parse_sectors, app
from quicklook.exceptions import NoDataError
from quicklook.tql import TessQuickLook, format_sector_ranges

runner = CliRunner()


@pytest.mark.parametrize(
    "sectors, text",
    [
        ([66, 67, 68], "66-68"),
        ([66, 67, 68, 93, 95], "66-68+93+95"),
        ([95, 66, 67, 67], "66-67+95"),
        ([5], "5"),
    ],
)
def test_format_sector_ranges_compacts_consecutive_runs(sectors, text):
    assert format_sector_ranges(sectors) == text


@pytest.mark.parametrize(
    "text, expected",
    [("66-68,93", [66, 67, 68, 93]), ("66, 67", [66, 67]), ("93,66-67,67", [66, 67, 93])],
)
def test_parse_sectors(text, expected):
    assert _parse_sectors(text) == expected


@pytest.mark.parametrize("bad", ["", "a", "68-66", "66-", "0", "-5", "66,,67"])
def test_parse_sectors_rejects_bad_input(bad):
    import typer

    with pytest.raises(typer.BadParameter):
        _parse_sectors(bad)


def _bare_ql(sector, requested):
    ql = object.__new__(TessQuickLook)
    ql.sector = sector
    ql._sectors_requested = requested
    return ql


def test_label_and_display_for_requested_subset():
    ql = _bare_ql([66, 67, 68, 93], [66, 67, 68, 93])
    assert ql._sector_label() == "66-68+93"
    assert ql._sector_display() == "66-68, 93 (4 sectors)"


def test_label_is_all_for_full_stitch():
    assert _bare_ql([1, 13, 27], None)._sector_label() == "all"


@pytest.mark.parametrize(
    "extra",
    [
        ["--sector", "27"],
        ["--sectors", "66-67", "--each-sector"],
        ["--each-pipeline"],
    ],
)
def test_cli_sectors_conflicts(extra):
    args = ["run", "--name", "TOI-1234", "--sectors", "66-67", "--exptime", "200"] + extra
    result = runner.invoke(app, args)
    assert result.exit_code == 2
    assert "--sectors" in result.output


def test_cli_sectors_requires_exptime():
    result = runner.invoke(app, ["run", "--name", "TOI-1234", "--sectors", "66-67"])
    assert result.exit_code == 2
    assert "--exptime" in result.output


def test_cli_sectors_needs_at_least_two_sectors():
    result = runner.invoke(
        app, ["run", "--name", "TOI-1234", "--sectors", "66", "--exptime", "200"]
    )
    assert result.exit_code == 2
    assert "at least two" in result.output.lower()


def test_cli_passes_sector_list_to_tessquicklook():
    with patch("quicklook.tql.TessQuickLook") as tql:
        runner.invoke(
            app,
            ["run", "--name", "TOI-1234", "--sectors", "66-68", "--exptime", "200", "--save"],
        )
    assert tql.call_args.kwargs["sector"] == [66, 67, 68]


class _FakeSearch:
    """Minimal stand-in for a lightkurve SearchResult."""

    def __init__(self, table):
        self.table = table

    def __len__(self):
        return len(self.table)

    def __getitem__(self, mask):
        return _FakeSearch(self.table[mask])


def _search(rows):
    return _FakeSearch(
        Table(
            {
                "provenance_name": [r[0] for r in rows],
                "mission": [f"TESS Sector {r[1]:02d}" for r in rows],
                "t_exptime": [float(r[2]) for r in rows],
            }
        )
    )


def _get_lc_instance():
    ql = object.__new__(TessQuickLook)
    ql.query_name = "TOI-173"
    ql.verbose = False
    ql.quality_bitmask = "default"
    return ql


def test_get_lc_reports_requested_sectors_missing_at_that_cadence():
    rows = [("QLP", 66, 200), ("QLP", 67, 200), ("QLP", 1, 1800)]
    with patch("quicklook.tql.lk.search_lightcurve", return_value=_search(rows)):
        with pytest.raises(NoDataError, match="68"):
            _get_lc_instance().get_lc(author="QLP", sector=[66, 67, 68], exptime=200)


def test_get_lc_only_downloads_requested_sectors():
    rows = [("QLP", s, 200) for s in (66, 67, 68, 93)] + [("QLP", 1, 1800)]
    captured = {}
    ql = _get_lc_instance()

    def fake_download(callable_):
        return callable_()

    class _LC:
        meta = {"AUTHOR": "QLP"}
        time = type("T", (), {"jd": np.arange(0, 1, 200 / 86400)})()

        def select_flux(self, *_):
            return self

    class _Stitchable:
        def stitch(self):
            return _LC()

    class _Result(_FakeSearch):
        def __getitem__(self, mask):
            res = _Result(self.table[mask])
            captured["sectors"] = sorted(int(m.split()[-1]) for m in res.table["mission"])
            return res

        def download_all(self, **_):
            return _Stitchable()

    fake = _Result(_search(rows).table)
    ql.flux_type, ql.exptime, ql.flatten_method = "pdcsap", None, "biweight"
    with (
        patch("quicklook.tql.lk.search_lightcurve", return_value=fake),
        patch.object(TessQuickLook, "_download_with_retry", side_effect=fake_download),
    ):
        try:
            ql.get_lc(author="QLP", sector=[66, 67, 68], exptime=200)
        except Exception:
            pass  # only the sector selection is under test
    assert captured["sectors"] == [66, 67, 68]
    assert ql.sector == [66, 67, 68]
