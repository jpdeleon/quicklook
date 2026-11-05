"""`--sector all`: stitch every sector of one pipeline/cadence into one run."""

from unittest.mock import patch

import pytest
from typer.testing import CliRunner

from quicklook.cli.app import app
from quicklook.tql import TessQuickLook

runner = CliRunner()


def _bare_ql(sector):
    ql = object.__new__(TessQuickLook)
    ql.sector = sector
    return ql


@pytest.mark.parametrize(
    "sector, label",
    [(5, "05"), (65, "65"), ([1, 13, 27], "all"), ((1, 2), "all")],
)
def test_sector_label_is_filename_safe(sector, label):
    assert _bare_ql(sector)._sector_label() == label


@pytest.mark.parametrize(
    "sector, expected",
    [(65, 65), ([1, 13, 27], 27), (None, None)],
)
def test_tpf_sector_uses_latest_sector_of_a_stitched_run(sector, expected):
    assert _bare_ql(sector)._tpf_sector == expected


def test_sector_display_summarizes_a_stitched_run():
    assert _bare_ql(65)._sector_display() == "65"
    assert _bare_ql([1, 13, 27])._sector_display() == "1-27 (3 sectors)"


def test_cli_sector_all_requires_exptime():
    result = runner.invoke(app, ["run", "--name", "TOI-1234", "--sector", "all"])
    assert result.exit_code == 2
    assert "--exptime" in result.output


@pytest.mark.parametrize("flag", ["--each-sector", "--each-pipeline"])
def test_cli_sector_all_conflicts_with_each_flags(flag):
    result = runner.invoke(
        app, ["run", "--name", "TOI-1234", "--sector", "all", "--exptime", "200", flag]
    )
    assert result.exit_code == 2
    assert "--sector all" in result.output


def test_cli_rejects_non_numeric_sector():
    result = runner.invoke(app, ["run", "--name", "TOI-1234", "--sector", "abc"])
    assert result.exit_code == 2


@pytest.mark.parametrize("arg, expected", [("all", "all"), ("ALL", "all"), ("27", 27)])
def test_cli_passes_parsed_sector_to_tessquicklook(arg, expected):
    with patch("quicklook.tql.TessQuickLook") as tql:
        runner.invoke(
            app,
            ["run", "--name", "TOI-1234", "--sector", arg, "--exptime", "200", "--save"],
        )
    assert tql.call_args.kwargs["sector"] == expected
