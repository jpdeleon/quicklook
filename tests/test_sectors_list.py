"""`--sectors 66-68,93`: stitch an explicit set of sectors (e.g. consecutive ones)."""

from unittest.mock import patch

import lightkurve as lk
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
        _t = np.arange(0, 1, 200 / 86400)
        time = type("T", (), {"jd": _t, "value": _t})()  # sorted, like a real stitch

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


def _raw_sector_lc(t0, sector, level, seed):
    """A sector downloaded at its own absolute flux level (aperture,
    background, pipeline drift differ sector to sector): ``flux`` and
    ``pdcsap_flux`` both hold the same un-normalized counts, as a real
    SPOC/PDCSAP product does before any selection or normalization."""
    rng = np.random.default_rng(seed)
    t = np.arange(t0, t0 + 25.0, 120 / 86400)
    raw = level * (1 + rng.normal(0, 5e-4, t.size))
    lc = lk.LightCurve(time=t, flux=raw, flux_err=np.full_like(raw, 5e-4 * level))
    lc.meta["SECTOR"] = sector
    lc.meta["AUTHOR"] = "SPOC"
    lc.meta["FRAMETIM"] = 2.0
    lc.meta["NUM_FRM"] = 60
    lc["pdcsap_flux"] = raw.copy()
    lc["pdcsap_flux_err"] = lc["flux_err"].copy()
    return lc


def test_get_lc_normalizes_every_sector_to_the_same_baseline():
    """TOI-1355 and TOI-216: multi-sector runs showed a different flux
    baseline per sector. get_lc() used to reselect the flux column *after*
    stitch() had already normalized a different column, which restored each
    sector's raw, un-normalized level instead of leaving every sector at 1.0.
    """
    rows = [("SPOC", 1, 120), ("SPOC", 2, 120)]
    lcs = lk.LightCurveCollection(
        [_raw_sector_lc(1000.0, 1, 1000.0, 1), _raw_sector_lc(1030.0, 2, 5000.0, 2)]
    )

    class _DownloadableResult(_FakeSearch):
        def __getitem__(self, mask):
            return _DownloadableResult(self.table[mask])

        def download_all(self, **_):
            return lcs

    fake = _DownloadableResult(_search(rows).table)
    ql = _get_lc_instance()
    ql.flux_type, ql.exptime, ql.flatten_method = "pdcsap", None, "biweight"
    ql.sigma_clip_raw = None

    with (
        patch("quicklook.tql.lk.search_lightcurve", return_value=fake),
        patch.object(TessQuickLook, "_download_with_retry", side_effect=lambda f: f()),
    ):
        result = ql.get_lc(author="SPOC", sector=[1, 2], exptime=120)

    in_sector_1 = np.asarray(result.time.value) < 1030.0
    assert np.nanmedian(result.flux.value[in_sector_1]) == pytest.approx(1.0, abs=0.02)
    assert np.nanmedian(result.flux.value[~in_sector_1]) == pytest.approx(1.0, abs=0.02)
