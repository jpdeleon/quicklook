"""Offline tests for quicklook.utils."""

import numpy as np
import pytest
from astropy.table import Table


def _fake_sector_table(rows):
    """Build an astropy Table matching Tesscut.get_sectors() output shape."""
    return Table(
        rows=rows,
        names=("sectorName", "sector", "camera", "ccd"),
        dtype=("U20", "int64", "int64", "int64"),
    )


def test_get_tglc_sectors_filters_synthetic_mast_codes(monkeypatch):
    """MAST sometimes returns synthetic sectorNames (e.g. tess-s1751-1-3 for
    TOI-5169) whose ``sector`` field is not a real observing sector. These
    must be dropped so each-sector queueing doesn't enqueue invalid jobs.
    Regression for the TOI-5169_s1751 bug.
    """
    from astropy.coordinates import SkyCoord
    from astroquery.mast import Tesscut
    from quicklook import utils

    fake_table = _fake_sector_table(
        [
            ("tess-s0045-4-2", 45, 4, 2),
            ("tess-s0046-2-1", 46, 2, 1),
            ("tess-s0072-3-3", 72, 3, 3),
            ("tess-s1751-1-3", 1751, 1, 3),
        ]
    )
    monkeypatch.setattr(SkyCoord, "from_name", classmethod(lambda cls, name: object()))
    monkeypatch.setattr(Tesscut, "get_sectors", staticmethod(lambda coordinates: fake_table))

    sectors = utils._get_tglc_sectors_via_tesscut("TOI-5169")
    assert sectors == [45, 46, 72]
    assert 1751 not in sectors


def test_get_tglc_sectors_keeps_real_sectors(monkeypatch):
    """All sectors in the legal [1, 300] range must be preserved and sorted."""
    from astropy.coordinates import SkyCoord
    from astroquery.mast import Tesscut
    from quicklook import utils

    fake_table = _fake_sector_table(
        [
            ("tess-s0050-1-1", 50, 1, 1),
            ("tess-s0001-2-2", 1, 2, 2),
            ("tess-s0300-3-3", 300, 3, 3),  # boundary inclusive
            ("tess-s0100-4-4", 100, 4, 4),
        ]
    )
    monkeypatch.setattr(SkyCoord, "from_name", classmethod(lambda cls, name: object()))
    monkeypatch.setattr(Tesscut, "get_sectors", staticmethod(lambda coordinates: fake_table))

    sectors = utils._get_tglc_sectors_via_tesscut("Dummy")
    assert sectors == [1, 50, 100, 300]


def test_get_tglc_sectors_drops_out_of_range(monkeypatch):
    """Sectors <1 or >300 must be filtered out."""
    from astropy.coordinates import SkyCoord
    from astroquery.mast import Tesscut
    from quicklook import utils

    fake_table = _fake_sector_table(
        [
            ("tess-s0000-1-1", 0, 1, 1),
            ("tess-s0017-2-2", 17, 2, 2),
            ("tess-s0301-3-3", 301, 3, 3),
            ("tess-s9999-4-4", 9999, 4, 4),
        ]
    )
    monkeypatch.setattr(SkyCoord, "from_name", classmethod(lambda cls, name: object()))
    monkeypatch.setattr(Tesscut, "get_sectors", staticmethod(lambda coordinates: fake_table))

    sectors = utils._get_tglc_sectors_via_tesscut("Dummy")
    assert sectors == [17]


def test_get_tglc_sectors_empty_table(monkeypatch):
    """No coverage -> empty list, no exception."""
    from astropy.coordinates import SkyCoord
    from astroquery.mast import Tesscut
    from quicklook import utils

    empty = Table(
        names=("sectorName", "sector", "camera", "ccd"),
        dtype=("U20", "int64", "int64", "int64"),
    )
    monkeypatch.setattr(SkyCoord, "from_name", classmethod(lambda cls, name: object()))
    monkeypatch.setattr(Tesscut, "get_sectors", staticmethod(lambda coordinates: empty))

    assert utils._get_tglc_sectors_via_tesscut("Dummy") == []
