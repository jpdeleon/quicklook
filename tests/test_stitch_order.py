"""Multi-sector stitching must return a time-ordered light curve.

MAST returned TOI-1227's sectors as [100, 99]; ``LightCurveCollection.stitch``
keeps that order, so the stitched light curve jumped back in time and the
rotation periodogram crashed (empty lightkurve frequency grid). Seen in the
benchmark suite (06_young).
"""

import lightkurve as lk
import numpy as np
import pytest

from quicklook.tql import TessQuickLook, _sort_by_time

CADENCE = 2 / 1440


def _sector_lc(start, sector, seed):
    rng = np.random.default_rng(seed)
    t = np.arange(start, start + 25.0, CADENCE)
    f = 1 + 1e-3 * np.sin(2 * np.pi * t / 3.0) + rng.normal(0, 5e-4, t.size)
    lc = lk.LightCurve(time=t, flux=f, flux_err=np.full_like(f, 5e-4))
    lc.meta["SECTOR"] = sector
    return lc


class _FakeSearch:
    """Stands in for a lightkurve SearchResult returning sectors 100, 99."""

    def __init__(self, lcs):
        self._lcs = lcs

    def download_all(self, quality_bitmask=None):
        return lk.LightCurveCollection(self._lcs)


def _bare_ql():
    ql = TessQuickLook.__new__(TessQuickLook)
    ql.quality_bitmask = "default"
    return ql


@pytest.mark.unit
def test_stitch_masked_sorts_sectors_returned_out_of_order():
    s100, s99 = _sector_lc(4074.0, 100, 1), _sector_lc(4046.0, 99, 2)

    lc = _bare_ql()._stitch_masked(_FakeSearch([s100, s99]), "SPOC")

    t = lc.time.value
    assert np.all(np.diff(t) > 0)
    assert len(lc) == len(s100) + len(s99)
    assert t[0] == pytest.approx(4046.0)


@pytest.mark.unit
def test_periodogram_works_on_out_of_order_stitch():
    # The crash in plot_tql: to_periodogram on the stitched (unsorted) light curve.
    lc = _bare_ql()._stitch_masked(
        _FakeSearch([_sector_lc(4074.0, 100, 1), _sector_lc(4046.0, 99, 2)]), "SPOC"
    )

    pg = lc.to_periodogram(method="lombscargle")

    assert pg.period_at_max_power.value == pytest.approx(3.0, rel=0.05)


@pytest.mark.unit
def test_sort_by_time_keeps_sorted_input_unchanged():
    lc = _sector_lc(1000.0, 1, 3)

    assert _sort_by_time(lc) is lc


@pytest.mark.unit
def test_sort_by_time_carries_flux_with_time():
    lc = lk.LightCurve(time=[3.0, 1.0, 2.0], flux=[30.0, 10.0, 20.0])

    out = _sort_by_time(lc)

    assert out.time.value.tolist() == [1.0, 2.0, 3.0]
    assert out.flux.value.tolist() == [10.0, 20.0, 30.0]
