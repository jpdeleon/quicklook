"""MAST TGLC HLSP light curves must honour quality_bitmask='hardest'.

lightkurve's TGLC reader overflows on 'hardest' because ``tess_flags`` is
int16 (see LIGHTKURVE_TGLC_BUG_REPORT.txt). quicklook downloads HLSP files
unmasked and applies the bitmask itself.
"""

import lightkurve as lk
import numpy as np
import pytest

from quicklook.tql import _apply_tglc_hlsp_mask, _download_bitmask


def _make_tglc_lc():
    n = 6
    lc = lk.LightCurve(time=np.arange(n, dtype=float), flux=np.ones(n))
    # bit 1 (attitude tweak, in default), bit 2048 (stray light, hardest only)
    lc["tess_flags"] = np.array([0, 1, 2048, 0, 0, 0], dtype=np.int16)
    lc["tglc_flags"] = np.array([0, 0, 0, 1, 0, 0], dtype=np.int16)
    lc.meta["AUTHOR"] = "TGLC"
    lc.meta["SECTOR"] = 11
    return lc


@pytest.mark.unit
def test_download_bitmask_is_none_for_tglc():
    assert _download_bitmask("TGLC", "hardest") == "none"
    assert _download_bitmask("tglc", "hard") == "none"


@pytest.mark.unit
def test_download_bitmask_passthrough_for_other_pipelines():
    assert _download_bitmask("SPOC", "hardest") == "hardest"
    assert _download_bitmask("QLP", "default") == "default"


@pytest.mark.unit
def test_hardest_mask_on_int16_flags_does_not_overflow():
    lc = _make_tglc_lc()

    out = _apply_tglc_hlsp_mask(lc, "hardest")

    # drops tess bit 1, tess bit 2048 and the TGLC background flag
    assert len(out) == 3
    assert out.meta["SECTOR"] == 11


@pytest.mark.unit
def test_default_mask_keeps_hardest_only_cadences():
    lc = _make_tglc_lc()

    out = _apply_tglc_hlsp_mask(lc, "default")

    assert 2.0 in np.asarray(out.time.value)
    assert len(out) == 4


@pytest.mark.unit
def test_none_mask_keeps_everything():
    lc = _make_tglc_lc()

    assert len(_apply_tglc_hlsp_mask(lc, "none")) == len(lc)
