import os
import numpy as np
import pytest
import requests
import lightkurve as lk
from astroquery.exceptions import RemoteServiceError

from astropy.table import Table

from quicklook.tglc import find_target_index, fit_lc, get_tglc_lc

_NETWORK_ERRORS = (
    RemoteServiceError,
    requests.exceptions.HTTPError,
    requests.exceptions.ConnectionError,
    ConnectionError,
    TimeoutError,
)


@pytest.fixture
def tglc_inputs():
    """Bright, isolated, single-target config for a TGLC smoke run.

    pi Men (TIC 261136679) has FFI coverage in sector 1. Size 50 is the
    minimum allowed cutout (enforced by get_tglc_lc).
    """
    return {
        "target_name": "TIC 261136679",
        "sector": 1,
        "size": 50,
        "limit_mag": 13,
        "verbose": False,
    }


@pytest.mark.network
def test_get_tglc_lc(tglc_inputs):
    """End-to-end TGLC ePSF pipeline smoke test."""
    try:
        lc = get_tglc_lc(**tglc_inputs)
    except _NETWORK_ERRORS as e:
        if os.getenv("CI", "false") == "true":
            pytest.xfail(f"Network failure (not a bug): {e}")
        raise

    assert isinstance(lc, lk.LightCurve)
    assert len(lc) > 0
    assert np.isfinite(np.nanmedian(lc.flux.value))
    assert lc.sector == tglc_inputs["sector"]


@pytest.mark.parametrize("bad_size", [0, 1, 25, 49, 100, 150, 200])
def test_get_tglc_lc_rejects_out_of_range_size(bad_size, tglc_inputs):
    """Sizes outside [50, 99] must raise ValueError before any network call."""
    inputs = {**tglc_inputs, "size": bad_size}
    with pytest.raises(ValueError, match="cutsize must be between 50 and 99"):
        get_tglc_lc(**inputs)


@pytest.mark.parametrize("good_size", [50, 75, 99])
def test_get_tglc_lc_accepts_in_range_size(good_size, tglc_inputs, monkeypatch):
    """Sizes 50-99 must pass validation. Stub the pipeline so this is offline."""
    inputs = {**tglc_inputs, "size": good_size}
    sentinel = RuntimeError("stub-pipeline-reached")

    def _boom(*args, **kwargs):
        raise sentinel

    # First thing get_tglc_lc does after validation is build a Source_cut.
    # If validation passes, we should hit this stub (not a ValueError).
    monkeypatch.setattr("quicklook.tglc.Source_cut", _boom)
    with pytest.raises(RuntimeError, match="stub-pipeline-reached"):
        get_tglc_lc(**inputs)


def test_fit_lc_prior_requires_coordinate_arrays():
    """When prior is set, fit_lc must be given x_all and y_all."""
    with pytest.raises(ValueError, match="requires x_all and y_all"):
        fit_lc(
            A=np.zeros((1, 1)),
            source=object(),
            star_info=None,
            x=0.0,
            y=0.0,
            star_num=0,
            e_psf=np.zeros((1, 1)),
            near_edge=False,
            prior=0.1,
        )


def test_get_tglc_lc_strips_hdulist_from_legacy_cache(tglc_inputs, tmp_path, monkeypatch):
    """A cache pickle written before the hdulist guard existed must still be
    memory-safe once loaded: get_tglc_lc must strip hdulist after pickle.load,
    not just on the fresh-Source_cut-construction path."""
    import pickle
    from types import SimpleNamespace

    cache_dir = str(tmp_path)
    slug = tglc_inputs["target_name"].replace(" ", "")
    sector_tag = f"s{tglc_inputs['sector']}"
    source_dir = os.path.join(cache_dir, "source")
    os.makedirs(source_dir, exist_ok=True)
    source_pkl = os.path.join(
        source_dir, f"source_{slug}_{sector_tag}_size{tglc_inputs['size']}.pkl"
    )

    # Simulate a legacy cache: a Source_cut-like object that still carries a
    # raw hdulist attribute (a real one would be multi-GB for a 200s cadence
    # sector; a sentinel is enough to prove it gets dropped).
    legacy_source = SimpleNamespace(hdulist=["stand-in-for-a-huge-raw-hdulist"])
    with open(source_pkl, "wb") as fh:
        pickle.dump(legacy_source, fh)

    captured = {}

    def _capture_and_boom(source, *args, **kwargs):
        captured["source"] = source
        raise RuntimeError("stub-pipeline-reached")

    # get_psf(source) is the first thing get_tglc_lc does after resolving
    # source (from cache or fresh build); stubbing it lets us inspect the
    # loaded source without running the rest of the pipeline.
    monkeypatch.setattr("quicklook.tglc.get_psf", _capture_and_boom)

    with pytest.raises(RuntimeError, match="stub-pipeline-reached"):
        get_tglc_lc(**tglc_inputs, cache_dir=cache_dir)

    assert captured["source"].hdulist is None


def test_fit_lc_prior_none_keeps_default_behavior():
    """prior=None (default) must not require x_all/y_all and must not raise
    the new ValueError. We trigger a different failure later (missing data)
    to confirm validation passed."""
    # Build a minimal source that has the attributes fit_lc touches before
    # the psf-fit block. We only need to reach past the prior-branch guard.
    with pytest.raises(Exception) as excinfo:
        fit_lc(
            A=np.zeros((1, 1)),
            source=object(),
            star_info=None,
            x=0.0,
            y=0.0,
            star_num=0,
            e_psf=np.zeros((1, 1)),
            near_edge=False,
        )
    # The failure must NOT be the prior-coord-array guard.
    assert "requires x_all and y_all" not in str(excinfo.value)


# --- target identification (regression for the gaia[0] bug) ------------------


def _fake_source(with_coords=True):
    """Stand-in Source_cut: gaia sorted by brightness, target is NOT row 0."""
    from types import SimpleNamespace

    gaia = Table(
        {
            "designation": ["bright-neighbour", "target", "faint"],
            "tess_mag": [8.5, 10.6, 14.0],
            "ra": [68.174, 68.416, 68.417],
            "dec": [-51.755, -51.957, -51.950],
            "sector_3_x": [10.0, 49.2, 50.5],
            "sector_3_y": [80.0, 48.9, 60.0],
        }
    )
    src = SimpleNamespace(gaia=gaia, size=99, sector=3)
    if with_coords:
        src.target_ra, src.target_dec = 68.4161, -51.9574
    return src


def test_find_target_index_matches_target_coordinates_not_brightest():
    assert find_target_index(_fake_source()) == 1


def test_find_target_index_falls_back_to_cutout_centre_for_legacy_cache():
    assert find_target_index(_fake_source(with_coords=False)) == 1


def test_find_target_index_empty_gaia_raises():
    src = _fake_source()
    src.gaia = src.gaia[:0]
    with pytest.raises(ValueError, match="No Gaia stars"):
        find_target_index(src)


@pytest.mark.parametrize("bad_prior", [0, -0.1, 0.0])
def test_get_tglc_lc_rejects_non_positive_prior(bad_prior, tglc_inputs):
    """prior must be a positive float or None; reject before any network call."""
    with pytest.raises(ValueError, match="prior must be a positive float or None"):
        get_tglc_lc(**tglc_inputs, prior=bad_prior)


@pytest.mark.parametrize("good_prior", [None, 0.1, 1.0])
def test_get_tglc_lc_accepts_valid_prior(good_prior, tglc_inputs, monkeypatch):
    """Valid priors pass validation and reach the pipeline (stubbed offline)."""

    def _boom(*args, **kwargs):
        raise RuntimeError("stub-pipeline-reached")

    monkeypatch.setattr("quicklook.tglc.Source_cut", _boom)
    with pytest.raises(RuntimeError, match="stub-pipeline-reached"):
        get_tglc_lc(**tglc_inputs, prior=good_prior, use_cache=False)


def _flagged_lc():
    """Light curve with get_tglc_lc-style TESS and TGLC flag columns."""
    lc = lk.TessLightCurve(time=np.arange(6.0), flux=np.ones(6))
    tess = np.array([0, 128, 2048, 0, 0, 16384], dtype=np.int32)  # manual, straylight, badcal
    lc["tess_flags"] = tess
    lc["quality"] = tess
    lc["tglc_flags"] = np.array([0, 0, 0, 1, 0, 0], dtype=np.int16)
    return lc


@pytest.mark.parametrize(
    "bitmask,expected",
    [
        ("none", [True] * 6),
        # default: manual exclude + bad calibration, not straylight; TGLC flag always
        ("default", [True, False, True, False, True, False]),
        ("hard", [True, False, False, False, True, False]),
        ("hardest", [True, False, False, False, True, False]),
        (128, [True, False, True, False, True, True]),
    ],
)
def test_tglc_quality_mask_combines_tess_bitmask_and_tglc_flag(bitmask, expected):
    from quicklook.tglc import tglc_quality_mask

    assert tglc_quality_mask(_flagged_lc(), bitmask).tolist() == expected


def test_tglc_quality_mask_without_flag_columns_keeps_everything():
    from quicklook.tglc import tglc_quality_mask

    lc = lk.LightCurve(time=np.arange(3.0), flux=np.ones(3))
    assert tglc_quality_mask(lc, "hard").tolist() == [True, True, True]
