import os
import numpy as np
import pytest
import requests
import lightkurve as lk
from astroquery.exceptions import RemoteServiceError

from quicklook.tglc import fit_lc, get_tglc_lc

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
