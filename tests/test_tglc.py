import os
import numpy as np
import pytest
import requests
import lightkurve as lk
from astroquery.exceptions import RemoteServiceError

from quicklook.tglc import get_tglc_lc

_NETWORK_ERRORS = (
    RemoteServiceError,
    requests.exceptions.HTTPError,
    requests.exceptions.ConnectionError,
    ConnectionError,
    TimeoutError,
)


@pytest.fixture
def tglc_inputs():
    """Bright, isolated, single-target config for a fast TGLC smoke run.

    pi Men (TIC 261136679) has short-cadence SPOC in sector 1 and enough
    FFI coverage to exercise the ePSF pipeline end-to-end. A 15x15 cutout
    is below TGLC's recommended 25x25 (the module will warn) but keeps the
    TESScut payload small enough to stay under MAST timeouts in CI.
    """
    return {
        "target_name": "TIC 261136679",
        "sector": 1,
        "size": 15,
        "limit_mag": 13,
        "verbose": False,
    }


@pytest.mark.network
def test_get_tglc_lc(tglc_inputs):
    """End-to-end TGLC ePSF pipeline smoke test.

    Regression coverage for three bugs introduced in commit efbe476:
    - lowercase 'designation' column from MAST gaiadr3
    - missing f-string prefix on sector_{N}_x/y column writes
    - swapped isinstance(Source, source) argument order
    """
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
