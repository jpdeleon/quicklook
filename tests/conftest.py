"""
Pytest configuration file for the quicklook package.
"""

import os
import pytest
import tempfile
from pathlib import Path


@pytest.fixture(scope="session")
def temp_dir():
    """Create a temporary directory for test outputs."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        yield tmp_dir


@pytest.fixture
def mock_light_curve():
    """Create a mock light curve for testing."""
    import numpy as np
    import lightkurve as lk

    # Create a simple sine wave light curve
    time = np.linspace(0, 27, 1000)
    flux = 1.0 + 0.1 * np.sin(2 * np.pi * time / 5.0)
    flux_err = np.ones_like(flux) * 0.01

    # Add a transit-like dip
    transit_mask = (time % 3.0 > 2.9) | (time % 3.0 < 0.1)
    flux[transit_mask] -= 0.02

    # Create a LightCurve object
    lc = lk.LightCurve(time=time, flux=flux, flux_err=flux_err)
    lc.meta = {
        "MISSION": "TESS",
        "SECTOR": 1,
        "AUTHOR": "SPOC",
        "EXPOSURE": 120,
    }

    return lc


@pytest.fixture
def skip_if_no_internet(request):
    """Skip a test if there is no internet connection."""
    import socket

    try:
        # Try to connect to a well-known external server
        socket.create_connection(("www.google.com", 80), timeout=1)
    except OSError:
        pytest.skip("No internet connection available")
