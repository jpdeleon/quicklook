"""Detrending must refuse light curves whose flux is not positive.

wotan treats non-positive flux as masked. When a whole sliding window is
non-positive, the biweight location estimator receives an empty array, its
centre becomes NaN, and the numba ``fastmath`` convergence loop never exits.
A zero-centred TASOC light curve (WASP-12, sector 3) hung ``ql run
--each-pipeline`` for 4.5 h this way in the benchmark suite (12_cross_pipeline).
"""

import lightkurve as lk
import numpy as np
import pytest

from quicklook.tql import (
    MAX_NONPOSITIVE_FLUX_FRAC,
    TessQuickLook,
    _check_positive_flux,
)

CADENCE = 2 / 1440


def _flux(offset, seed=0, n=2000):
    rng = np.random.default_rng(seed)
    return offset + rng.normal(0, 1e-3, n)


def _bare_ql(flux, method="biweight"):
    ql = TessQuickLook.__new__(TessQuickLook)
    t = np.arange(flux.size) * CADENCE
    ql.raw_lc = lk.LightCurve(time=t, flux=flux, flux_err=np.full_like(flux, 1e-3))
    ql.flatten_method = method
    ql.window_length = 0.5
    ql.verbose = False
    return ql


@pytest.mark.unit
def test_normalised_flux_passes():
    _check_positive_flux(_flux(1.0))


@pytest.mark.unit
def test_zero_centred_flux_raises():
    with pytest.raises(ValueError, match="zero-centred or negative"):
        _check_positive_flux(_flux(0.0))


@pytest.mark.unit
def test_negative_median_flux_raises():
    with pytest.raises(ValueError, match="median"):
        _check_positive_flux(_flux(-5.0))


@pytest.mark.unit
def test_few_negative_outliers_are_tolerated():
    flux = _flux(1.0)
    n_bad = int(0.5 * MAX_NONPOSITIVE_FLUX_FRAC * flux.size)
    flux[:n_bad] = -1.0

    _check_positive_flux(flux)


@pytest.mark.unit
def test_too_many_nonpositive_points_raise():
    flux = _flux(1.0)
    n_bad = int(2 * MAX_NONPOSITIVE_FLUX_FRAC * flux.size)
    flux[:n_bad] = 0.0

    with pytest.raises(ValueError, match="non-positive"):
        _check_positive_flux(flux)


@pytest.mark.unit
def test_nans_are_ignored():
    flux = _flux(1.0)
    flux[::3] = np.nan

    _check_positive_flux(flux)


@pytest.mark.unit
def test_all_nan_flux_raises():
    with pytest.raises(ValueError, match="finite"):
        _check_positive_flux(np.full(100, np.nan))


@pytest.mark.unit
@pytest.mark.parametrize("method", ["biweight", "notch", "locor"])
def test_flatten_raw_lc_rejects_zero_centred_flux_before_detrending(method):
    ql = _bare_ql(_flux(0.0), method=method)

    with pytest.raises(ValueError, match="zero-centred or negative"):
        ql.flatten_raw_lc()
