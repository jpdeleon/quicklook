"""Tests for the local-vs-MAST TGLC benchmark (quicklook.benchmark_tglc)."""

import os

import lightkurve as lk
import numpy as np
import pytest

from quicklook.benchmark_tglc import (
    COLUMNS,
    compare_light_curves,
    compare_series,
    mast_good_mask,
    match_cadences,
    point_to_point_scatter,
    robust_std,
)

CADENCE = 1800 / 86400


def _make_lc(time, flux_by_col, flags=None):
    """Build a TessLightCurve carrying TGLC-style photometry columns."""
    lc = lk.TessLightCurve(time=time, flux=flux_by_col["cal_aper_flux"])
    for col, values in flux_by_col.items():
        lc[col] = values
    if flags is not None:
        lc["tess_flags"] = flags
        lc["tglc_flags"] = np.zeros(len(time), dtype=np.int16)
    return lc


def _synthetic_pair(noise_local=1e-4, seed=0, n=500):
    """MAST-like and local-like light curves sharing a common signal."""
    rng = np.random.default_rng(seed)
    t = 1325.0 + np.arange(n) * CADENCE
    signal = 1 + 0.01 * np.sin(2 * np.pi * t / 2.0)
    signal[200:210] -= 0.02  # a transit-like dip
    shared = signal + rng.normal(0, 1e-3, n)
    mast = {c: shared * (5000 if "cal" not in c else 1) for c in COLUMNS}
    local = {
        c: (shared + rng.normal(0, noise_local, n)) * (5100 if "cal" not in c else 1)
        for c in COLUMNS
    }
    return t, local, mast


# --- robust statistics -----------------------------------------------------


def test_robust_std_matches_gaussian_sigma():
    x = np.random.default_rng(1).normal(0, 2.0, 20000)
    assert robust_std(x) == pytest.approx(2.0, rel=0.03)


def test_robust_std_ignores_nan_and_empty():
    assert np.isnan(robust_std(np.array([np.nan, np.nan])))
    assert robust_std(np.array([1.0, np.nan, 1.0, 1.0])) == 0.0


def test_point_to_point_scatter_insensitive_to_trend():
    rng = np.random.default_rng(2)
    noise = rng.normal(0, 1e-3, 5000)
    trend = np.linspace(0, 1, 5000)
    assert point_to_point_scatter(noise + trend) == pytest.approx(1e-3, rel=0.05)


def test_point_to_point_scatter_too_short_is_nan():
    assert np.isnan(point_to_point_scatter(np.array([1.0, 2.0])))


# --- cadence matching ------------------------------------------------------


def test_match_cadences_pairs_nearest_within_tolerance():
    t_local = np.array([0.0, 1.0, 2.0, 3.0])
    t_mast = np.array([1.0001, 2.9999, 10.0])
    il, im = match_cadences(t_local, t_mast, tol=0.01)
    assert il.tolist() == [1, 3]
    assert im.tolist() == [0, 1]


def test_match_cadences_handles_unsorted_local_times():
    t_local = np.array([3.0, 0.0, 2.0, 1.0])
    t_mast = np.array([2.0, 0.0])
    il, im = match_cadences(t_local, t_mast, tol=0.01)
    assert len(il) == 2
    assert np.allclose(t_local[il], t_mast[im])


def test_match_cadences_uses_each_local_cadence_once():
    il, im = match_cadences(np.array([0.0, 5.0]), np.array([0.0, 0.001]), tol=0.01)
    assert len(il) == len(set(il.tolist())) == 1


@pytest.mark.parametrize(
    "t_local,t_mast", [([], [1.0]), ([1.0], []), ([], [])], ids=["no-local", "no-mast", "both"]
)
def test_match_cadences_empty_inputs(t_local, t_mast):
    il, im = match_cadences(np.array(t_local), np.array(t_mast))
    assert il.size == im.size == 0


def test_match_cadences_single_local_cadence():
    il, im = match_cadences(np.array([1.0]), np.array([1.0, 1.5]), tol=0.01)
    assert il.tolist() == [0] and im.tolist() == [0]


# --- series comparison -----------------------------------------------------


def test_compare_series_identical_up_to_scale():
    x = 1 + np.random.default_rng(3).normal(0, 1e-3, 1000)
    m = compare_series(3.0 * x, x, column="aperture_flux")
    assert m.pearson_r == pytest.approx(1.0)
    assert m.median_ratio == pytest.approx(3.0)
    assert m.resid_std == pytest.approx(0.0, abs=1e-12)
    assert m.scatter_ratio == pytest.approx(1.0)


def test_compare_series_independent_noise_is_uncorrelated():
    rng = np.random.default_rng(4)
    a, b = 1 + rng.normal(0, 1e-3, 3000), 1 + rng.normal(0, 1e-3, 3000)
    m = compare_series(a, b)
    assert abs(m.pearson_r) < 0.1
    assert m.resid_to_scatter == pytest.approx(np.sqrt(2), rel=0.1)


def test_compare_series_ignores_nonfinite_cadences():
    a = np.array([1.0, 2.0, np.nan, 4.0, 5.0])
    b = np.array([1.0, 2.0, 3.0, np.inf, 5.0])
    assert compare_series(a, b).n_compared == 3


def test_compare_series_all_nan_returns_nan_metrics():
    m = compare_series(np.ones(10), np.full(10, np.nan), column="psf_flux")
    assert m.n_compared == 0
    assert np.isnan(m.pearson_r) and np.isnan(m.resid_to_scatter)


def test_compare_series_rejects_shape_mismatch():
    with pytest.raises(ValueError, match="shape mismatch"):
        compare_series(np.ones(3), np.ones(4))


# --- full light-curve comparison --------------------------------------------


def test_mast_good_mask_combines_tess_and_tglc_flags():
    lc = lk.LightCurve(time=np.arange(4.0), flux=np.ones(4))
    lc["tess_flags"] = np.array([0, 1, 0, 0], dtype=np.int16)
    lc["tglc_flags"] = np.array([0, 0, 1, 0], dtype=np.int16)
    assert mast_good_mask(lc).tolist() == [True, False, False, True]


def test_compare_light_curves_passes_for_reproduced_curve():
    t, local, mast = _synthetic_pair()
    flags = np.zeros(len(t), dtype=np.int16)
    result = compare_light_curves(_make_lc(t, local), _make_lc(t, mast, flags), "T", 1)
    assert result.passed, result.failures
    assert result.n_matched == len(t)
    assert result.columns["aperture_flux"]["median_ratio"] == pytest.approx(5100 / 5000, rel=1e-3)


def test_compare_light_curves_fails_for_different_star():
    t, local, _ = _synthetic_pair(seed=5)
    _, _, other = _synthetic_pair(seed=6)
    other = {c: np.random.default_rng(7).normal(1, 1e-3, len(t)) for c in other}
    result = compare_light_curves(_make_lc(t, local), _make_lc(t, other), "T", 1)
    assert not result.passed
    assert any("pearson_r" in f for f in result.failures)


def test_compare_light_curves_excludes_flagged_mast_cadences():
    t, local, mast = _synthetic_pair()
    flags = np.zeros(len(t), dtype=np.int16)
    flags[:100] = 128
    result = compare_light_curves(_make_lc(t, local), _make_lc(t, mast, flags))
    assert result.n_mast == len(t) - 100
    assert result.n_matched == len(t) - 100


def test_compare_light_curves_skips_unavailable_psf_columns():
    """MAST PSF columns can be all-NaN; that must not count as a failure."""
    t, local, mast = _synthetic_pair()
    mast["cal_psf_flux"] = np.full(len(t), np.nan)
    mast["psf_flux"] = np.full(len(t), np.nan)
    result = compare_light_curves(_make_lc(t, local), _make_lc(t, mast))
    assert result.passed, result.failures
    assert result.columns["cal_psf_flux"]["n_compared"] == 0


def test_compare_light_curves_fails_when_nothing_comparable():
    t, local, mast = _synthetic_pair()
    for c in COLUMNS:
        mast[c] = np.full(len(t), np.nan)
    result = compare_light_curves(_make_lc(t, local), _make_lc(t, mast))
    assert not result.passed
    assert any("no calibrated column" in f for f in result.failures)


def test_compare_light_curves_fails_on_poor_time_overlap():
    t, local, mast = _synthetic_pair()
    shifted = _make_lc(t + 30.0, mast)
    result = compare_light_curves(_make_lc(t, local), shifted)
    assert result.n_matched == 0
    assert any("overlap" in f for f in result.failures)


def test_benchmark_result_to_dict_is_json_serialisable():
    import json

    t, local, mast = _synthetic_pair()
    result = compare_light_curves(_make_lc(t, local), _make_lc(t, mast), "T", 2)
    json.dumps(result.to_dict(), default=float)


# --- end-to-end benchmark against MAST ----------------------------------------


@pytest.mark.network
@pytest.mark.slow
def test_local_tglc_reproduces_mast_hlsp(tmp_path):
    """Full benchmark: local get_tglc_lc vs MAST TGLC HLSP for TOI-270 s3."""
    from quicklook.benchmark_tglc import run_benchmark

    try:
        result, _, _ = run_benchmark(
            "TIC 259377017", 3, size=99, cache_dir=os.path.expanduser("~/.tglc")
        )
    except (ConnectionError, TimeoutError, LookupError) as e:
        if os.getenv("CI", "false") == "true":
            pytest.xfail(f"Network failure (not a bug): {e}")
        raise
    assert result.passed, result.failures
