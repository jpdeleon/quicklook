"""Tests for pure helpers in quicklook.utils."""

import numpy as np
import pytest

from quicklook import utils


# --- mag_to_flux -----------------------------------------------------------


def test_mag_to_flux_zero_mag_returns_unit_flux():
    assert utils.mag_to_flux(np.array([0.0]))[0] == pytest.approx(1.0)


def test_mag_to_flux_known_value():
    """A 5-mag difference is exactly a 100x flux ratio (2.5 dex)."""
    flux = utils.mag_to_flux(np.array([5.0]))
    assert flux[0] == pytest.approx(0.01)


def test_mag_to_flux_propagates_uncertainty():
    flux, ferr = utils.mag_to_flux(np.array([0.0]), np.array([0.1]))
    # df/dm = -0.4 * ln(10) * f  →  |df| ≈ 0.921 * f * dm = 0.0921 at f=1
    assert flux[0] == pytest.approx(1.0)
    assert ferr[0] == pytest.approx(0.0921, abs=1e-4)


def test_mag_to_flux_no_err_returns_single_array():
    """Without mag_err, the function returns just the flux array (not a tuple)."""
    out = utils.mag_to_flux(np.array([1.0, 2.0, 3.0]))
    assert isinstance(out, np.ndarray)
    assert out.shape == (3,)


# --- get_cartersian_distance ----------------------------------------------


def test_distance_zero_for_same_point():
    assert utils.get_cartersian_distance(1, 2, 1, 2) == 0


def test_distance_pythagorean_triple():
    assert utils.get_cartersian_distance(0, 0, 3, 4) == pytest.approx(5.0)


def test_distance_is_symmetric():
    assert utils.get_cartersian_distance(2, 5, -1, 1) == utils.get_cartersian_distance(-1, 1, 2, 5)


# --- is_point_inside_mask --------------------------------------------------


def _square():
    """Closed square border (last point = first).

    Vertices listed clockwise — required by ``is_point_inside_mask``'s
    signed-winding implementation, which expects a +360° total turn for
    interior points.
    """
    return [(0, 0), (0, 10), (10, 10), (10, 0), (0, 0)]


def test_point_inside_square():
    assert utils.is_point_inside_mask(_square(), (5, 5)) is True


def test_point_outside_square():
    assert utils.is_point_inside_mask(_square(), (20, 20)) is False


def test_point_far_outside_square():
    assert utils.is_point_inside_mask(_square(), (-50, 5)) is False


# --- make_round_mask -------------------------------------------------------


def test_round_mask_shape_matches_image():
    img = np.zeros((11, 11))
    img[5, 5] = 10.0  # brightest at centre
    mask = utils.make_round_mask(img, radius=2)
    assert mask.shape == img.shape
    assert mask.dtype == bool


def test_round_mask_contains_centre_pixel():
    img = np.zeros((11, 11))
    img[5, 5] = 10.0
    mask = utils.make_round_mask(img, radius=2)
    assert mask[5, 5] is np.True_ or mask[5, 5] == True  # noqa: E712


def test_round_mask_excludes_far_corner():
    img = np.zeros((11, 11))
    img[5, 5] = 10.0
    mask = utils.make_round_mask(img, radius=2)
    assert not mask[0, 0]
    assert not mask[10, 10]


# --- make_square_mask ------------------------------------------------------


def test_square_mask_shape_and_dtype():
    img = np.zeros((9, 9))
    img[4, 4] = 5.0
    mask = utils.make_square_mask(img, size=1)
    assert mask.shape == img.shape
    assert mask.dtype == bool


def test_square_mask_covers_center_block():
    img = np.zeros((9, 9))
    img[4, 4] = 5.0
    mask = utils.make_square_mask(img, size=1)
    # 3×3 block around (4,4) — rows/cols 3..5 inclusive — should all be True
    assert mask[3:6, 3:6].all()
    # corners should be False
    assert not mask[0, 0]
    assert not mask[8, 8]


# --- get_params_from_exofop ------------------------------------------------


def test_get_params_returns_none_when_section_missing():
    assert utils.get_params_from_exofop({}, name="planet_parameters") is None
    assert utils.get_params_from_exofop({"planet_parameters": []}, name="planet_parameters") is None


def test_get_params_picks_latest_by_pdate_when_idx_none():
    info = {
        "planet_parameters": [
            {"pdate": "2024-01-01", "per": 1.0},
            {"pdate": "2025-06-30", "per": 2.0},  # latest
            {"pdate": "2024-12-15", "per": 3.0},
        ]
    }
    out = utils.get_params_from_exofop(info, name="planet_parameters")
    assert out["per"] == 2.0


def test_get_params_picks_latest_by_sdate_for_stellar():
    info = {
        "stellar_parameters": [
            {"sdate": "2023-05-01", "srad": 1.0},
            {"sdate": "2026-01-01", "srad": 1.2},  # latest
        ]
    }
    out = utils.get_params_from_exofop(info, name="stellar_parameters")
    assert out["srad"] == 1.2


def test_get_params_explicit_idx_overrides_date_sort():
    info = {
        "planet_parameters": [
            {"pdate": "2024-01-01", "per": 1.0},
            {"pdate": "2025-06-30", "per": 2.0},
        ]
    }
    out = utils.get_params_from_exofop(info, name="planet_parameters", idx=0)
    assert out["per"] == 1.0


def test_get_params_falls_back_to_index_0_when_dates_unparseable():
    info = {
        "planet_parameters": [
            {"pdate": "not-a-date", "per": 9.99},
            {"pdate": "also-bad", "per": 0.0},
        ]
    }
    out = utils.get_params_from_exofop(info, name="planet_parameters")
    assert out["per"] == 9.99


def test_get_params_handles_bad_idx_safely():
    info = {"planet_parameters": [{"per": 1.0}]}
    out = utils.get_params_from_exofop(info, name="planet_parameters", idx=99)
    assert out is None  # IndexError caught, logged, None returned
