import pytest
import numpy as np
import os
import sys
from unittest.mock import patch, MagicMock
from matplotlib.figure import Figure
from quicklook.tql import TessQuickLook
import lightkurve as lk
import astropy.units as u
import requests
from astroquery.exceptions import RemoteServiceError

# Network errors that indicate MAST/remote service issues, not code bugs.
_NETWORK_ERRORS = (RemoteServiceError, requests.exceptions.HTTPError, ConnectionError, TimeoutError)


@pytest.fixture
def planet_inputs():
    """Fixture for a known exoplanet target"""
    return {
        "target_name": "WASP-21",
        "sector": 56,
        "flux_type": "pdcsap",
        "pipeline": "SPOC",
        "verbose": False,  # Reduce output during tests
    }


@pytest.fixture
def eb_inputs():
    """Fixture for a known eclipsing binary target"""
    return {
        "target_name": "TIC 144539611",
        "sector": 4,
        "flux_type": "pdcsap",
        "pipeline": "SPOC",
        "verbose": False,  # Reduce output during tests
    }


@pytest.fixture
def variable_star_inputs():
    """Fixture for a known variable star target"""
    return {
        "target_name": "TIC 277539431",  # RR Lyrae variable
        "sector": 12,
        "flux_type": "pdcsap",
        "pipeline": "SPOC",
        "verbose": False,
    }


@pytest.mark.network
def test_tql_planet(planet_inputs):
    """Test TessQuickLook with a known exoplanet target"""
    try:
        ql = TessQuickLook(**planet_inputs)
        assert isinstance(ql, TessQuickLook)
        fig = ql.plot_tql()
        assert isinstance(fig, Figure)

        # Check that key attributes were set correctly
        assert ql.target_name == planet_inputs["target_name"]
        assert ql.sector == planet_inputs["sector"]
        assert ql.flux_type == planet_inputs["flux_type"].lower()
        assert ql.pipeline == planet_inputs["pipeline"].lower()

        # Check that light curves were created
        assert isinstance(ql.raw_lc, lk.LightCurve)
        assert isinstance(ql.flat_lc, lk.LightCurve)
        assert isinstance(ql.trend_lc, lk.LightCurve)

        # Check that TLS was run
        assert hasattr(ql, "tls_results")
        assert ql.tls_results.period > 0
    except _NETWORK_ERRORS as e:
        # CI-specific soft failure — MAST outages are not code bugs
        if os.getenv("CI", "false") == "true":
            pytest.xfail(f"Network failure (not a bug): {e}")
        else:
            raise  # fail locally


@pytest.mark.network
def test_tql_eb(eb_inputs):
    """Test TessQuickLook with a known eclipsing binary target"""
    try:
        ql = TessQuickLook(**eb_inputs)
        assert isinstance(ql, TessQuickLook)
        fig = ql.plot_tql()
        assert isinstance(fig, Figure)

        # Check that key attributes were set correctly
        assert ql.target_name == eb_inputs["target_name"]
        assert ql.sector == eb_inputs["sector"]
        assert ql.flux_type == eb_inputs["flux_type"].lower()
        assert ql.pipeline == eb_inputs["pipeline"].lower()

        # Check that light curves were created
        assert isinstance(ql.raw_lc, lk.LightCurve)
        assert isinstance(ql.flat_lc, lk.LightCurve)
        assert isinstance(ql.trend_lc, lk.LightCurve)

        # Check that TLS was run
        assert hasattr(ql, "tls_results")
        assert ql.tls_results.period > 0
    except _NETWORK_ERRORS as e:
        # CI-specific soft failure — MAST outages are not code bugs
        if os.getenv("CI", "false") == "true":
            pytest.xfail(f"Network failure (not a bug): {e}")
        else:
            raise  # fail locally


@pytest.mark.network
def test_tql_variable_star(variable_star_inputs):
    """Test TessQuickLook with a known variable star target"""
    try:
        ql = TessQuickLook(**variable_star_inputs)
        assert isinstance(ql, TessQuickLook)
        fig = ql.plot_tql()
        assert isinstance(fig, Figure)

        # Check that key attributes were set correctly
        assert ql.target_name == variable_star_inputs["target_name"]
        assert ql.sector == variable_star_inputs["sector"]
        assert ql.flux_type == variable_star_inputs["flux_type"].lower()
        assert ql.pipeline == variable_star_inputs["pipeline"].lower()
    except _NETWORK_ERRORS as e:
        # CI-specific soft failure — MAST outages are not code bugs
        if os.getenv("CI", "false") == "true":
            pytest.xfail(f"Network failure (not a bug): {e}")
        else:
            raise  # fail locally


@pytest.mark.skipif(os.getenv("CI") == "true", reason="Skip performance tests in CI")
@pytest.mark.benchmark
def test_tql_runtime(benchmark, planet_inputs):
    def run_ql():
        ql = TessQuickLook(**planet_inputs)
        _ = ql.plot_tql()

    benchmark(run_ql)


def test_with_mock_light_curve(mock_light_curve, planet_inputs):
    """Test TessQuickLook with a mock light curve and no external services."""
    inputs = planet_inputs.copy()
    inputs["tls_use_threads"] = 1
    exofop_data = {
        "basic_info": {
            "star_names": "WASP-21, TIC 234523599, Gaia DR3 2345395591154370176",
            "tic_id": "234523599",
        },
        "coordinates": {"ra": 349.8985, "dec": -10.3270},
        "planet_parameters": [],
    }

    with (
        patch("quicklook.tql.get_exofop_json", return_value=exofop_data),
        patch.object(TessQuickLook, "get_simbad_obj_type", return_value=None),
        patch.object(TessQuickLook, "get_lc", return_value=mock_light_curve),
        patch.object(TessQuickLook, "check_output_file_exists", return_value=None),
    ):
        ql = TessQuickLook(**inputs)

    # Manually set the required attributes
    ql.pipeline = inputs["pipeline"].lower()
    ql.sector = inputs["sector"]

    # Check that the light curve was set correctly
    assert ql.raw_lc is mock_light_curve

    # Check that the flattened light curve was created
    assert isinstance(ql.flat_lc, lk.LightCurve)
    assert isinstance(ql.trend_lc, lk.LightCurve)

    # Check that TLS was run
    assert hasattr(ql, "tls_results")
    assert hasattr(ql.tls_results, "period")

    # Test basic attributes
    assert ql.target_name == inputs["target_name"]
    assert ql.flux_type == inputs["flux_type"].lower()
    assert ql.pipeline == inputs["pipeline"].lower()
    assert ql.sector == inputs["sector"]


@pytest.mark.parametrize("flatten_method", ["notch", "locor"])
def test_with_mock_light_curve_notch_locor(mock_light_curve, planet_inputs, flatten_method):
    """flatten_raw_lc() dispatches to quicklook.notch_locor for these methods
    instead of wotan, including when the GLS-based auto period/window falls
    back gracefully (get_lc is mocked here, so self.pipeline isn't set yet
    when flatten_raw_lc() runs -- init_gls() should fail soft, not crash)."""
    inputs = planet_inputs.copy()
    inputs["tls_use_threads"] = 1
    inputs["flatten_method"] = flatten_method
    exofop_data = {
        "basic_info": {
            "star_names": "WASP-21, TIC 234523599, Gaia DR3 2345395591154370176",
            "tic_id": "234523599",
        },
        "coordinates": {"ra": 349.8985, "dec": -10.3270},
        "planet_parameters": [],
    }

    with (
        patch("quicklook.tql.get_exofop_json", return_value=exofop_data),
        patch.object(TessQuickLook, "get_simbad_obj_type", return_value=None),
        patch.object(TessQuickLook, "get_lc", return_value=mock_light_curve),
        patch.object(TessQuickLook, "check_output_file_exists", return_value=None),
    ):
        ql = TessQuickLook(**inputs)

    assert isinstance(ql.flat_lc, lk.LightCurve)
    assert isinstance(ql.trend_lc, lk.LightCurve)
    assert np.isfinite(ql.flat_lc.flux.value).all()
    assert ql.window_length > 0


def _bare_ql(raw_lc, toi_epoch, toi_period, toi_dur, sector=1, all_sectors=None):
    """Build a TessQuickLook without running __init__ (no network).

    Only the attributes touched by ``get_transit_mask`` and
    ``_warn_no_transits_in_sector`` are populated.
    """
    ql = TessQuickLook.__new__(TessQuickLook)
    ql.raw_lc = raw_lc
    ql.toi_epoch = toi_epoch
    ql.toi_period = toi_period
    ql.toi_dur = toi_dur
    ql.sector = sector
    ql.all_sectors = all_sectors if all_sectors is not None else [sector]
    ql.ephem_source = "TFOP"
    return ql


def test_get_transit_mask_empty_when_transit_outside_baseline(mock_light_curve):
    """A long-period ephemeris whose transit misses the sector yields no mask.

    Mirrors TOI-2074 (P=177.6 d): the transit does not fall in a single
    ~27-day sector, so the mask must be all-False rather than raising.
    """
    ql = _bare_ql(
        mock_light_curve,
        toi_epoch=np.array((1000.0, 0.01)),  # far outside the 0-27 d baseline
        toi_period=np.array((177.58, 0.005)),
        toi_dur=np.array((0.2, 0.01)),
    )

    tmask = ql.get_transit_mask()

    assert tmask.dtype == bool
    assert len(tmask) == len(mock_light_curve.time)
    assert tmask.sum() == 0


def test_get_transit_mask_flags_transit_within_baseline(mock_light_curve):
    """A short-period ephemeris inside the baseline flags in-transit cadences."""
    ql = _bare_ql(
        mock_light_curve,
        toi_epoch=np.array((5.0, 0.01)),  # within the 0-27 d baseline
        toi_period=np.array((3.0, 0.005)),
        toi_dur=np.array((0.2, 0.01)),
    )

    tmask = ql.get_transit_mask()

    assert tmask.sum() > 0


def test_warn_no_transits_in_sector_does_not_raise(mock_light_curve):
    """Zero predicted transits warn (with guidance) instead of raising."""
    ql = _bare_ql(
        mock_light_curve,
        toi_epoch=np.array((1000.0, 0.01)),
        toi_period=np.array((177.58, 0.005)),
        toi_dur=np.array((0.2, 0.01)),
        sector=77,
        all_sectors=[16, 23, 50, 77],
    )

    with patch("quicklook.tql.logger.warning") as mock_warn:
        ql._warn_no_transits_in_sector()  # must not raise

    mock_warn.assert_called_once()
    msg = mock_warn.call_args[0][0]
    assert "sector 77" in msg
    assert "Nearest predicted transit" in msg
    assert "[16, 23, 50, 77]" in msg


def test_format_available_sectors_single_sector():
    assert TessQuickLook._format_available_sectors([56]) == "56"


def test_format_available_sectors_wraps_short_lists():
    summary = TessQuickLook._format_available_sectors([14, 15, 19, 22, 26, 40])

    assert summary == "14, 15, 19, 22\n26, 40"


def test_format_available_sectors_counts_long_lists():
    summary = TessQuickLook._format_available_sectors(list(range(1, 12)))

    assert summary == "11 sectors available"


def test_format_available_sectors_lists_at_max_listed_boundary():
    summary = TessQuickLook._format_available_sectors(list(range(1, 9)))

    assert summary == "1, 2, 3, 4\n5, 6, 7, 8"


def test_format_available_sectors_counts_one_past_max_listed():
    summary = TessQuickLook._format_available_sectors(list(range(1, 10)))

    assert summary == "9 sectors available"


# --- Gaia RUWE extraction for the summary panel ----------------------------


def _gaia_table(rows):
    import pandas as pd

    return pd.DataFrame(rows)


def test_extract_ruwe_returns_target_value_matched_on_source_id():
    sources = _gaia_table(
        [
            {"source_id": 111, "ruwe": 2.31},
            {"source_id": 222, "ruwe": 0.98},
        ]
    )
    assert TessQuickLook._extract_ruwe(sources, 222) == pytest.approx(0.98)


def test_extract_ruwe_nan_when_no_row_matches():
    sources = _gaia_table([{"source_id": 111, "ruwe": 1.05}])
    assert np.isnan(TessQuickLook._extract_ruwe(sources, 999))


def test_extract_ruwe_nan_when_ruwe_column_absent():
    sources = _gaia_table([{"source_id": 111, "phot_g_mean_mag": 12.0}])
    assert np.isnan(TessQuickLook._extract_ruwe(sources, 111))


def test_extract_ruwe_nan_when_sources_is_none():
    assert np.isnan(TessQuickLook._extract_ruwe(None, 111))


def test_extract_ruwe_nan_when_gaiaid_not_coercible():
    sources = _gaia_table([{"source_id": 111, "ruwe": 1.4}])
    assert np.isnan(TessQuickLook._extract_ruwe(sources, "not-an-id"))


# --- summary panel: RUWE + Num. Transits rows ------------------------------
#
# _summary_sections needs a lot of object state but no network, so build a bare
# TessQuickLook and stub the ExoFOP stellar-parameter parse.

from quicklook.tql import _TLSResult  # noqa: E402

_STAR_PARAMS = {
    "srad": 1.1, "srad_e": 0.05, "mass": 1.0, "mass_e": 0.05,
    "teff": 5800, "teff_e": 100, "logg": 4.4, "logg_e": 0.1,
    "dist": 120.0, "dist_e": 5.0,
}  # fmt: skip


def _bare_summary_ql(distinct_transit_count=3, gaia_ruwe=1.42, n_planets=1):
    from astropy.coordinates import SkyCoord

    ql = TessQuickLook.__new__(TessQuickLook)
    ql.pipeline = "qlp"
    ql.ephem_source = "TFOP"
    ql.sector = 56
    ql.all_sectors = [56]
    ql.gaiaid = 12345
    ql.gaia_ruwe = gaia_ruwe
    ql.Prot_ls = 3.2
    ql.nearby_star_sep = None
    ql.simbad_obj_type = None
    ql.target_coord = SkyCoord(ra=350.0, dec=-10.0, unit="deg")
    ql.toi_period = ql.toi_epoch = ql.toi_dur = ql.toi_depth = ql.toi_rp = None
    lc = lk.LightCurve(time=[1, 2, 3], flux=[1.0, 1.0, 1.0])
    lc.meta = {"FLUX_ORIGIN": "pdcsap"}
    ql.raw_lc = lc
    tls_values = {
        "SDE": 21.4, "period": 8.213, "period_uncertainty": 0.002,
        "T0": 1650.1, "duration": 0.12, "depth": 0.997, "rp_rs": 0.05,
        "odd_even_mismatch": 0.3,
    }  # fmt: skip
    if distinct_transit_count is not None:
        tls_values["distinct_transit_count"] = distinct_transit_count
    ql.tls_results = _TLSResult(tls_values)
    ql.iterative_tls_results = [ql.tls_results] * n_planets
    ql.exofop_data = {"magnitudes": [{"band": "T", "value": 11.2}]}
    return ql


def test_summary_adds_transits_after_odd_even_and_ruwe_in_stellar():
    with patch("quicklook.tql.get_params_from_exofop", return_value=_STAR_PARAMS):
        sections = _bare_summary_ql(distinct_transit_count=3, gaia_ruwe=1.42)._summary_sections()

    candidate = sections[0][1]
    stellar = sections[1][1]
    cand_labels = [label for label, _ in candidate]
    # Num. Transits sits immediately after Odd-Even in Candidate Properties.
    assert cand_labels[cand_labels.index("Odd-Even") + 1] == "Num. Transits"
    assert dict(candidate)["Num. Transits"] == "3"
    assert "Sector" not in cand_labels
    assert dict(candidate)["Available sectors"] == "56"
    # RUWE appears in the Stellar Properties column.
    assert dict(stellar)["RUWE"] == "1.42"


def test_summary_omits_transits_and_ruwe_when_unavailable():
    """GPU TLS results lack distinct_transit_count and RUWE may be missing —
    both rows must be skipped rather than rendered as errors/NaN."""
    with patch("quicklook.tql.get_params_from_exofop", return_value=_STAR_PARAMS):
        sections = _bare_summary_ql(
            distinct_transit_count=None, gaia_ruwe=np.nan
        )._summary_sections()

    assert "Num. Transits" not in dict(sections[0][1])
    assert "RUWE" not in dict(sections[1][1])


def test_summary_adds_num_planets_row_when_multiple():
    """A multi-planet iterative search surfaces its planet count."""
    with patch("quicklook.tql.get_params_from_exofop", return_value=_STAR_PARAMS):
        sections = _bare_summary_ql(n_planets=3)._summary_sections()

    assert dict(sections[0][1])["Num. Planets"] == "3"


def test_summary_omits_num_planets_row_by_default():
    """Single-planet runs must not show a planet-count row."""
    with patch("quicklook.tql.get_params_from_exofop", return_value=_STAR_PARAMS):
        sections = _bare_summary_ql()._summary_sections()

    assert "Num. Planets" not in dict(sections[0][1])


# @pytest.mark.parametrize("pg_method", ["gls", "lombscargle", "bls"])
# def test_different_periodogram_methods(planet_inputs, pg_method):
#     """Test TessQuickLook with different periodogram methods"""
#     inputs = planet_inputs.copy()
#     inputs["pg_method"] = pg_method

#     ql = TessQuickLook(**inputs)
#     # fig = ql.plot_tql()
#     # assert isinstance(fig, Figure)

#     # Check that the periodogram method was set correctly
#     assert ql.pg_method == pg_method


# @pytest.mark.parametrize(
#     "flatten_method", ["biweight", "median", "mean", "cosine"]
# )
# def test_different_flatten_methods(planet_inputs, flatten_method):
#     """Test TessQuickLook with different flatten methods"""
#     inputs = planet_inputs.copy()
#     inputs["flatten_method"] = flatten_method

#     ql = TessQuickLook(**inputs)
#     # fig = ql.plot_tql()
#     # assert isinstance(fig, Figure)

#     # Check that the flatten method was set correctly
#     assert ql.flatten_method == flatten_method


# def test_custom_period_limits(planet_inputs):
#     """Test TessQuickLook with custom period limits"""
#     inputs = planet_inputs.copy()
#     inputs["Porb_limits"] = (1.0, 5.0)  # Search between 1 and 5 days

#     ql = TessQuickLook(**inputs)
#     fig = ql.plot_tql()
#     assert isinstance(fig, Figure)

#     # Check that the period limits were set correctly
#     assert ql.Porb_min == 1.0
#     assert ql.Porb_max == 5.0

#     # Check that TLS results are within the specified range
#     assert 1.0 <= ql.tls_results.period <= 5.0


# def test_get_transit_mask(planet_inputs):
#     """Test the get_transit_mask method"""
#     ql = TessQuickLook(**planet_inputs)

#     # Get the transit mask
#     tmask = ql.get_transit_mask()

#     # Check that the mask is a boolean array
#     assert isinstance(tmask, np.ndarray)
#     assert tmask.dtype == bool

#     # Check that the mask has the same length as the light curve
#     assert len(tmask) == len(ql.raw_lc.time)


# def test_flatten_raw_lc(planet_inputs):
#     """Test the flatten_raw_lc method"""
#     ql = TessQuickLook(**planet_inputs)

#     # Get the flattened light curve
#     flat_lc, trend_lc = ql.flatten_raw_lc()

#     # Check that the returned objects are light curves
#     assert isinstance(flat_lc, lk.LightCurve)
#     assert isinstance(trend_lc, lk.LightCurve)

#     # Check that the flattened light curve has the same length as the raw light curve
#     assert len(flat_lc.time) == len(ql.raw_lc.time)

#     # Check that the trend light curve has the same length as the raw light curve
#     assert len(trend_lc.time) == len(ql.raw_lc.time)


# @pytest.mark.parametrize("window_length", [0.1, 0.5, 1.0])
# def test_different_window_lengths(planet_inputs, window_length):
#     """Test TessQuickLook with different window lengths for flattening"""
#     inputs = planet_inputs.copy()
#     inputs["window_length"] = window_length

#     ql = TessQuickLook(**inputs)
#     # fig = ql.plot_tql()
#     # assert isinstance(fig, Figure)

#     # Check that the window length was set correctly
#     assert ql.window_length == window_length


# @pytest.mark.parametrize("edge_cutoff", [0.0, 0.1, 0.2])
# def test_different_edge_cutoffs(planet_inputs, edge_cutoff):
#     """Test TessQuickLook with different edge cutoffs for flattening"""
#     inputs = planet_inputs.copy()
#     inputs["edge_cutoff"] = edge_cutoff

#     ql = TessQuickLook(**inputs)
#     # fig = ql.plot_tql()
#     # assert isinstance(fig, Figure)

#     # Check that the edge cutoff was set correctly
#     assert ql.edge_cutoff == edge_cutoff


# @pytest.mark.parametrize("gp_kernel", ["matern", "squared_exp", "periodic"])
# def test_different_gp_kernels(planet_inputs, gp_kernel):
#     """Test TessQuickLook with different GP kernels for flattening"""
#     inputs = planet_inputs.copy()
#     inputs["gp_kernel"] = gp_kernel
#     inputs["flatten_method"] = "gp"  # Use GP flattening

#     ql = TessQuickLook(**inputs)
#     # fig = ql.plot_tql()
#     # assert isinstance(fig, Figure)

#     # Check that the GP kernel was set correctly
#     assert ql.gp_kernel == gp_kernel


# @pytest.mark.parametrize("gp_kernel_size", [0.1, 0.5, 1.0])
# def test_different_gp_kernel_sizes(planet_inputs, gp_kernel_size):
#     """Test TessQuickLook with different GP kernel sizes for flattening"""
#     inputs = planet_inputs.copy()
#     inputs["gp_kernel_size"] = gp_kernel_size
#     inputs["flatten_method"] = "gp"  # Use GP flattening

#     ql = TessQuickLook(**inputs)
#     # fig = ql.plot_tql()
#     # assert isinstance(fig, Figure)

#     # Check that the GP kernel size was set correctly
#     assert ql.gp_kernel_size == gp_kernel_size


# @pytest.mark.parametrize("sigma_clip_raw", [(3, 3), (5, 5), None])
# def test_different_sigma_clip_raw(planet_inputs, sigma_clip_raw):
#     """Test TessQuickLook with different sigma clipping for raw light curve"""
#     inputs = planet_inputs.copy()
#     inputs["sigma_clip_raw"] = sigma_clip_raw

#     ql = TessQuickLook(**inputs)
#     # fig = ql.plot_tql()
#     # assert isinstance(fig, Figure)

#     # Check that the sigma clip raw was set correctly
#     assert ql.sigma_clip_raw == sigma_clip_raw


# @pytest.mark.parametrize("sigma_clip_flat", [(3, 3), (5, 5), None])
# def test_different_sigma_clip_flat(planet_inputs, sigma_clip_flat):
#     """Test TessQuickLook with different sigma clipping for flattened light curve"""
#     inputs = planet_inputs.copy()
#     inputs["sigma_clip_flat"] = sigma_clip_flat

#     ql = TessQuickLook(**inputs)
#     # fig = ql.plot_tql()
#     # assert isinstance(fig, Figure)

#     # Check that the sigma clip flat was set correctly
#     assert ql.sigma_clip_flat == sigma_clip_flat


# def test_custom_ephem(planet_inputs):
#     """Test TessQuickLook with custom ephemeris"""
#     inputs = planet_inputs.copy()
#     # Custom ephemeris: [epoch, period, duration, depth, target_name, source]
#     custom_ephem = [2459000.0, 3.0, 0.1, 0.01, "WASP-21", "custom"]
#     inputs["custom_ephem"] = custom_ephem

#     ql = TessQuickLook(**inputs)
#     fig = ql.plot_tql()
#     assert isinstance(fig, Figure)

#     # Check that the custom ephemeris was set correctly
#     assert ql.custom_ephem == custom_ephem


# def test_mask_ephem(planet_inputs):
#     """Test TessQuickLook with transit masking"""
#     inputs = planet_inputs.copy()
#     inputs["mask_ephem"] = True

#     ql = TessQuickLook(**inputs)
#     # fig = ql.plot_tql()
#     # assert isinstance(fig, Figure)

#     # Check that mask_ephem was set correctly
#     assert ql.mask_ephem is True


# @patch("lightkurve.search_lightcurve")
# def test_get_lc_mocked(mock_search, planet_inputs):
#     """Test get_lc method with mocked lightkurve search"""
#     # Create a mock search result
#     mock_search_result = MagicMock()
#     mock_search_result.table = MagicMock()
#     mock_search_result.table.to_pandas.return_value = MagicMock()

#     # Create a mock light curve
#     mock_lc = MagicMock(spec=lk.LightCurve)
#     mock_lc.meta = {"SECTOR": 56, "EXPOSURE": 120}
#     mock_lc.time = MagicMock()
#     mock_lc.time.jd = np.linspace(2459000, 2459030, 1000)

#     # Set up the mock search to return our mock search result
#     mock_search.return_value = mock_search_result
#     mock_search_result.download = MagicMock(return_value=mock_lc)

#     # Create a TessQuickLook instance with our mocked search
#     inputs = planet_inputs.copy()
#     with patch.object(TessQuickLook, "get_lc", return_value=mock_lc):
#         with patch.object(
#             TessQuickLook, "flatten_raw_lc", return_value=(mock_lc, mock_lc)
#         ):
#             with patch.object(TessQuickLook, "run_tls"):
#                 with patch.object(
#                     TessQuickLook,
#                     "get_transit_mask",
#                     return_value=np.zeros(1000, dtype=bool),
#                 ):
#                     ql = TessQuickLook(**inputs)

#                     # Check that the light curve was set correctly
#                     assert ql.raw_lc is mock_lc


# @patch("lightkurve.search_targetpixelfile")
# def test_get_tpf_mocked(mock_search, planet_inputs):
#     """Test get_tpf method with mocked lightkurve search"""
#     # Create a mock search result
#     mock_search_result = MagicMock()
#     mock_search_result.author = ["SPOC"]
#     mock_search_result.exptime = [120 * u.s]

#     # Create a mock TPF
#     mock_tpf = MagicMock()
#     mock_tpf.meta = {"SECTOR": 56}

#     # Set up the mock search to return our mock search result
#     mock_search.return_value = mock_search_result
#     mock_search_result.download = MagicMock(return_value=mock_tpf)

#     # Create a TessQuickLook instance
#     inputs = planet_inputs.copy()
#     ql = TessQuickLook(**inputs)

#     # Replace the get_tpf method with our mocked version
#     with patch.object(TessQuickLook, "get_tpf", return_value=mock_tpf):
#         tpf = ql.get_tpf(sector=56, author="SPOC")

#         # Check that the TPF was returned correctly
#         assert tpf is mock_tpf


# def test_init_gls(planet_inputs):
#     """Test init_gls method"""
#     ql = TessQuickLook(**planet_inputs)

#     # Initialize GLS
#     gls = ql.init_gls()

#     # Check that GLS was initialized correctly
#     assert hasattr(gls, "power")
#     assert hasattr(gls, "best")


# def test_run_tls(planet_inputs):
#     """Test run_tls method"""
#     ql = TessQuickLook(**planet_inputs)

#     # Run TLS
#     ql.run_tls()

#     # Check that TLS results were set correctly
#     assert hasattr(ql, "tls_results")
#     assert hasattr(ql.tls_results, "period")
#     assert hasattr(ql.tls_results, "power")
#     assert hasattr(ql.tls_results, "SDE")


# def test_append_tls_results(planet_inputs):
#     """Test append_tls_results method"""
#     ql = TessQuickLook(**planet_inputs)

#     # Append TLS results
#     ql.append_tls_results()

#     # Check that TLS results were appended correctly
#     assert hasattr(ql.tls_results, "simbad_obj")


# @pytest.mark.parametrize("archival_survey", ["dss1", "poss2ukstu_red"])
# def test_different_archival_surveys(planet_inputs, archival_survey):
#     """Test TessQuickLook with different archival surveys"""
#     inputs = planet_inputs.copy()
#     inputs["archival_survey"] = archival_survey

#     ql = TessQuickLook(**inputs)

#     # Check that the archival survey was set correctly
#     assert ql.archival_survey == archival_survey


# @patch("lightkurve.search_lightcurve")
# def test_error_handling_no_lightcurve(mock_search, planet_inputs):
#     """Test error handling when no light curve is found"""
#     # Create an empty search result
#     mock_search_result = MagicMock()
#     mock_search_result.__len__.return_value = 0

#     # Set up the mock search to return our empty search result
#     mock_search.return_value = mock_search_result

#     # Create a TessQuickLook instance with our mocked search
#     inputs = planet_inputs.copy()

#     # Check that the correct error is raised
#     with pytest.raises(SystemExit):
#         with patch.object(sys, "exit") as mock_exit:
#             _ = TessQuickLook(**inputs)
#             mock_exit.assert_called_once()


# @patch("lightkurve.search_targetpixelfile")
# def test_error_handling_no_tpf(mock_search, planet_inputs):
#     """Test error handling when no TPF is found"""
#     # Create an empty search result
#     mock_search_result = MagicMock()
#     mock_search_result.__len__.return_value = 0

#     # Set up the mock search to return our empty search result
#     mock_search.return_value = mock_search_result

#     # Create a TessQuickLook instance
#     inputs = planet_inputs.copy()
#     ql = TessQuickLook(**inputs)

#     # Check that the correct error is raised
#     with pytest.raises(SystemExit):
#         with patch.object(sys, "exit") as mock_exit:
#             ql.get_tpf(sector=56, author="SPOC")
#             mock_exit.assert_called_once()


# def test_save_output_files(planet_inputs, tmp_path):
#     """Test saving output files"""
#     # Create a temporary directory for output files
#     outdir = tmp_path / "output"
#     outdir.mkdir()

#     # Create inputs with save options enabled
#     inputs = planet_inputs.copy()
#     inputs["savefig"] = True
#     inputs["savetls"] = True
#     inputs["outdir"] = str(outdir)

#     # Create a TessQuickLook instance
#     ql = TessQuickLook(**inputs)

#     # Plot and save
#     ql.plot_tql()

#     # Check that files were created
#     _ = (
#         outdir
#         / f"{ql.target_name}_s{ql.sector}_{ql.flux_type}_{ql.exptime_suffix}.png"
#     )
#     _ = (
#         outdir
#         / f"{ql.target_name}_s{ql.sector}_{ql.flux_type}_{ql.exptime_suffix}_tls.h5"
#     )

#     # Note: In a real test, we would check if these files exist
#     # But since we're running in a test environment, we'll just check the attributes
#     assert ql.savefig is True
#     assert ql.savetls is True
#     assert ql.outdir == str(outdir)


# @pytest.mark.parametrize(
#     "pipeline", ["SPOC", "QLP", "TESS-SPOC", "CDIPS", "PATHOS"]
# )
# def test_different_pipelines(planet_inputs, pipeline):
#     """Test TessQuickLook with different pipelines"""
#     inputs = planet_inputs.copy()
#     inputs["pipeline"] = pipeline

#     # We need to mock the search and download since different pipelines
#     # might not be available for our test target
#     with patch.object(TessQuickLook, "get_lc") as mock_get_lc:
#         # Create a mock light curve
#         mock_lc = MagicMock(spec=lk.LightCurve)
#         mock_lc.meta = {"SECTOR": 56, "EXPOSURE": 120}
#         mock_lc.time = MagicMock()
#         mock_lc.time.jd = np.linspace(2459000, 2459030, 1000)
#         mock_lc.flux = MagicMock()
#         mock_lc.flux_err = MagicMock()

#         # Set up the mock to return our mock light curve
#         mock_get_lc.return_value = mock_lc

#         # Create a TessQuickLook instance with our mocked get_lc
#         with patch.object(
#             TessQuickLook, "flatten_raw_lc", return_value=(mock_lc, mock_lc)
#         ):
#             with patch.object(TessQuickLook, "run_tls"):
#                 with patch.object(
#                     TessQuickLook,
#                     "get_transit_mask",
#                     return_value=np.zeros(1000, dtype=bool),
#                 ):
#                     ql = TessQuickLook(**inputs)

#                     # Check that the pipeline was set correctly
#                     assert ql.pipeline == pipeline.lower()


# def test_repr_method(planet_inputs):
#     """Test the __repr__ method"""
#     ql = TessQuickLook(**planet_inputs)

#     # Get the string representation
#     repr_str = repr(ql)

#     # Check that the string contains key information
#     assert ql.target_name in repr_str
#     assert str(ql.sector) in repr_str
#     assert ql.pipeline in repr_str
