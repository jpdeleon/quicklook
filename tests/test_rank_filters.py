"""Tests for the automatic vetting filters in ``quicklook.cli.app._apply_rank_filters``.

These metrics (odd-even mismatch, secondary eclipse significance, per-transit
depth/duration consistency, Gaia RUWE) are computed by ``TessQuickLook`` per
candidate but were previously dropped before reaching this filter. A missing
or NaN metric must never exclude a candidate -- only a metric that is present
and out of bounds should.
"""

import numpy as np
import pandas as pd
import pytest

from quicklook.cli.app import _apply_rank_filters


def _base_row(**overrides):
    row = {
        "TOI": np.nan,
        "depth": 5.0,
        "SDE_tls": 10.0,
        "Prot_gls": 13.5,
        "Porb_tls": 4.0,  # neither ratio's fractional part is near 0/1 (avoids the rotation-harmonic filter)
        "simbad_object": None,
        "exptime": 30,
        "per_transit_count": 15,
        "odd_even_mismatch": 0.5,
        "secondary_sde": 1.0,
        "depth_variance_ratio": 0.1,
        "duration_consistency_ratio": 1.0,
        "gaia_ruwe": 1.0,
    }
    row.update(overrides)
    return row


def _make_df(rows):
    return pd.DataFrame(rows)


def test_clean_candidate_survives_all_filters():
    df = _make_df([_base_row()])
    result = _apply_rank_filters(df)
    assert len(result) == 1


@pytest.mark.parametrize(
    "column,bad_value",
    [
        ("odd_even_mismatch", 5.0),  # > default max_odd_even_sigma=3.0
        ("secondary_sde", 8.0),  # > default max_secondary_sde=5.0
        ("depth_variance_ratio", 0.9),  # > default max_depth_variance_ratio=0.5
        ("gaia_ruwe", 2.0),  # > default max_ruwe=1.4
    ],
)
def test_flagged_metric_excludes_candidate_when_a_clean_one_remains(column, bad_value):
    df = _make_df([_base_row(), _base_row(**{column: bad_value})])
    result = _apply_rank_filters(df)
    assert len(result) == 1
    assert result.iloc[0][column] != bad_value


@pytest.mark.parametrize("bad_ratio", [0.1, 5.0])  # outside default (0.3, 3.0) bounds
def test_duration_ratio_outside_bounds_excludes_candidate(bad_ratio):
    df = _make_df([_base_row(), _base_row(duration_consistency_ratio=bad_ratio)])
    result = _apply_rank_filters(df)
    assert len(result) == 1
    assert result.iloc[0]["duration_consistency_ratio"] == pytest.approx(1.0)


@pytest.mark.parametrize(
    "column",
    [
        "odd_even_mismatch",
        "secondary_sde",
        "depth_variance_ratio",
        "duration_consistency_ratio",
        "gaia_ruwe",
    ],
)
def test_nan_metric_does_not_exclude_candidate(column):
    df = _make_df([_base_row(**{column: np.nan})])
    result = _apply_rank_filters(df)
    assert len(result) == 1


def test_missing_vetting_columns_are_backward_compatible():
    """Older CSVs summarized before this change lack these columns entirely."""
    row = _base_row()
    for column in (
        "odd_even_mismatch",
        "secondary_sde",
        "depth_variance_ratio",
        "duration_consistency_ratio",
        "gaia_ruwe",
    ):
        del row[column]
    df = _make_df([row])
    result = _apply_rank_filters(df)
    assert len(result) == 1


def test_thresholds_are_configurable():
    df = _make_df([_base_row(), _base_row(gaia_ruwe=2.0)])
    result = _apply_rank_filters(df, max_ruwe=3.0)
    assert len(result) == 2
