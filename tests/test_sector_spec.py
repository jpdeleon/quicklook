"""Sector choices shared by the CLI and the web GUI."""

import pytest

from quicklook.sector_spec import (
    is_multi_sector,
    parse_sector_ranges,
    parse_sector_spec,
)


@pytest.mark.parametrize(
    "value, expected",
    [
        (None, -1),
        ("", -1),
        ("  ", -1),
        ("-1", -1),
        (-1, -1),
        ("57", 57),
        (57, 57),
        ("ALL", "all"),
        (" all ", "all"),
        ("57,58", [57, 58]),
        ("57, 58, 77-79", [57, 58, 77, 78, 79]),
        ("84,57,57", [57, 84]),
        ("57-57", 57),
        ([57, 58], [57, 58]),
    ],
)
def test_parse_sector_spec(value, expected):
    assert parse_sector_spec(value) == expected


@pytest.mark.parametrize("value", ["abc", "57;58", "0", "-2", -5, "60-57", "57,", "all,57", 0])
def test_parse_sector_spec_rejects_bad_input(value):
    with pytest.raises(ValueError):
        parse_sector_spec(value)


@pytest.mark.parametrize("value", ["1-1025", "1-1024,2000", "1-10000000000"])
def test_oversized_sector_selection(value):
    with pytest.raises(ValueError, match="at most"):
        parse_sector_ranges(value)


def test_overlapping_ranges_at_limit():
    assert parse_sector_ranges("1-1024,100-200") == list(range(1, 1025))


def test_parse_sector_ranges_matches_cli_format():
    assert parse_sector_ranges("66-68,93") == [66, 67, 68, 93]


@pytest.mark.parametrize(
    "spec, expected", [(-1, False), (57, False), ("all", True), ([1, 2], True)]
)
def test_is_multi_sector(spec, expected):
    assert is_multi_sector(spec) is expected
