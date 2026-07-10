"""Tests for target-name sanitization in quicklook.cli.ql.

The sanitized name is interpolated straight into filesystem paths (the
per-target ``.log`` file and the figure/H5 output names), so names that
could escape their directory must be rejected outright.
"""

import pytest

from quicklook.cli.ql import sanitize_target_name
from quicklook.exceptions import InvalidInputError


# --- normalization ---------------------------------------------------------


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("  TOI-1234  ", "TOI-1234"),
        ("TOI 1234", "TOI-1234"),
        ("toi 1234", "toi-1234"),
        ("TOI1234", "TOI-1234"),
    ],
)
def test_toi_names_normalize_to_single_hyphen(raw, expected):
    assert sanitize_target_name(raw) == expected


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("WASP 21", "WASP21"),
        ("  HD 189733  ", "HD189733"),
        ("TIC 123456789", "TIC123456789"),
        ("Gaia DR2 12345", "GaiaDR212345"),
    ],
)
def test_non_toi_names_drop_spaces(raw, expected):
    assert sanitize_target_name(raw) == expected


# --- path traversal rejection ----------------------------------------------


@pytest.mark.parametrize(
    "malicious",
    [
        "../../../../home/jerome/.bashrc",
        "..",
        "../TOI-1234",
        "TOI-1234/../../etc/passwd",
        "/etc/passwd",
        "subdir/TOI-1234",
        "TOI-1234/",
        "back\\slash",
        "nul\x00byte",
        ".hidden",
    ],
)
def test_rejects_names_that_escape_the_output_directory(malicious):
    with pytest.raises(InvalidInputError):
        sanitize_target_name(malicious)


def test_traversal_rejected_even_when_spaces_would_hide_it():
    """Spaces are stripped before validation, so '. . / . .' style input
    must not slip a separator through."""
    with pytest.raises(InvalidInputError):
        sanitize_target_name(". . / . . / etc / passwd")


@pytest.mark.parametrize("blank", ["", "   ", "\t"])
def test_rejects_empty_name(blank):
    with pytest.raises(InvalidInputError):
        sanitize_target_name(blank)


def test_legitimate_name_survives_validation():
    """Guard against an over-eager rejection rule breaking real targets."""
    assert sanitize_target_name("TOI-1234.01") == "TOI-1234.01"
