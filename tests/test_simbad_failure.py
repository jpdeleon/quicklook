"""SIMBAD is optional metadata: an unreachable service must not abort a run."""

from unittest.mock import patch

import pytest

from quicklook.tql import TessQuickLook


def _bare_ql():
    ql = object.__new__(TessQuickLook)
    ql.star_names = ["TOI-137", "TIC 1234"]
    ql.target_name = "TOI-137"
    ql.verbose = False
    return ql


@pytest.mark.parametrize("failing", ["add_votable_fields", "query_object"])
def test_query_simbad_returns_none_when_service_unavailable(failing):
    # Arrange
    ql = _bare_ql()
    with patch("quicklook.tql.Simbad") as simbad:
        getattr(simbad, failing).side_effect = RuntimeError("Connection failed")

        # Act
        result = ql.query_simbad()

    # Assert
    assert result is None


def test_get_simbad_obj_type_is_none_when_service_unavailable():
    # Arrange
    ql = _bare_ql()
    with patch("quicklook.tql.Simbad") as simbad:
        simbad.add_votable_fields.side_effect = RuntimeError("Connection failed")

        # Act
        result = ql.get_simbad_obj_type()

    # Assert
    assert result is None


def test_query_simbad_returns_first_resolved_result():
    # Arrange
    ql = _bare_ql()
    with patch("quicklook.tql.Simbad") as simbad:
        simbad.query_object.side_effect = [None, "table"]

        # Act
        result = ql.query_simbad()

    # Assert
    assert result == "table"
