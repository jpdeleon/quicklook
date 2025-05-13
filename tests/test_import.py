"""
Test that the package can be imported correctly.
This is a simple test to ensure the package is installed properly.
"""

import pytest


def test_import_quicklook():
    """Test that the quicklook package can be imported."""
    try:
        import quicklook

        assert quicklook is not None
    except ImportError as e:
        pytest.fail(f"Failed to import quicklook: {e}")


def test_import_tql():
    """Test that the TessQuickLook class can be imported."""
    try:
        from quicklook import TessQuickLook

        assert TessQuickLook is not None
    except ImportError as e:
        pytest.fail(f"Failed to import TessQuickLook: {e}")


def test_import_gls():
    """Test that the Gls class can be imported."""
    try:
        from quicklook import Gls

        assert Gls is not None
    except ImportError as e:
        pytest.fail(f"Failed to import Gls: {e}")


def test_import_utils():
    """Test that the utils module can be imported."""
    try:
        from quicklook import utils

        assert utils is not None
    except ImportError as e:
        pytest.fail(f"Failed to import utils: {e}")


def test_import_plot():
    """Test that the plot module can be imported."""
    try:
        from quicklook import plot

        assert plot is not None
    except ImportError as e:
        pytest.fail(f"Failed to import plot: {e}")
