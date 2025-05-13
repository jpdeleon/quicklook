"""
Compatibility module for different Python versions.
"""

import os
import sys
import importlib
import pathlib
from pathlib import Path


def get_data_path(package_name):
    """
    Get the path to the data directory for a package.

    This is a compatibility function that works across different Python versions.
    It provides similar functionality to importlib.resources.files introduced in Python 3.9.

    Parameters
    ----------
    package_name : str
        The name of the package

    Returns
    -------
    pathlib.Path
        Path object pointing to the package directory
    """
    # Try to use importlib.resources.files (Python 3.9+)
    if sys.version_info >= (3, 9):
        try:
            from importlib.resources import files

            return files(package_name)
        except (ImportError, ModuleNotFoundError):
            pass

    # Try to use importlib_resources backport
    try:
        from importlib_resources import files

        return files(package_name)
    except (ImportError, ModuleNotFoundError):
        pass

    # Fallback for older Python versions
    try:
        # Get the package spec
        spec = importlib.util.find_spec(package_name)
        if spec is not None and spec.origin is not None:
            # If it's a directory (package), return the directory path
            if os.path.isdir(spec.origin):
                return Path(spec.origin)
            # If it's a file (module), return the parent directory
            return Path(spec.origin).parent
    except (ImportError, AttributeError):
        pass

    # Last resort fallback
    package_path = os.path.dirname(sys.modules[package_name].__file__)
    return Path(package_path)
