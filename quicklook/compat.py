"""
Compatibility module for different Python versions.
"""

import sys
from pathlib import Path


def get_data_path(package_name: str) -> Path:
    """
    Get the path to the data directory for a package, compatible across Python versions.
    """
    # Python ≥3.9
    if sys.version_info >= (3, 9):
        try:
            from importlib.resources import files, as_file

            data_path = files(package_name).joinpath("data")
            # as_file ensures we get a real filesystem path
            return Path(as_file(data_path).__enter__())
        except Exception:
            pass

    # Backport for older Python
    try:
        from importlib_resources import files, as_file

        data_path = files(package_name).joinpath("data")
        return Path(as_file(data_path).__enter__())
    except Exception:
        pass

    # Fallback using __file__
    import importlib.util
    import os

    spec = importlib.util.find_spec(package_name)
    if spec and spec.origin:
        return Path(spec.origin).parent.joinpath("data")

    # Last resort
    return Path(sys.modules[package_name].__file__).parent.joinpath("data")
