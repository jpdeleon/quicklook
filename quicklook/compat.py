"""
Compatibility module for different Python versions.
"""

from pathlib import Path
from importlib.resources import files, as_file


def get_data_path(package_name: str) -> Path:
    """
    Get the path to the data directory for a package.
    """
    try:
        data_path = files(package_name).joinpath("data")
        return Path(as_file(data_path).__enter__())
    except Exception:
        pass

    # Fallback using __file__
    import sys
    import importlib.util

    spec = importlib.util.find_spec(package_name)
    if spec and spec.origin:
        return Path(spec.origin).parent.joinpath("data")

    return Path(sys.modules[package_name].__file__).parent.joinpath("data")
