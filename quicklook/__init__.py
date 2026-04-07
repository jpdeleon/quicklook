"""TESS QuickLook -- lazy top-level imports.

Heavy dependencies (matplotlib, astropy, etc.) are only loaded when the
user accesses a name that lives in one of the submodules.  This keeps CLI
startup fast (~0.1 s instead of ~5 s).

Interactive / notebook usage is unchanged:
    >>> from quicklook import TessQuickLook, Gls, plot
"""

import importlib as _importlib

_SUBMODULES = {"plot", "measure", "gls", "utils", "tql", "inject", "h5io"}


def __getattr__(name: str):
    # Direct submodule access  (e.g. quicklook.plot, quicklook.h5io)
    if name in _SUBMODULES:
        return _importlib.import_module(f".{name}", __name__)

    # Public name from a submodule  (e.g. quicklook.TessQuickLook)
    # Try each submodule until found.
    for submod_name in _SUBMODULES:
        try:
            mod = _importlib.import_module(f".{submod_name}", __name__)
        except ImportError:
            continue
        try:
            return getattr(mod, name)
        except AttributeError:
            continue

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    names = list(_SUBMODULES)
    for submod_name in _SUBMODULES:
        try:
            mod = _importlib.import_module(f".{submod_name}", __name__)
        except ImportError:
            continue
        public = getattr(mod, "__all__", None)
        if public is None:
            public = [n for n in dir(mod) if not n.startswith("_")]
        names.extend(public)
    return names
