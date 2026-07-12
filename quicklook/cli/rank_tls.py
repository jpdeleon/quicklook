#!/usr/bin/env python
"""CLI redirect to the unified Typer app (``ql rank-tls``)."""

import sys as _sys


def main():
    """Redirect to the unified Typer CLI (``ql rank-tls``)."""
    _sys.argv = [_sys.argv[0], "rank-tls"] + _sys.argv[1:]
    from quicklook.cli.app import app

    app()


if __name__ == "__main__":
    main()
