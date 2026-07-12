#!/usr/bin/env python
"""Target-name sanitization and CLI redirect to the unified Typer app."""

import sys as _sys

from quicklook.exceptions import InvalidInputError


def sanitize_target_name(name: str) -> str:
    """
    Sanitize the user input target name for TessQuickLook.

    Rules:
    - Strip leading/trailing spaces
    - If the target starts with 'TOI', ensure it is formatted as 'TOI-xxxx'
    - Otherwise, remove spaces

    Raises
    ------
    InvalidInputError
        If the sanitized name is empty or could escape the directory it is
        interpolated into.
    """
    name = name.strip()
    if name[:3].lower() == "toi":
        name = name.replace(" ", "-")
        if not name.upper().startswith("TOI-"):
            name = "TOI-" + name[3:].lstrip("-")
    else:
        name = name.replace(" ", "")

    if not name:
        raise InvalidInputError("Target name is empty")
    if any(char in name for char in ("/", "\\", "\x00")) or ".." in name or name.startswith("."):
        raise InvalidInputError(f"Invalid target name: {name!r}")
    return name


def main():
    """Redirect to the unified Typer CLI (``ql run``)."""
    _sys.argv = [_sys.argv[0], "run"] + _sys.argv[1:]
    from quicklook.cli.app import app

    app()


if __name__ == "__main__":
    main()
