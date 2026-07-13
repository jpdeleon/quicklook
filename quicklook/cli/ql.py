#!/usr/bin/env python
"""Target-name sanitization shared by the CLI and web GUI.

The user-facing commands live in :mod:`quicklook.cli.app` (the unified Typer
app exposed as the ``quicklook`` and ``ql`` console scripts).
"""

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
