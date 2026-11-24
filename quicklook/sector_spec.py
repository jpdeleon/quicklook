"""Parse the sector choices shared by the CLI and the web GUI.

A sector spec is one of:

- ``-1``: the latest available sector
- a sector number, e.g. ``57``
- ``"all"``: stitch every sector of one pipeline and cadence
- a list of at least two sectors, e.g. ``"57,58,77-84"``, stitched together

``TessQuickLook(sector=...)`` accepts all four forms.
"""

from __future__ import annotations

import re
from typing import Union

LATEST_SECTOR = -1
MAX_SECTORS = 1024
SectorSpec = Union[int, str, list]
_RANGE = re.compile(r"(\d+)(?:-(\d+))?")


def parse_sector_ranges(text: str) -> list:
    """Sorted, unique sectors from text such as ``"66-68,93"``.

    Raises ``ValueError`` on anything but positive numbers and ascending ranges.
    """
    sectors = set()
    for part in str(text).split(","):
        part = part.strip()
        match = _RANGE.fullmatch(part)
        if not match:
            raise ValueError(
                f"sectors must look like '66-68,93' (numbers and ranges), got {text!r}"
            )
        first = int(match.group(1))
        last = int(match.group(2)) if match.group(2) else first
        if first < 1 or last < first:
            raise ValueError(f"invalid sector range {part!r}")
        if last - first + 1 > MAX_SECTORS:
            raise ValueError(f"select at most {MAX_SECTORS} sectors")
        for sector in range(first, last + 1):
            sectors.add(sector)
            if len(sectors) > MAX_SECTORS:
                raise ValueError(f"select at most {MAX_SECTORS} sectors")
    return sorted(sectors)


def parse_sector_spec(value) -> SectorSpec:
    """``-1``, a sector number, ``"all"`` or a list of >= 2 sectors from user input.

    Empty input means the latest sector. A list that names a single sector
    collapses to that number. Raises ``ValueError`` on anything else.
    """
    if value is None or isinstance(value, bool):
        return LATEST_SECTOR
    if isinstance(value, int):
        return _single(value)
    if isinstance(value, (list, tuple)):
        value = ",".join(str(v) for v in value)
    text = str(value).strip().lower()
    if text in ("", str(LATEST_SECTOR)):
        return LATEST_SECTOR
    if text == "all":
        return "all"
    sectors = parse_sector_ranges(text)
    return sectors[0] if len(sectors) == 1 else sectors


def _single(sector: int) -> int:
    if sector == LATEST_SECTOR or sector >= 1:
        return sector
    raise ValueError(f"sector must be -1 (latest) or a positive number, got {sector}")


def is_multi_sector(spec: SectorSpec) -> bool:
    """True when ``spec`` stitches several sectors (``"all"`` or a list)."""
    return spec == "all" or isinstance(spec, list)
