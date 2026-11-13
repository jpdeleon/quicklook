"""Cached scan of ``*_tls.h5`` results for the GUI's /tls-summary page.

Opening every result file on each page load made the page take over a
minute on a cold disk cache (hundreds of files, about 1 GB of light-curve
arrays, on spinning RAID). This module:

* reads only the datasets the summary table needs (``SUMMARY_DATASET_KEYS``)
  instead of the full light curves, and
* keeps the extracted values in the ``tls_summary_cache`` table in jobs.db,
  keyed by path, mtime and size, so unchanged files are never reopened.

Any cache failure (locked or corrupt database) falls back to reading the
files directly, so the page still works, just more slowly.
"""

from __future__ import annotations

import json
import os
import re
import sqlite3
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Protocol

import numpy as np
from loguru import logger

from quicklook import h5io
from quicklook.pipelines import HLSP_PIPELINES

# Bump whenever ``extract_fields`` changes what it stores, so rows cached by
# an older version are re-extracted instead of being served stale.
FIELDS_VERSION = 2
_VERSION_KEY = "_v"

TLS_GLOB = "*_tls.h5"

# Array-valued entries the table reads. Scalars live in the JSON attribute,
# which h5io.load always returns, so only datasets need to be listed here.
SUMMARY_DATASET_KEYS = (
    "toiid",
    "ticid",
    "gaiaid",
    "pipeline",
    "flux_type",
    "exptime",
    "sector",
    "period",
    "duration",
    "rp_rs",
    "SDE",
    "snr",
    "FAP",
    "odd_even_mismatch",
    "distinct_transit_count",
    "depth",
    "Prot_gls",
    "amp_gls",
    "power_gls",
    "per_transit_count",
    "simbad_obj",
)

_SPOC_FLUX_TYPES = {"pdcsap", "sap"}
# Typical SPOC 2-minute cadence, used as a hint when exptime was not saved.
_SPOC_SHORT_CADENCE_S = 120


class SummaryCache(Protocol):
    """The subset of ``JobsDB`` the scan relies on."""

    def tls_cache_lookup(self, directory: str) -> dict[str, tuple[float, int, dict[str, Any]]]: ...

    def tls_cache_store(self, entries: list[tuple[str, float, int, dict[str, Any]]]) -> None: ...

    def tls_cache_prune(self, directory: str, keep: set[str]) -> None: ...


@dataclass(frozen=True)
class SummaryEntry:
    path: Path
    mtime: float
    fields: dict[str, Any]


@dataclass(frozen=True)
class SummaryScan:
    entries: list[SummaryEntry] = field(default_factory=list)
    skipped: list[dict[str, str]] = field(default_factory=list)


def parse_tls_filename(stem: str) -> dict[str, str | None]:
    """Best-effort parse of pipeline, flux_type, cadence from a TLS filename.

    Filenames follow `{name}_s{NN}_{lctype}_{c}c[_mask_ephem][_{suffix}]_tls`
    (see TessQuickLook.check_output_file_exists). ``lctype`` is the flux type
    when pipeline=="spoc" (pdcsap/sap), else the pipeline name itself.
    Returns a dict of fallbacks; values are None when not confidently
    recoverable.
    """
    if stem.endswith("_tls"):
        stem = stem[: -len("_tls")]
    tokens = stem.split("_")
    pipeline = flux_type = cadence = None
    sec_idx = next((i for i, t in enumerate(tokens) if re.fullmatch(r"s\d{1,3}", t)), None)
    cad_idx = next((i for i, t in enumerate(tokens) if re.fullmatch(r"[sl]c", t)), None)
    if sec_idx is not None and cad_idx is not None and cad_idx > sec_idx + 1:
        lctype = tokens[sec_idx + 1]
        if lctype in _SPOC_FLUX_TYPES:
            pipeline, flux_type = "spoc", lctype
        elif lctype in HLSP_PIPELINES:
            pipeline = lctype
        cadence = "short" if tokens[cad_idx] == "sc" else "long"
    return {"pipeline": pipeline, "flux_type": flux_type, "cadence": cadence}


def coerce_scalar(val: Any) -> Any:
    """Unwrap 0-d/1-element arrays and numpy scalars into plain Python values."""
    if val is None:
        return None
    if isinstance(val, np.ndarray):
        if val.ndim == 0:
            return val.item()
        if val.size == 1:
            return val.flat[0].item()
        return val
    if isinstance(val, np.generic):
        return val.item()
    if isinstance(val, bytes):
        try:
            return val.decode()
        except UnicodeDecodeError:
            return None
    return val


def _first(seq: Any) -> Any:
    try:
        return coerce_scalar(seq[0]) if seq is not None and len(seq) > 0 else None
    except (TypeError, IndexError, KeyError):
        return coerce_scalar(seq)


def _join_counts(ptc: Any) -> str:
    if ptc is None:
        return ""
    try:
        return ",".join(str(int(x)) for x in ptc)
    except (TypeError, ValueError):
        return str(ptc)


def _json_safe(val: Any) -> Any:
    """Make a value storable as JSON in the cache (numpy -> Python, else str)."""
    val = coerce_scalar(val)
    if isinstance(val, np.ndarray):
        return val.tolist()
    try:
        json.dumps(val)
    except (TypeError, ValueError):
        return str(val)
    return val


def extract_fields(data: dict[str, Any], stem: str) -> dict[str, Any]:
    """Derive the data-dependent table columns for one TLS result file."""
    fallback = parse_tls_filename(stem)
    pipeline = coerce_scalar(data.get("pipeline")) or fallback["pipeline"]
    flux_type = coerce_scalar(data.get("flux_type")) or fallback["flux_type"]
    exptime = coerce_scalar(data.get("exptime"))
    if exptime is None and fallback["cadence"] == "short" and pipeline == "spoc":
        # A hint, not a hard claim: older files did not save exptime.
        exptime = _SPOC_SHORT_CADENCE_S

    depth = coerce_scalar(data.get("depth"))
    duration_days = coerce_scalar(data.get("duration"))
    fields = {
        "toi": coerce_scalar(data.get("toiid")),
        "tic": coerce_scalar(data.get("ticid")),
        "gaiaid": coerce_scalar(data.get("gaiaid")),
        "pipeline": pipeline,
        "flux_type": flux_type,
        "exptime": exptime,
        "sector": coerce_scalar(data.get("sector")),
        "Porb": coerce_scalar(data.get("period")),
        "duration_hr": duration_days * 24.0 if duration_days is not None else None,
        "rp_rs": coerce_scalar(data.get("rp_rs")),
        "SDE": coerce_scalar(data.get("SDE")),
        "snr": coerce_scalar(data.get("snr")),
        "FAP": coerce_scalar(data.get("FAP")),
        "odd_even": coerce_scalar(data.get("odd_even_mismatch")),
        "n_tr": coerce_scalar(data.get("distinct_transit_count")),
        "depth_ppt": (1 - depth) * 1e3 if depth is not None else None,
        "Prot_gls": _first(data.get("Prot_gls")),
        "amp_gls": _first(data.get("amp_gls")),
        "power_gls": _first(data.get("power_gls")),
        "per_transit_count": _join_counts(data.get("per_transit_count")),
        "simbad_obj": coerce_scalar(data.get("simbad_obj")),
        # Single-event search (absent unless --single-event-search was used).
        "event_class": coerce_scalar(data.get("event_class")),
        "n_events": coerce_scalar(data.get("n_events")),
        "p_min": coerce_scalar(data.get("p_min")),
    }
    return {key: _json_safe(val) for key, val in fields.items()}


def _cache_call(fn: Callable[..., Any], *args: Any, default: Any = None) -> Any:
    """Run a cache operation; on database errors log and return *default*."""
    try:
        return fn(*args)
    except sqlite3.Error as e:
        logger.warning(f"TLS summary cache unavailable ({e}); reading files directly")
        return default


def _list_files(directory: Path) -> list[tuple[Path, os.stat_result]]:
    files = []
    for path in directory.glob(TLS_GLOB):
        try:
            files.append((path, path.stat()))
        except OSError:
            continue  # removed between glob and stat
    return files


def _cached_fields(hit: Any, st: os.stat_result) -> dict[str, Any] | None:
    if hit is None:
        return None
    mtime, size, fields = hit
    if mtime != st.st_mtime or size != st.st_size or fields.get(_VERSION_KEY) != FIELDS_VERSION:
        return None
    return {k: v for k, v in fields.items() if k != _VERSION_KEY}


def scan(directory: Path, cache: SummaryCache) -> SummaryScan:
    """Summarize every TLS result in *directory*, newest first.

    Only files that are new or changed since the last scan (or whose cached
    row predates ``FIELDS_VERSION``) are opened.
    """
    directory = Path(directory).resolve()
    dir_key = str(directory)
    files = _list_files(directory)
    cached = _cache_call(cache.tls_cache_lookup, dir_key, default={})

    entries: list[SummaryEntry] = []
    skipped: list[dict[str, str]] = []
    fresh: list[tuple[str, float, int, dict[str, Any]]] = []
    for path, st in files:
        fields = _cached_fields(cached.get(str(path)), st)
        if fields is None:
            try:
                data = h5io.load(str(path), keys=SUMMARY_DATASET_KEYS)
            except Exception as e:  # h5py raises OSError, KeyError, ValueError, ...
                skipped.append({"name": path.name, "error": str(e)})
                continue
            fields = extract_fields(data, path.stem)
            fresh.append(
                (str(path), st.st_mtime, st.st_size, {**fields, _VERSION_KEY: FIELDS_VERSION})
            )
        entries.append(SummaryEntry(path=path, mtime=st.st_mtime, fields=fields))

    _cache_call(cache.tls_cache_store, fresh)
    _cache_call(cache.tls_cache_prune, dir_key, {str(p) for p, _ in files})
    entries.sort(key=lambda e: e.mtime, reverse=True)
    return SummaryScan(entries=entries, skipped=skipped)
