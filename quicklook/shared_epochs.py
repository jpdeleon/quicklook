"""Single events shared by unrelated stars: a batch-level systematics check.

Scattered light near orbit gaps and momentum dumps dim many stars at the same
time, and each target's single-event search sees it as a transit. In the
young-host TGLC run (issue #42, step 6) 13 of 31 mono/duo hosts in sector 14
had events within hours of each other. Comparing event times across a batch of
quicklook results (``*_tls.h5`` written with ``--single-event-search``)
catches this; a single target cannot.

An event is *shared* when events from at least ``min_targets`` distinct targets
(itself included) fall within ``tol`` days of it. A mono/duo target whose
events are all shared is re-labelled ``systematic`` in ``event_class_vetted``.
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import numpy as np
import pandas as pd
from loguru import logger

from quicklook import h5io

SHARED_TOL_DAYS = 0.25
SHARED_MIN_TARGETS = 3
SINGLE_EVENT_CLASSES = ("mono", "duo")
_EVENT_KEYS = ("event_t0s", "event_snrs", "event_depths", "event_durations")
EVENT_COLUMNS = ["ticid", "sector", "event_class", "t0", "snr", "depth", "duration", "file"]


def _sector_label(sector) -> str:
    if sector is None:
        return ""
    values = np.atleast_1d(sector)
    return "+".join(str(int(s)) for s in values if np.isfinite(s))


def load_events(paths: Iterable) -> pd.DataFrame:
    """One row per single event in the given ``*_tls.h5`` files.

    Files without single-event results contribute no rows. Unreadable files
    are skipped with a warning.
    """
    rows = []
    for path in paths:
        try:
            data = h5io.load(path, keys=_EVENT_KEYS)
        except (OSError, KeyError, ValueError) as e:
            logger.warning(f"Skipping {path}: {e}")
            continue
        t0s = np.atleast_1d(np.asarray(data.get("event_t0s", []), dtype=float))
        if t0s.size == 0:
            continue
        cols = {
            k: np.resize(np.asarray(data.get(f"event_{k}s", np.nan), dtype=float), t0s.size)
            for k in ("snr", "depth", "duration")
        }
        for i, t0 in enumerate(t0s):
            rows.append(
                {
                    "ticid": data.get("ticid"),
                    "sector": _sector_label(data.get("sector")),
                    "event_class": data.get("event_class", ""),
                    "t0": float(t0),
                    "snr": float(cols["snr"][i]),
                    "depth": float(cols["depth"][i]),
                    "duration": float(cols["duration"][i]),
                    "file": str(path),
                }
            )
    return pd.DataFrame(rows, columns=EVENT_COLUMNS)


def flag_shared_events(
    events: pd.DataFrame,
    tol: float = SHARED_TOL_DAYS,
    min_targets: int = SHARED_MIN_TARGETS,
) -> pd.DataFrame:
    """Copy of ``events`` with ``n_targets_nearby`` and ``shared`` columns.

    ``n_targets_nearby`` counts distinct targets (this one included) with an
    event within ``tol`` days.
    """
    out = events.copy()
    if out.empty:
        out["n_targets_nearby"] = pd.Series(dtype=int)
        out["shared"] = pd.Series(dtype=bool)
        return out
    t = out["t0"].to_numpy(dtype=float)
    tic = out["ticid"].to_numpy()
    order = np.argsort(t)
    ts, tics = t[order], tic[order]
    lo = np.searchsorted(ts, ts - tol, side="left")
    hi = np.searchsorted(ts, ts + tol, side="right")
    counts = np.array([len(set(tics[a:b])) for a, b in zip(lo, hi)])
    nearby = np.empty_like(counts)
    nearby[order] = counts
    out["n_targets_nearby"] = nearby
    out["shared"] = nearby >= min_targets
    return out


def vet_targets(flagged: pd.DataFrame) -> pd.DataFrame:
    """Per-target summary with ``event_class_vetted``.

    Mono/duo targets whose events are all shared become ``systematic``;
    every other class is kept.
    """
    cols = ["ticid", "sector", "event_class", "n_events", "n_shared", "event_class_vetted"]
    if flagged.empty:
        return pd.DataFrame(columns=cols)
    grouped = flagged.groupby(["ticid", "sector", "event_class"], dropna=False)
    summary = grouped.agg(n_events=("t0", "size"), n_shared=("shared", "sum")).reset_index()
    summary["n_shared"] = summary["n_shared"].astype(int)
    all_shared = summary["n_shared"] == summary["n_events"]
    single = summary["event_class"].isin(SINGLE_EVENT_CLASSES)
    summary["event_class_vetted"] = np.where(
        single & all_shared, "systematic", summary["event_class"]
    )
    return summary[cols]


def shared_epochs_report(
    input_dir,
    tol: float = SHARED_TOL_DAYS,
    min_targets: int = SHARED_MIN_TARGETS,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(events, targets) tables for every ``*_tls.h5`` under ``input_dir``."""
    paths = sorted(Path(input_dir).glob("*_tls.h5"))
    if not paths:
        raise FileNotFoundError(f"No *_tls.h5 files found in {input_dir}")
    events = flag_shared_events(load_events(paths), tol=tol, min_targets=min_targets)
    return events, vet_targets(events)


__all__ = [
    "SHARED_MIN_TARGETS",
    "SHARED_TOL_DAYS",
    "flag_shared_events",
    "load_events",
    "shared_epochs_report",
    "vet_targets",
]
