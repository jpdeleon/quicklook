"""Tiny persistence layer for the GUI job-history SQLite database.

One module owns the schema, the version pragma, and the CRUD it needs.
The app code talks to ``JobsDB`` instead of cracking open sqlite3
inline, so future schema changes ship from one place — bump
``SCHEMA_VERSION``, add a migration branch in ``_migrate``, and the
next process start handles it.
"""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path
from typing import Any

SCHEMA_VERSION = 1


def _connect(path: Path) -> sqlite3.Connection:
    return sqlite3.connect(str(path))


def _migrate(conn: sqlite3.Connection) -> None:
    current = conn.execute("PRAGMA user_version").fetchone()[0]
    if current == SCHEMA_VERSION:
        return
    if current == 0:
        conn.execute(
            """CREATE TABLE IF NOT EXISTS job_history (
                name         TEXT PRIMARY KEY,
                status       TEXT NOT NULL,
                error        TEXT DEFAULT '',
                params       TEXT DEFAULT '{}',
                submitted_at REAL,
                finished_at  REAL,
                step_times   TEXT DEFAULT '{}'
            )"""
        )
        conn.execute(f"PRAGMA user_version = {SCHEMA_VERSION}")
        conn.commit()
        return
    # Future versions: add elif current == N branches that migrate
    # forward to SCHEMA_VERSION. Never silently downgrade.
    raise RuntimeError(
        f"jobs.db user_version={current} is newer than SCHEMA_VERSION={SCHEMA_VERSION}; "
        "refusing to run against an unknown future schema."
    )


class JobsDB:
    """Thin wrapper around the job_history table. Stateless per call."""

    def __init__(self, path: Path):
        self.path = path
        self.init()

    def init(self) -> None:
        conn = _connect(self.path)
        try:
            _migrate(conn)
        finally:
            conn.close()

    def save(
        self,
        name: str,
        status: str,
        error: str = "",
        params: dict[str, Any] | None = None,
        submitted_at: float | None = None,
        finished_at: float | None = None,
        step_times: dict[str, float] | None = None,
    ) -> None:
        conn = _connect(self.path)
        try:
            conn.execute(
                "INSERT OR REPLACE INTO job_history "
                "(name, status, error, params, submitted_at, finished_at, step_times) "
                "VALUES (?, ?, ?, ?, ?, ?, ?)",
                (
                    name,
                    status,
                    error,
                    json.dumps(params or {}),
                    submitted_at,
                    finished_at,
                    json.dumps(step_times or {}),
                ),
            )
            conn.commit()
        finally:
            conn.close()

    def recent_step_times(self, limit: int = 20) -> list[dict[str, float]]:
        """Return parsed step_times dicts from the most recent successful jobs."""
        conn = _connect(self.path)
        try:
            rows = conn.execute(
                "SELECT step_times FROM job_history WHERE status='done' "
                "ORDER BY finished_at DESC LIMIT ?",
                (limit,),
            ).fetchall()
        finally:
            conn.close()
        return [json.loads(st) if st else {} for (st,) in rows]

    def delete(self, name: str) -> None:
        conn = _connect(self.path)
        try:
            conn.execute("DELETE FROM job_history WHERE name = ?", (name,))
            conn.commit()
        finally:
            conn.close()
