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

SCHEMA_VERSION = 2

# Seconds to wait on a locked database before raising; the GUI server is
# threaded, so a job finishing while /tls-summary writes cache rows can
# briefly contend for the write lock.
_BUSY_TIMEOUT_S = 10.0

# A cache entry: (absolute path, mtime, size in bytes, extracted fields).
TLSCacheEntry = tuple[str, float, int, dict[str, Any]]


def _connect(path: Path) -> sqlite3.Connection:
    return sqlite3.connect(str(path), timeout=_BUSY_TIMEOUT_S)


def _create_job_history(conn: sqlite3.Connection) -> None:
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


def _create_tls_summary_cache(conn: sqlite3.Connection) -> None:
    # One row per *_tls.h5 file. ``directory`` is stored separately so a
    # page load can fetch one folder's rows with an indexed equality match.
    conn.execute(
        """CREATE TABLE IF NOT EXISTS tls_summary_cache (
            path      TEXT PRIMARY KEY,
            directory TEXT NOT NULL,
            mtime     REAL NOT NULL,
            size      INTEGER NOT NULL,
            fields    TEXT NOT NULL
        )"""
    )
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_tls_summary_cache_dir ON tls_summary_cache (directory)"
    )


def _migrate(conn: sqlite3.Connection) -> None:
    current = conn.execute("PRAGMA user_version").fetchone()[0]
    if current == SCHEMA_VERSION:
        return
    if current > SCHEMA_VERSION:
        # Never silently downgrade.
        raise RuntimeError(
            f"jobs.db user_version={current} is newer than SCHEMA_VERSION={SCHEMA_VERSION}; "
            "refusing to run against an unknown future schema."
        )
    # Step forward one version at a time so every older database follows the
    # same path. Add a new ``if current == N`` block for each future bump.
    if current == 0:
        _create_job_history(conn)
        current = 1
    if current == 1:
        _create_tls_summary_cache(conn)
        current = 2
    conn.execute(f"PRAGMA user_version = {current}")
    conn.commit()


def _parent_dir(path: str) -> str:
    return str(Path(path).parent)


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

    # -- TLS summary cache --------------------------------------------------

    def tls_cache_lookup(self, directory: str) -> dict[str, tuple[float, int, dict[str, Any]]]:
        """Return ``{path: (mtime, size, fields)}`` for files directly in *directory*."""
        conn = _connect(self.path)
        try:
            rows = conn.execute(
                "SELECT path, mtime, size, fields FROM tls_summary_cache WHERE directory = ?",
                (directory,),
            ).fetchall()
        finally:
            conn.close()
        return {path: (mtime, size, json.loads(fields)) for path, mtime, size, fields in rows}

    def tls_cache_store(self, entries: list[TLSCacheEntry]) -> None:
        """Insert or replace cache rows in a single transaction."""
        if not entries:
            return
        conn = _connect(self.path)
        try:
            conn.executemany(
                "INSERT OR REPLACE INTO tls_summary_cache "
                "(path, directory, mtime, size, fields) VALUES (?, ?, ?, ?, ?)",
                [
                    (path, _parent_dir(path), mtime, size, json.dumps(fields))
                    for path, mtime, size, fields in entries
                ],
            )
            conn.commit()
        finally:
            conn.close()

    def tls_cache_prune(self, directory: str, keep: set[str]) -> None:
        """Delete rows for files in *directory* whose path is not in *keep*."""
        cached = self.tls_cache_lookup(directory)
        stale = [(p,) for p in cached if p not in keep]
        if not stale:
            return
        conn = _connect(self.path)
        try:
            conn.executemany("DELETE FROM tls_summary_cache WHERE path = ?", stale)
            conn.commit()
        finally:
            conn.close()
