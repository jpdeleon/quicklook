"""Tests for the GUI job-history SQLite persistence layer."""

import sqlite3
import pytest

from quicklook.app.jobs_db import JobsDB, SCHEMA_VERSION


def test_init_creates_schema(tmp_path):
    db = JobsDB(tmp_path / "jobs.db")
    conn = sqlite3.connect(str(db.path))
    try:
        version = conn.execute("PRAGMA user_version").fetchone()[0]
        tables = [
            r[0]
            for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()
        ]
    finally:
        conn.close()
    assert version == SCHEMA_VERSION
    assert "job_history" in tables


def test_save_and_read_round_trip(tmp_path):
    db = JobsDB(tmp_path / "jobs.db")
    db.save(
        name="TIC-100_s1",
        status="done",
        error="",
        params={"pipeline": "spoc", "sector": 1},
        submitted_at=1.0,
        finished_at=10.0,
        step_times={"Transit search (TLS)": 7.5, "Done": 0.1},
    )
    times = db.recent_step_times(limit=5)
    assert times == [{"Transit search (TLS)": 7.5, "Done": 0.1}]


def test_recent_step_times_only_returns_done(tmp_path):
    db = JobsDB(tmp_path / "jobs.db")
    db.save("A", "done", finished_at=1.0, step_times={"x": 1.0})
    db.save("B", "cancelled", finished_at=2.0, step_times={"x": 2.0})
    db.save("C", "error", finished_at=3.0, step_times={"x": 3.0})
    db.save("D", "done", finished_at=4.0, step_times={"x": 4.0})
    times = db.recent_step_times(limit=10)
    # Should include only A and D, ordered by finished_at DESC.
    assert times == [{"x": 4.0}, {"x": 1.0}]


def test_recent_step_times_respects_limit(tmp_path):
    db = JobsDB(tmp_path / "jobs.db")
    for i in range(5):
        db.save(f"N{i}", "done", finished_at=float(i), step_times={"x": float(i)})
    times = db.recent_step_times(limit=2)
    assert len(times) == 2
    # Highest finished_at first.
    assert times[0]["x"] == 4.0
    assert times[1]["x"] == 3.0


def test_save_is_idempotent_per_name(tmp_path):
    """Same name overwrites — uses INSERT OR REPLACE."""
    db = JobsDB(tmp_path / "jobs.db")
    db.save("X", "queued", submitted_at=1.0, finished_at=None, step_times={"a": 1.0})
    db.save("X", "done", submitted_at=1.0, finished_at=2.0, step_times={"b": 2.0})
    times = db.recent_step_times(limit=5)
    assert times == [{"b": 2.0}]


def test_save_with_defaults(tmp_path):
    """error, params, step_times all default safely (no None crashes)."""
    db = JobsDB(tmp_path / "jobs.db")
    db.save("Y", "done", finished_at=1.0)
    conn = sqlite3.connect(str(db.path))
    try:
        row = conn.execute(
            "SELECT error, params, step_times FROM job_history WHERE name=?",
            ("Y",),
        ).fetchone()
    finally:
        conn.close()
    assert row == ("", "{}", "{}")


def test_delete_removes_row(tmp_path):
    db = JobsDB(tmp_path / "jobs.db")
    db.save("Z", "done", finished_at=1.0, step_times={"x": 1.0})
    db.delete("Z")
    assert db.recent_step_times() == []


def test_delete_missing_name_is_noop(tmp_path):
    db = JobsDB(tmp_path / "jobs.db")
    # Must not raise even though "nope" was never inserted.
    db.delete("nope")


def test_migrate_rejects_future_version(tmp_path):
    path = tmp_path / "future.db"
    conn = sqlite3.connect(str(path))
    try:
        conn.execute(f"PRAGMA user_version = {SCHEMA_VERSION + 1}")
        conn.commit()
    finally:
        conn.close()
    with pytest.raises(RuntimeError, match="newer than SCHEMA_VERSION"):
        JobsDB(path)


def test_migrate_is_idempotent(tmp_path):
    path = tmp_path / "jobs.db"
    JobsDB(path)  # creates schema
    JobsDB(path)  # second init must not error or alter the schema
    conn = sqlite3.connect(str(path))
    try:
        assert conn.execute("PRAGMA user_version").fetchone()[0] == SCHEMA_VERSION
    finally:
        conn.close()


def _tables(path):
    conn = sqlite3.connect(str(path))
    try:
        return {
            r[0]
            for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()
        }
    finally:
        conn.close()


def test_init_creates_tls_summary_cache_table(tmp_path):
    db = JobsDB(tmp_path / "jobs.db")
    assert "tls_summary_cache" in _tables(db.path)


def test_migrate_v1_to_current_keeps_job_history(tmp_path):
    path = tmp_path / "v1.db"
    conn = sqlite3.connect(str(path))
    try:
        conn.execute(
            "CREATE TABLE job_history (name TEXT PRIMARY KEY, status TEXT NOT NULL, "
            "error TEXT DEFAULT '', params TEXT DEFAULT '{}', submitted_at REAL, "
            "finished_at REAL, step_times TEXT DEFAULT '{}')"
        )
        conn.execute("INSERT INTO job_history (name, status) VALUES ('old', 'done')")
        conn.execute("PRAGMA user_version = 1")
        conn.commit()
    finally:
        conn.close()

    JobsDB(path)

    conn = sqlite3.connect(str(path))
    try:
        assert conn.execute("PRAGMA user_version").fetchone()[0] == SCHEMA_VERSION
        assert conn.execute("SELECT name FROM job_history").fetchall() == [("old",)]
    finally:
        conn.close()
    assert "tls_summary_cache" in _tables(path)


def test_tls_cache_store_and_lookup_round_trip(tmp_path):
    db = JobsDB(tmp_path / "jobs.db")
    db.tls_cache_store([("/d/a_tls.h5", 10.0, 100, {"SDE": 9.5, "toi": None})])
    got = db.tls_cache_lookup("/d")
    assert got == {"/d/a_tls.h5": (10.0, 100, {"SDE": 9.5, "toi": None})}


def test_tls_cache_store_replaces_existing_entry(tmp_path):
    db = JobsDB(tmp_path / "jobs.db")
    db.tls_cache_store([("/d/a_tls.h5", 10.0, 100, {"SDE": 1})])
    db.tls_cache_store([("/d/a_tls.h5", 20.0, 200, {"SDE": 2})])
    assert db.tls_cache_lookup("/d")["/d/a_tls.h5"] == (20.0, 200, {"SDE": 2})


def test_tls_cache_lookup_is_scoped_to_one_directory(tmp_path):
    db = JobsDB(tmp_path / "jobs.db")
    db.tls_cache_store(
        [
            ("/d/a_tls.h5", 1.0, 1, {}),
            ("/d/sub/b_tls.h5", 1.0, 1, {}),
            ("/dx/c_tls.h5", 1.0, 1, {}),
        ]
    )
    assert set(db.tls_cache_lookup("/d")) == {"/d/a_tls.h5"}


def test_tls_cache_prune_drops_only_missing_files_in_directory(tmp_path):
    db = JobsDB(tmp_path / "jobs.db")
    db.tls_cache_store(
        [
            ("/d/keep_tls.h5", 1.0, 1, {}),
            ("/d/gone_tls.h5", 1.0, 1, {}),
            ("/other/x_tls.h5", 1.0, 1, {}),
        ]
    )
    db.tls_cache_prune("/d", keep={"/d/keep_tls.h5"})
    assert set(db.tls_cache_lookup("/d")) == {"/d/keep_tls.h5"}
    assert set(db.tls_cache_lookup("/other")) == {"/other/x_tls.h5"}
