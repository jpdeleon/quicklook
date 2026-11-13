"""Tests for the cached TLS-summary scan behind the /tls-summary page."""

import os
import sqlite3

import numpy as np
import pytest

from quicklook import h5io
from quicklook.app import tls_summary
from quicklook.app.jobs_db import JobsDB


def _write_tls(directory, name, sde=10.0, mtime=None):
    path = directory / f"{name}_s12_pdcsap_sc_tls.h5"
    h5io.save(
        path,
        {
            "toiid": 700.01,
            "period": 2.5,
            "duration": 0.1,
            "depth": 0.999,
            "SDE": sde,
            "Prot_gls": np.array([3.2, 0.1]),
            "per_transit_count": np.array([5, 6]),
            "flux_raw": np.ones(5_000),
        },
    )
    if mtime is not None:
        os.utime(path, (mtime, mtime))
    return path


@pytest.fixture
def load_counter(monkeypatch):
    """Count h5io.load calls made by the summary scan."""
    calls = []
    real_load = h5io.load

    def counting_load(path, keys=None):
        calls.append(os.path.basename(str(path)))
        return real_load(path, keys=keys)

    monkeypatch.setattr(tls_summary.h5io, "load", counting_load)
    return calls


def test_extract_fields_derives_table_values():
    fields = tls_summary.extract_fields(
        {
            "toiid": 700.01,
            "period": 2.5,
            "duration": 0.1,
            "depth": 0.999,
            "Prot_gls": np.array([3.2, 0.1]),
            "per_transit_count": np.array([5, 6]),
        },
        "TOI700_s12_pdcsap_sc_tls",
    )
    assert fields["toi"] == pytest.approx(700.01)
    assert fields["Porb"] == pytest.approx(2.5)
    assert fields["duration_hr"] == pytest.approx(2.4)
    assert fields["depth_ppt"] == pytest.approx(1.0)
    assert fields["Prot_gls"] == pytest.approx(3.2)
    assert fields["per_transit_count"] == "5,6"
    assert fields["pipeline"] == "spoc"  # filename fallback
    assert fields["flux_type"] == "pdcsap"
    assert fields["exptime"] == 120  # short-cadence SPOC hint


def test_scan_loads_only_summary_keys(tmp_path, monkeypatch):
    _write_tls(tmp_path, "A")
    seen_keys = []
    real_load = h5io.load

    def spy(path, keys=None):
        seen_keys.append(keys)
        return real_load(path, keys=keys)

    monkeypatch.setattr(tls_summary.h5io, "load", spy)
    tls_summary.scan(tmp_path, JobsDB(tmp_path / "jobs.db"))
    assert seen_keys and seen_keys[0] is not None
    assert "flux_raw" not in seen_keys[0]


def test_second_scan_reads_no_files(tmp_path, load_counter):
    _write_tls(tmp_path, "A")
    _write_tls(tmp_path, "B")
    db = JobsDB(tmp_path / "jobs.db")

    first = tls_summary.scan(tmp_path, db)
    assert len(load_counter) == 2
    second = tls_summary.scan(tmp_path, db)

    assert len(load_counter) == 2  # served entirely from the cache
    assert [e.path.name for e in second.entries] == [e.path.name for e in first.entries]
    assert second.entries[0].fields == first.entries[0].fields


def test_modified_file_is_reloaded(tmp_path, load_counter):
    path = _write_tls(tmp_path, "A", sde=10.0, mtime=1_000_000)
    db = JobsDB(tmp_path / "jobs.db")
    tls_summary.scan(tmp_path, db)

    _write_tls(tmp_path, "A", sde=42.0, mtime=2_000_000)
    result = tls_summary.scan(tmp_path, db)

    assert load_counter.count(path.name) == 2
    assert result.entries[0].fields["SDE"] == pytest.approx(42.0)


def test_deleted_file_is_dropped_and_pruned(tmp_path, load_counter):
    keep = _write_tls(tmp_path, "A")
    gone = _write_tls(tmp_path, "B")
    db = JobsDB(tmp_path / "jobs.db")
    tls_summary.scan(tmp_path, db)

    gone.unlink()
    result = tls_summary.scan(tmp_path, db)

    assert [e.path.name for e in result.entries] == [keep.name]
    assert set(db.tls_cache_lookup(str(tmp_path.resolve()))) == {str(keep.resolve())}


def test_entries_are_sorted_newest_first(tmp_path):
    _write_tls(tmp_path, "OLD", mtime=1_000_000)
    _write_tls(tmp_path, "NEW", mtime=2_000_000)
    result = tls_summary.scan(tmp_path, JobsDB(tmp_path / "jobs.db"))
    assert [e.path.name.split("_")[0] for e in result.entries] == ["NEW", "OLD"]


def test_unreadable_file_is_skipped_and_not_cached(tmp_path, load_counter):
    bad = tmp_path / "broken_s1_pdcsap_sc_tls.h5"
    bad.write_bytes(b"not an hdf5 file")
    db = JobsDB(tmp_path / "jobs.db")

    result = tls_summary.scan(tmp_path, db)

    assert result.entries == []
    assert [s["name"] for s in result.skipped] == [bad.name]
    assert db.tls_cache_lookup(str(tmp_path.resolve())) == {}


def test_stale_fields_version_forces_reload(tmp_path, load_counter, monkeypatch):
    _write_tls(tmp_path, "A")
    db = JobsDB(tmp_path / "jobs.db")
    tls_summary.scan(tmp_path, db)

    monkeypatch.setattr(tls_summary, "FIELDS_VERSION", tls_summary.FIELDS_VERSION + 1)
    tls_summary.scan(tmp_path, db)

    assert len(load_counter) == 2


class _BrokenCache:
    def tls_cache_lookup(self, directory):
        raise sqlite3.OperationalError("database is locked")

    def tls_cache_store(self, entries):
        raise sqlite3.OperationalError("database is locked")

    def tls_cache_prune(self, directory, keep):
        raise sqlite3.OperationalError("database is locked")


def test_cache_failure_falls_back_to_direct_reads(tmp_path):
    _write_tls(tmp_path, "A")
    result = tls_summary.scan(tmp_path, _BrokenCache())
    assert len(result.entries) == 1
    assert result.entries[0].fields["SDE"] == pytest.approx(10.0)


def test_route_renders_cached_rows(tmp_path, monkeypatch):
    import quicklook.app.app as app_module

    _write_tls(tmp_path, "TOI700")
    (tmp_path / "TOI700_s12_pdcsap_sc.png").write_bytes(b"png")
    monkeypatch.setattr(app_module, "OUTPUT_DIR", tmp_path)
    monkeypatch.setattr(app_module, "_jobs_db", JobsDB(tmp_path / "jobs.db"))

    with app_module.app.test_client() as client:
        first = client.get("/tls-summary")
        second = client.get("/tls-summary")

    assert first.status_code == 200
    assert b"TOI700_s12_pdcsap_sc_tls.h5" in first.data
    assert b"/static/outputs/TOI700_s12_pdcsap_sc.png" in first.data
    assert first.data == second.data
