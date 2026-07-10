"""Tests for cancel-flow improvements.

Covers:
- ``TessQuickLook._check_cancel`` raises InterruptedError when cancel is signalled
- The race-free finalize of a queued job cancelled before pickup
"""

import time
import pytest


def test_check_cancel_raises_when_signalled():
    """When cancel_check returns truthy, _check_cancel must raise."""
    from quicklook import tql as tql_mod

    qlook = tql_mod.TessQuickLook.__new__(tql_mod.TessQuickLook)
    qlook.cancel_check = lambda: True
    with pytest.raises(InterruptedError, match="cancelled"):
        qlook._check_cancel()


def test_check_cancel_includes_phase_label():
    from quicklook import tql as tql_mod

    qlook = tql_mod.TessQuickLook.__new__(tql_mod.TessQuickLook)
    qlook.cancel_check = lambda: True
    with pytest.raises(InterruptedError, match=r"before TLS"):
        qlook._check_cancel("before TLS transit search")


def test_check_cancel_noop_when_check_returns_false():
    from quicklook import tql as tql_mod

    qlook = tql_mod.TessQuickLook.__new__(tql_mod.TessQuickLook)
    qlook.cancel_check = lambda: False
    qlook._check_cancel()  # must not raise


def test_check_cancel_noop_when_no_check_provided():
    from quicklook import tql as tql_mod

    qlook = tql_mod.TessQuickLook.__new__(tql_mod.TessQuickLook)
    qlook.cancel_check = None
    qlook._check_cancel("anywhere")  # must not raise


def test_cancelled_while_queued_finalizes_history(monkeypatch):
    """A job whose cancel_event was set before the worker picked it up
    must be finalized: status="cancelled", finished_at set, history row
    persisted. Previously this silently returned without recording."""
    from threading import Event
    from quicklook.app import app as app_mod

    saved = {}

    def _fake_save(**kwargs):
        saved.update(kwargs)

    monkeypatch.setattr(app_mod, "_save_job_history", _fake_save)

    cancel_event = Event()
    cancel_event.set()  # cancelled while queued

    name = "TIC-TEST_s99"
    submitted = time.time() - 1
    with app_mod.jobs_lock:
        app_mod.jobs[name] = {
            "status": "queued",
            "log_file": "",
            "error": "",
            "cancel_event": cancel_event,
            "submitted_at": submitted,
            "params": {"pipeline": "spoc"},
            "step_times": {},
        }

    try:
        app_mod.run_quicklook_background(name=name, cancel_event=cancel_event)

        with app_mod.jobs_lock:
            info = app_mod.jobs[name]
        assert info["status"] == "cancelled"
        assert "finished_at" in info
        assert info["finished_at"] >= submitted

        assert saved["name"] == name
        assert saved["status"] == "cancelled"
        assert saved["submitted_at"] == submitted
        assert saved["finished_at"] == info["finished_at"]
    finally:
        with app_mod.jobs_lock:
            app_mod.jobs.pop(name, None)


def test_queued_job_not_marked_running_when_status_already_cancelled(monkeypatch):
    """If the cancel route flipped status to 'cancelled' before the worker
    pulled the job (without yet setting cancel_event, hypothetically),
    the worker must not overwrite that status with 'running'."""
    from threading import Event
    from quicklook.app import app as app_mod

    monkeypatch.setattr(app_mod, "_save_job_history", lambda **kw: None)

    cancel_event = Event()  # NOT set — we're testing the status path
    name = "TIC-TEST_s100"
    with app_mod.jobs_lock:
        app_mod.jobs[name] = {
            "status": "cancelled",  # pre-flipped by cancel route
            "log_file": "",
            "error": "",
            "cancel_event": cancel_event,
            "submitted_at": time.time(),
            "params": {},
            "step_times": {},
        }

    try:
        app_mod.run_quicklook_background(name=name, cancel_event=cancel_event)
        with app_mod.jobs_lock:
            assert app_mod.jobs[name]["status"] == "cancelled"
    finally:
        with app_mod.jobs_lock:
            app_mod.jobs.pop(name, None)
