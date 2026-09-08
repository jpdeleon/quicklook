"""Tests for the unified Typer CLI (quicklook.cli.app)."""

import sys

import pytest
from click import unstyle
from typer.testing import CliRunner
from quicklook.cli.app import app

runner = CliRunner()


def test_app_help_succeeds():
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "Quick look analysis" in result.output


def test_app_no_args_shows_help():
    result = runner.invoke(app, [])
    assert result.exit_code == 2
    assert "Commands" in result.output or "Usage" in result.output


def test_run_help_succeeds():
    result = runner.invoke(app, ["run", "--help"])
    assert result.exit_code == 0
    assert "Target name" in result.output


def test_run_requires_name():
    result = runner.invoke(app, ["run"])
    assert result.exit_code == 2
    assert "Missing" in result.output


def test_read_tls_help_succeeds():
    result = runner.invoke(app, ["read-tls", "--help"])
    assert result.exit_code == 0
    # Typer/Click render the required-argument placeholder differently across
    # versions (e.g. "INPUT_DIR" vs "{input_dir}"); match case-insensitively.
    assert "input_dir" in result.output.lower()
    assert "Parameter to sort by" in result.output


def test_read_tls_fails_without_input():
    result = runner.invoke(app, ["read-tls"])
    assert result.exit_code == 2
    assert "Missing" in result.output


def test_rank_tls_help_succeeds():
    result = runner.invoke(app, ["rank-tls", "--help"])
    assert result.exit_code == 0
    assert "input_dir" in result.output.lower()
    assert "Minimum SDE" in result.output


def test_rank_tls_fails_without_input():
    result = runner.invoke(app, ["rank-tls"])
    assert result.exit_code == 2
    assert "Missing" in result.output


def test_run_with_mutually_exclusive_flags():
    result = runner.invoke(
        app,
        ["run", "--name", "TOI-1234", "--each-sector", "--each-pipeline"],
    )
    assert result.exit_code == 2
    assert "mutually exclusive" in result.output.lower()


def test_run_with_invalid_target_name():
    result = runner.invoke(app, ["run", "--name", "../etc/passwd"])
    assert result.exit_code == 2
    assert "error" in result.output.lower() or "invalid" in result.output.lower()


def test_run_default_flag_values():
    result = runner.invoke(app, ["run", "--name", "TOI-1234", "--help"])
    assert result.exit_code == 0
    assert "biweight" in result.output
    assert "periodic_auto" in result.output


@pytest.mark.skipif(
    sys.platform != "linux",
    reason="the DISPLAY-based headless guard only applies on Linux; on macOS/Windows "
    "matplotlib doesn't need $DISPLAY, so this run would fall through to the real, "
    "network-touching pipeline instead of aborting",
)
def test_run_aborts_on_headless_without_save():
    import os

    had_display = os.environ.pop("DISPLAY", None)
    try:
        result = runner.invoke(app, ["run", "--name", "TOI-6109"])
        assert result.exit_code == 1
        assert "headless" in result.output.lower()
    finally:
        if had_display is not None:
            os.environ["DISPLAY"] = had_display


def test_run_with_save_skips_headless_check(monkeypatch):
    """With --save, a headless run reaches the analysis pipeline."""
    import quicklook.tql as tql_mod

    calls = {}

    class FakeQuickLook:
        def __init__(self, **kwargs):
            calls.update(kwargs)

        def plot_tql(self):
            calls["plot_tql"] = True

    monkeypatch.setattr(tql_mod, "TessQuickLook", FakeQuickLook)
    monkeypatch.delenv("DISPLAY", raising=False)

    result = runner.invoke(app, ["run", "--name", "TOI-6109", "--save"])

    assert result.exit_code == 0, result.output
    assert calls["savefig"] is True
    assert calls["savetls"] is True
    assert calls["plot_tql"] is True


def test_gui_help_succeeds():
    result = runner.invoke(app, ["gui", "--help"])
    output = unstyle(result.output)
    assert result.exit_code == 0
    assert "web" in output.lower() or "gui" in output.lower()
    assert "--host" in output
    assert "--port" in output


def test_gui_invokes_run_gui_with_options(monkeypatch):
    """`quicklook gui --host ... --port ...` forwards to run_gui without
    actually starting the blocking Flask server."""
    import quicklook.app.app as gui_module

    calls = {}

    def fake_run_gui(host="127.0.0.1", port=5000, debug=None):
        calls.update(host=host, port=port, debug=debug)

    monkeypatch.setattr(gui_module, "run_gui", fake_run_gui)
    result = runner.invoke(app, ["gui", "--host", "0.0.0.0", "--port", "8080"])
    assert result.exit_code == 0, result.output
    assert calls["host"] == "0.0.0.0"
    assert calls["port"] == 8080


def test_gui_debug_flag_forwarded(monkeypatch):
    import quicklook.app.app as gui_module

    calls = {}

    def fake_run_gui(host="127.0.0.1", port=5000, debug=None):
        calls.update(host=host, port=port, debug=debug)

    monkeypatch.setattr(gui_module, "run_gui", fake_run_gui)
    result = runner.invoke(app, ["gui", "--debug"])
    assert result.exit_code == 0, result.output
    assert calls["debug"] is True
