"""Tests for the unified Typer CLI (quicklook.cli.app)."""

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
    assert "INPUT_DIR" in result.output
    assert "Parameter to sort by" in result.output


def test_read_tls_fails_without_input():
    result = runner.invoke(app, ["read-tls"])
    assert result.exit_code == 2
    assert "Missing" in result.output


def test_rank_tls_help_succeeds():
    result = runner.invoke(app, ["rank-tls", "--help"])
    assert result.exit_code == 0
    assert "INPUT_DIR" in result.output
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


def test_run_with_save_skips_headless_check():
    """With --save the headless check is bypassed. The pipeline may still fail
    for other reasons, but it shouldn't exit with headless-abort code 1."""
    import os
    from loguru import logger

    logger.disable("quicklook")
    had_display = os.environ.pop("DISPLAY", None)
    try:
        result = runner.invoke(
            app, ["run", "--name", "TOI-6109", "--save", "--overwrite", "--verbose"]
        )
        assert result.exit_code != 1, f"unexpected headless abort: {result.output}"
    finally:
        if had_display is not None:
            os.environ["DISPLAY"] = had_display
        logger.enable("quicklook")
