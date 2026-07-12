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
