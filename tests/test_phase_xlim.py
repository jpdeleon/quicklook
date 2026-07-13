import os

import pytest


def test_phase_window_uses_custom_delta():
    from quicklook.plot import _phase_window

    assert _phase_window(0, 0.02, 0.1) == pytest.approx((-0.1, 0.1))
    assert _phase_window(0.5, 0.02, 0.1) == pytest.approx((0.4, 0.6))


def test_phase_window_defaults_to_auto_delta():
    from quicklook.plot import _phase_window

    assert _phase_window(0, 0.02, None) == pytest.approx((-0.02, 0.02))
    assert _phase_window(0.5, 0.02, None) == pytest.approx((0.48, 0.52))


@pytest.mark.parametrize("bad_value", [0, -0.1, 1.1])
def test_phase_xlim_rejects_invalid_delta(bad_value):
    from quicklook.plot import _phase_window

    with pytest.raises(ValueError, match="phase_xlim"):
        _phase_window(0, 0.02, bad_value)


def test_cli_phase_xlim_flag_wired(monkeypatch):
    """The unified Typer option reaches ``TessQuickLook``."""
    from typer.testing import CliRunner

    import quicklook.tql as tql_mod
    from quicklook.cli.app import app

    calls = {}

    class FakeQuickLook:
        def __init__(self, **kwargs):
            calls.update(kwargs)

        def plot_tql(self):
            return None

    monkeypatch.setattr(tql_mod, "TessQuickLook", FakeQuickLook)

    result = CliRunner().invoke(app, ["run", "--name", "TOI-1234", "--phase-xlim", "0.1", "--save"])

    assert result.exit_code == 0, result.output
    assert calls["phase_xlim"] == pytest.approx(0.1)


def test_gui_phase_xlim_input_wired():
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    html = open(os.path.join(here, "quicklook/app/templates/index.html")).read()
    app_py = open(os.path.join(here, "quicklook/app/app.py")).read()

    assert 'name="phase_xlim"' in html
    assert "'phase_xlim'" in html
    assert '"phase_xlim": safe_float(form.get("phase_xlim") or None, None)' in app_py
    assert 'phase_xlim=kwargs.get("phase_xlim")' in app_py
