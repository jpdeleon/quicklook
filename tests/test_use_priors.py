"""Offline tests for the --use_priors flag.

Verifies that ExoFOP stellar parameters are forwarded to
``transitleastsquares.power(...)`` only when ``use_star_priors=True``,
and that missing/non-finite parameters are silently skipped.
"""

import sys
import types
import pytest


def _install_fake_logger(monkeypatch):
    """Quiet loguru import inside tql._stellar_prior_kwargs."""
    # tql.py imports `from loguru import logger`; nothing to stub here, but
    # the helper does logger.warning/info — harmless when loguru is real.
    return None


class _RecordingTLS:
    """Stand-in for ``transitleastsquares.transitleastsquares``.

    Captures the kwargs passed to ``.power(...)`` so the test can assert
    on them, and returns a dummy results object with the few attributes
    that ``run_tls`` reads downstream (none, in this isolated test).
    """

    last_power_kwargs = None

    def __init__(self, *args, **kwargs):
        pass

    def power(self, **kwargs):
        _RecordingTLS.last_power_kwargs = dict(kwargs)
        raise _StopAfterPower("captured power() kwargs")


class _StopAfterPower(Exception):
    """Sentinel to halt run_tls right after power() is invoked."""


class _AttrDict(dict):
    """Dict that also exposes keys as attributes, mimicking TLS's Bunch."""

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError as e:
            raise AttributeError(name) from e


def _build_qlook(monkeypatch, use_star_priors, star_params):
    """Build a minimally initialized TessQuickLook for run_tls testing.

    Bypasses __init__ entirely — we never need the network, the ExoFOP
    fetch, or a real light curve. We only need the attributes that
    ``run_tls`` and ``_stellar_prior_kwargs`` touch.
    """
    import numpy as np
    from quicklook import tql as tql_mod

    monkeypatch.setattr(tql_mod, "tls", _RecordingTLS)
    monkeypatch.setattr(
        tql_mod,
        "get_params_from_exofop",
        lambda *a, **kw: dict(star_params),
    )

    qlook = tql_mod.TessQuickLook.__new__(tql_mod.TessQuickLook)
    qlook.exofop_data = {"stub": True}
    qlook.verbose = False
    qlook.use_star_priors = use_star_priors
    qlook.tls_use_threads = None
    qlook.Porb_min = 0.5
    qlook.Porb_max = 10.0

    # Minimal flat_lc stub with the attributes run_tls reads.
    class _Col:
        def __init__(self, arr):
            self.value = arr

    class _FlatLC:
        time = _Col(np.linspace(0, 27, 100))
        flux = _Col(np.ones(100))
        flux_err = _Col(np.full(100, 0.001))

    qlook.flat_lc = _FlatLC()
    return qlook


def test_use_priors_off_passes_no_stellar_kwargs(monkeypatch):
    """Default behavior: no stellar kwargs reach TLS.power()."""
    qlook = _build_qlook(
        monkeypatch,
        use_star_priors=False,
        star_params={"srad": 1.2, "srad_e": 0.05, "mass": 1.1, "mass_e": 0.04},
    )
    with pytest.raises(_StopAfterPower):
        qlook.run_tls()

    kw = _RecordingTLS.last_power_kwargs
    assert kw is not None
    assert "R_star" not in kw and "M_star" not in kw
    assert "R_star_min" not in kw and "M_star_max" not in kw
    # Period bounds still forwarded.
    assert kw["period_min"] == 0.5
    assert kw["period_max"] == 10.0


def test_use_priors_on_forwards_stellar_kwargs(monkeypatch):
    """With the flag on, R_star/M_star + ±1σ bounds reach TLS.power()."""
    qlook = _build_qlook(
        monkeypatch,
        use_star_priors=True,
        star_params={"srad": 1.2, "srad_e": 0.05, "mass": 1.1, "mass_e": 0.04},
    )
    with pytest.raises(_StopAfterPower):
        qlook.run_tls()

    kw = _RecordingTLS.last_power_kwargs
    assert kw["R_star"] == pytest.approx(1.2)
    assert kw["R_star_min"] == pytest.approx(1.15)
    assert kw["R_star_max"] == pytest.approx(1.25)
    assert kw["M_star"] == pytest.approx(1.1)
    assert kw["M_star_min"] == pytest.approx(1.06)
    assert kw["M_star_max"] == pytest.approx(1.14)


def test_use_priors_clips_R_star_min_to_tls_floor(monkeypatch):
    """R_star_min must not go below TLS's 0.13 R_sun lower bound."""
    qlook = _build_qlook(
        monkeypatch,
        use_star_priors=True,
        star_params={"srad": 0.15, "srad_e": 0.10, "mass": 0.5, "mass_e": 0.5},
    )
    with pytest.raises(_StopAfterPower):
        qlook.run_tls()

    kw = _RecordingTLS.last_power_kwargs
    assert kw["R_star_min"] == pytest.approx(0.13)
    assert kw["M_star_min"] == pytest.approx(0.1)


def test_use_priors_skips_missing_values(monkeypatch):
    """A non-finite value for one param must not block the others."""
    qlook = _build_qlook(
        monkeypatch,
        use_star_priors=True,
        star_params={
            "srad": float("nan"),
            "srad_e": 0.05,
            "mass": 1.0,
            "mass_e": None,  # no error -> central value only, no bounds
        },
    )
    with pytest.raises(_StopAfterPower):
        qlook.run_tls()

    kw = _RecordingTLS.last_power_kwargs
    assert "R_star" not in kw  # nan srad dropped
    assert "R_star_min" not in kw
    assert kw["M_star"] == pytest.approx(1.0)
    assert "M_star_min" not in kw  # no usable error
    assert "M_star_max" not in kw


def test_use_priors_handles_missing_exofop(monkeypatch):
    """If ExoFOP lookup raises, run_tls still runs without stellar kwargs."""
    import numpy as np
    from quicklook import tql as tql_mod

    monkeypatch.setattr(tql_mod, "tls", _RecordingTLS)

    def _raise(*a, **kw):
        raise KeyError("stellar_parameters")

    monkeypatch.setattr(tql_mod, "get_params_from_exofop", _raise)

    qlook = tql_mod.TessQuickLook.__new__(tql_mod.TessQuickLook)
    qlook.exofop_data = {}
    qlook.verbose = False
    qlook.use_star_priors = True
    qlook.tls_use_threads = None
    qlook.Porb_min = 0.5
    qlook.Porb_max = 10.0

    class _Col:
        def __init__(self, a):
            self.value = a

    class _FlatLC:
        time = _Col(np.linspace(0, 27, 100))
        flux = _Col(np.ones(100))
        flux_err = _Col(np.full(100, 0.001))

    qlook.flat_lc = _FlatLC()

    with pytest.raises(_StopAfterPower):
        qlook.run_tls()

    kw = _RecordingTLS.last_power_kwargs
    assert "R_star" not in kw and "M_star" not in kw


def test_cli_use_priors_flag_wired():
    """Regression: --use_priors flag and TessQuickLook plumbing must be live
    (the flag was previously commented-out in cli/ql.py)."""
    import quicklook.cli.ql as ql_mod

    src = open(ql_mod.__file__).read()
    assert '"--use_priors"' in src
    assert "use_star_priors=args.use_priors" in src


def test_gui_use_priors_checkbox_wired():
    """Regression: GUI checkbox + Flask plumbing for use_priors must be live.

    Three things must agree: the HTML checkbox exists and defaults to checked,
    the form parser maps it to a ``use_priors`` kwarg, and the Flask worker
    forwards it as ``use_star_priors`` to TessQuickLook.
    """
    import os

    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    html = open(os.path.join(here, "quicklook/app/templates/index.html")).read()
    app_py = open(os.path.join(here, "quicklook/app/app.py")).read()

    assert 'name="use_priors"' in html
    assert 'id="use_priors" checked' in html  # default-on
    assert "'use_priors'" in html  # in PERSIST_CHECKS
    assert '"use_priors": _is_truthy(form.get("use_priors"))' in app_py
    assert 'use_star_priors=kwargs.get("use_priors"' in app_py


def test_h5_persists_use_star_priors_flag(monkeypatch):
    """append_tls_results must record use_star_priors in the dict that gets
    serialized to .h5, plus the prior kwargs used when the flag is on."""
    import numpy as np
    from quicklook import tql as tql_mod

    monkeypatch.setattr(tql_mod, "tls", _RecordingTLS)
    monkeypatch.setattr(
        tql_mod,
        "get_params_from_exofop",
        lambda *a, **kw: {"srad": 1.0, "srad_e": 0.1, "mass": 0.9, "mass_e": 0.05},
    )

    qlook = tql_mod.TessQuickLook.__new__(tql_mod.TessQuickLook)
    qlook.exofop_data = {"stub": True}
    qlook.verbose = False
    qlook.use_star_priors = True
    qlook.tls_results = _AttrDict({"depth": 0.999})  # supports .depth like TLS's Bunch
    # Minimal attrs touched by append_tls_results.
    qlook.pipeline = "spoc"
    qlook.flux_type = "pdcsap"
    qlook.exptime = 120
    qlook.cadence = "short"
    qlook.Porb_min = 0.5
    qlook.Porb_max = 10.0
    qlook.toi_period = None
    qlook.toi_epoch = None
    qlook.toi_dur = None
    qlook.toi_depth = None
    qlook.gaiaid = 1
    qlook.ticid = 2
    qlook.toiid = None
    qlook.sector = 1
    qlook.flatten_method = "biweight"
    qlook.window_length = 0.5
    qlook.gls = None
    qlook.simbad_obj_type = None

    class _Col:
        def __init__(self, a):
            self.value = a

    class _LC:
        time = _Col(np.array([0.0]))
        flux = _Col(np.array([1.0]))
        flux_err = _Col(np.array([0.001]))

    qlook.raw_lc = _LC()
    qlook.flat_lc = _LC()

    # Stub get_toi_radius (called inside append_tls_results).
    monkeypatch.setattr(tql_mod.TessQuickLook, "get_toi_radius", lambda self: None)

    qlook.append_tls_results()

    assert qlook.tls_results["use_star_priors"] is True
    assert "star_prior_kwargs" in qlook.tls_results
    sp = qlook.tls_results["star_prior_kwargs"]
    assert sp["R_star"] == pytest.approx(1.0)
    assert sp["M_star"] == pytest.approx(0.9)


def test_use_priors_on_logs_use_priors_true(monkeypatch, capsys):
    """run_tls must log the use_priors=True banner and the resolved prior
    values, irrespective of self.verbose."""
    from loguru import logger as loguru_logger

    sink_lines = []
    sink_id = loguru_logger.add(sink_lines.append, level="INFO")
    try:
        qlook = _build_qlook(
            monkeypatch,
            use_star_priors=True,
            star_params={"srad": 1.2, "srad_e": 0.05, "mass": 1.1, "mass_e": 0.04},
        )
        with pytest.raises(_StopAfterPower):
            qlook.run_tls()
    finally:
        loguru_logger.remove(sink_id)

    text = "".join(sink_lines)
    assert "use_priors=True" in text
    assert "Using ExoFOP stellar priors" in text
    assert "R_star" in text and "M_star" in text


def test_use_priors_off_logs_use_priors_false(monkeypatch):
    """run_tls must announce Sun-like defaults when use_priors is off."""
    from loguru import logger as loguru_logger

    sink_lines = []
    sink_id = loguru_logger.add(sink_lines.append, level="INFO")
    try:
        qlook = _build_qlook(
            monkeypatch,
            use_star_priors=False,
            star_params={},
        )
        with pytest.raises(_StopAfterPower):
            qlook.run_tls()
    finally:
        loguru_logger.remove(sink_id)

    text = "".join(sink_lines)
    assert "use_priors=False" in text
    assert "Sun-like defaults" in text


def test_h5_records_false_when_priors_disabled(monkeypatch):
    """When use_star_priors=False, the h5 must record False and omit prior kwargs."""
    import numpy as np
    from quicklook import tql as tql_mod

    qlook = tql_mod.TessQuickLook.__new__(tql_mod.TessQuickLook)
    qlook.exofop_data = {"stub": True}
    qlook.verbose = False
    qlook.use_star_priors = False
    qlook.tls_results = _AttrDict({"depth": 0.999})
    qlook.pipeline = "spoc"
    qlook.flux_type = "pdcsap"
    qlook.exptime = 120
    qlook.cadence = "short"
    qlook.Porb_min = 0.5
    qlook.Porb_max = 10.0
    qlook.toi_period = None
    qlook.toi_epoch = None
    qlook.toi_dur = None
    qlook.toi_depth = None
    qlook.gaiaid = 1
    qlook.ticid = 2
    qlook.toiid = None
    qlook.sector = 1
    qlook.flatten_method = "biweight"
    qlook.window_length = 0.5
    qlook.gls = None
    qlook.simbad_obj_type = None

    class _Col:
        def __init__(self, a):
            self.value = a

    class _LC:
        time = _Col(np.array([0.0]))
        flux = _Col(np.array([1.0]))
        flux_err = _Col(np.array([0.001]))

    qlook.raw_lc = _LC()
    qlook.flat_lc = _LC()
    monkeypatch.setattr(tql_mod.TessQuickLook, "get_toi_radius", lambda self: None)

    qlook.append_tls_results()

    assert qlook.tls_results["use_star_priors"] is False
    assert "star_prior_kwargs" not in qlook.tls_results
