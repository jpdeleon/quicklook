"""GLS-based detrending tuning must run during ``__init__``.

``flatten_raw_lc()`` calls ``init_gls()`` to pick LOCoR's rotation period, the
rotation-aware biweight window and spline ``max_splines``. ``init_gls()``
reads ``Porb_min``/``Porb_max``, which ``__init__`` used to set only *after*
flattening, so every one of those tunings raised AttributeError, was caught,
and silently fell back to a default (LOCoR period 1 d). Seen in all LOCoR runs
of the benchmark suite (06_young).
"""

from unittest.mock import patch

import lightkurve as lk
import numpy as np
import pytest

from quicklook.tql import TessQuickLook, _gls_period

P_ROT = 5.0


@pytest.fixture
def spotted_lc():
    rng = np.random.default_rng(0)
    t = np.arange(0, 27, 10 / 1440)
    f = 1.0 + 0.02 * np.sin(2 * np.pi * t / P_ROT) + rng.normal(0, 1e-3, t.size)
    lc = lk.LightCurve(time=t, flux=f, flux_err=np.full_like(f, 1e-3))
    lc.meta = {"MISSION": "TESS", "SECTOR": 1, "AUTHOR": "SPOC", "EXPOSURE": 600}
    return lc


def _run(lc, flatten_method):
    exofop = {
        "basic_info": {"star_names": "Test star, TIC 1, Gaia DR3 1", "tic_id": "1"},
        "coordinates": {"ra": 10.0, "dec": -10.0},
        "planet_parameters": [],
    }

    def fake_get_lc(self, *args, **kwargs):
        self.pipeline = "spoc"
        return lc

    with (
        patch("quicklook.tql.get_exofop_json", return_value=exofop),
        patch.object(TessQuickLook, "get_simbad_obj_type", return_value=None),
        patch.object(TessQuickLook, "get_lc", autospec=True, side_effect=fake_get_lc),
        patch.object(TessQuickLook, "check_output_file_exists", return_value=None),
    ):
        return TessQuickLook(
            target_name="TIC 1",
            sector=1,
            flux_type="pdcsap",
            pipeline="SPOC",
            flatten_method=flatten_method,
            tls_use_threads=1,
            verbose=False,
        )


@pytest.mark.unit
def test_locor_uses_gls_rotation_period(spotted_lc):
    ql = _run(spotted_lc, "locor")

    assert ql.window_length == pytest.approx(P_ROT, rel=0.1)


@pytest.mark.unit
def test_biweight_rotation_aware_window_uses_gls(spotted_lc):
    # No TOI duration, P_rot = 5 d > 2 d -> rotation-aware window of 1 d
    # (the toi_dur-free default would have been 0.5 d).
    ql = _run(spotted_lc, "biweight")

    assert ql.window_length == pytest.approx(1.0)


@pytest.mark.unit
def test_tls_period_limits_unchanged(spotted_lc):
    ql = _run(spotted_lc, "biweight")
    span = ql.flat_lc.time.value.max() - ql.flat_lc.time.value.min()

    assert ql.Porb_min == pytest.approx(0.1)
    assert ql.Porb_max <= span


class _FakeGls:
    def __init__(self, period):
        self.best = {"P": period}


@pytest.mark.unit
def test_gls_period_reads_best_peak():
    assert _gls_period(_FakeGls(4.9)) == pytest.approx(4.9)


@pytest.mark.unit
@pytest.mark.parametrize("bad", [np.nan, np.inf, 0.0, -1.0])
def test_gls_period_rejects_unusable_peak(bad):
    with pytest.raises(ValueError, match="no usable period"):
        _gls_period(_FakeGls(bad))


@pytest.mark.unit
def test_locor_period_ignores_trend_longer_than_baseline_fraction():
    # A strong 60 d trend dominates an unrestricted GLS; LOCoR must still get
    # the 2.9 d rotation period (capped at baseline / LOCOR_MIN_CYCLES).
    rng = np.random.default_rng(1)
    t = np.arange(0, 54, 10 / 1440)
    f = (
        1.0
        + 0.05 * np.sin(2 * np.pi * t / 60.0)
        + 0.01 * np.sin(2 * np.pi * t / 2.9)
        + rng.normal(0, 1e-3, t.size)
    )
    lc = lk.LightCurve(time=t, flux=f, flux_err=np.full_like(f, 1e-3))
    lc.meta = {"MISSION": "TESS", "SECTOR": 1, "AUTHOR": "SPOC", "EXPOSURE": 600}

    ql = _run(lc, "locor")

    assert ql.window_length == pytest.approx(2.9, rel=0.1)


@pytest.mark.unit
def test_locor_falls_back_when_gls_peak_sits_at_the_cap():
    # Only a monotonic ramp: the capped GLS peaks at its edge, which is not a
    # rotation period, so LOCoR keeps its 1 d default.
    rng = np.random.default_rng(2)
    t = np.arange(0, 54, 10 / 1440)
    f = 1.0 + 1e-3 * t + rng.normal(0, 1e-3, t.size)
    lc = lk.LightCurve(time=t, flux=f, flux_err=np.full_like(f, 1e-3))
    lc.meta = {"MISSION": "TESS", "SECTOR": 1, "AUTHOR": "SPOC", "EXPOSURE": 600}

    ql = _run(lc, "locor")

    assert ql.window_length == pytest.approx(1.0)
