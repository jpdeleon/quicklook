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


def _run(lc, flatten_method, **kwargs):
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
            **kwargs,
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
@pytest.mark.parametrize(
    ("options", "expected_days"),
    [
        ({}, 0.625),
        ({"window_t14_multiplier": 3}, 0.375),
        ({"window_t14_multiplier": 5, "window_length": 0.8}, 0.8),
    ],
)
def test_duration_scaled_window(spotted_lc, options, expected_days):
    ql = _run(spotted_lc, "biweight", custom_ephem=(1, 0.01, 5, 0.01, 3, 0.1), **options)
    assert ql.window_length == pytest.approx(expected_days)


@pytest.mark.unit
@pytest.mark.parametrize("value", [0, -1, float("nan"), float("inf")])
def test_invalid_duration_multiplier(value):
    from quicklook.exceptions import InvalidInputError

    with pytest.raises(InvalidInputError, match="finite and positive"):
        TessQuickLook("TIC 1", window_t14_multiplier=value)


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


@pytest.mark.unit
@pytest.mark.parametrize(
    ("amplitude", "period", "expected"),
    [(0.02, 5.0, 3), (0.02, 5.1, 5), (0.001, 2, 5), (0.02, 8.9, 5), (0.02, np.nan, 5)],
)
def test_adaptive_multiplier(spotted_lc, amplitude, period, expected):
    from quicklook.tql import _adaptive_multiplier

    lc = spotted_lc.copy()
    lc.flux = 1 + amplitude * np.sin(2 * np.pi * lc.time.value / 2)
    with patch("quicklook.tql.Gls", return_value=_FakeGls(period)):
        multiplier, reason = _adaptive_multiplier(lc)
    assert multiplier == expected
    assert reason


@pytest.mark.unit
@pytest.mark.parametrize(
    ("ephem", "options", "method", "days"),
    [
        (None, {}, "notch", 0.5),
        ((1, 0.01, 5, 0.01, 3, 0.1), {}, "biweight", 0.375),
        ((1, 0.01, 5, 0.01, 3, 0.1), {"window_t14_multiplier": 5}, "biweight", 0.625),
        (None, {"window_length": 0.75}, "biweight", 0.75),
        ((1, 0.01, 5, 0.01, 3, 0.1), {"flatten_method": "biweight"}, "biweight", 0.625),
    ],
)
def test_adaptive_constructor(spotted_lc, ephem, options, method, days):
    with patch("quicklook.tql._adaptive_multiplier", return_value=(3, "fast and variable")) as rule:
        ql = _run(
            spotted_lc,
            options.get("flatten_method"),
            adaptive=True,
            custom_ephem=ephem,
            **{k: v for k, v in options.items() if k != "flatten_method"},
        )
    assert ql.flatten_method == method
    assert ql.window_length == pytest.approx(days)
    assert rule.call_count == int(ephem is not None and not options)
    assert ql.adaptive_reason == (
        "explicit detrending option"
        if options
        else "fast and variable"
        if ephem
        else "no usable T14"
    )
