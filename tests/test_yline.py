import lightkurve as lk
import numpy as np
import pytest
from quicklook.plot import plot_odd_even_transit, plot_secondary_eclipse


class AttrDict(dict):
    def __getattr__(self, item):
        try:
            return self[item]
        except KeyError:
            raise AttributeError(item)


@pytest.fixture
def dummy_lc():
    # plot_odd_even_transit/plot_secondary_eclipse expect a lightkurve
    # LightCurve (attribute-style column access, .fold(), time_original
    # tracking), not a plain astropy TimeSeries.
    time = np.linspace(0, 10, 100)
    flux = np.ones(100)
    return lk.LightCurve(time=time, flux=flux)


@pytest.mark.parametrize(
    "depth_input,expected_yline",
    [
        (0.005, 0.995),  # depth fraction
        (0.995, 0.995),  # depth below 1 (relative flux)
        (5.0, 0.995),  # depth in ppt
        ((0.005, 0.001), 0.995),  # tuple of depth & err
        ((0.995, 0.001), 0.995),  # tuple of rel flux & err
        (np.array([0.005, 0.001]), 0.995),  # array
    ],
)
def test_yline_odd_even(dummy_lc, depth_input, expected_yline):
    tls_results = AttrDict(
        {
            "depth": depth_input,
            "duration": 0.1,
            "period": 2.0,
            "T0": 0.0,
            "model_folded_phase": np.linspace(0, 1, 50),
            "model_folded_model": np.ones(50),
        }
    )
    fold_lc = dummy_lc.fold(period=2.0, epoch_time=0.0)
    ax = plot_odd_even_transit(fold_lc, tls_results)
    # find horizontal lines
    hlines = [
        line
        for line in ax.get_lines()
        if len(line.get_xdata()) == 2 and line.get_ydata()[0] == line.get_ydata()[1]
    ]
    assert len(hlines) > 0
    y_values = [line.get_ydata()[0] for line in hlines]
    assert np.isclose(y_values[0], expected_yline)


@pytest.mark.parametrize(
    "depth_input,expected_yline",
    [
        (0.005, 0.995),
        (0.995, 0.995),
        (5.0, 0.995),
    ],
)
def test_yline_secondary(dummy_lc, depth_input, expected_yline):
    tls_results = AttrDict(
        {
            "depth": depth_input,
            "duration": 0.1,
            "period": 2.0,
            "T0": 0.0,
        }
    )
    tmask = np.zeros(len(dummy_lc), dtype=bool)
    ax = plot_secondary_eclipse(dummy_lc, tls_results, tmask)
    hlines = [
        line
        for line in ax.get_lines()
        if len(line.get_xdata()) == 2 and line.get_ydata()[0] == line.get_ydata()[1]
    ]
    assert len(hlines) > 0
    y_values = [line.get_ydata()[0] for line in hlines]
    assert np.isclose(y_values[0], expected_yline)
