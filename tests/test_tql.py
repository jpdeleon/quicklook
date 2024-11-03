import pytest
from quicklook.tql import TessQuickLook
from matplotlib.figure import Figure


@pytest.fixture
def planet_inputs():
    # Provide sample data for testing
    return {
        "target_name": "WASP-21",
        "sector": 56,
        "flux_type": "pdcsap",
        "pipeline": "SPOC",
    }


def test_tql_planet(planet_inputs):
    ql = TessQuickLook(**planet_inputs)
    assert isinstance(ql, TessQuickLook)
    fig = ql.plot_tql()
    assert isinstance(fig, Figure)


@pytest.fixture
def eb_inputs():
    # Provide sample data for testing
    return {
        "target_name": "TIC 144539611",
        "sector": 4,
        "flux_type": "pdcsap",
        "pipeline": "SPOC",
    }


def test_tql_eb(eb_inputs):
    ql = TessQuickLook(**eb_inputs)
    assert isinstance(ql, TessQuickLook)
    fig = ql.plot_tql()
    assert isinstance(fig, Figure)
