import pytest
from quicklook.tql import TessQuickLook
from matplotlib.figure import Figure


@pytest.fixture
def sample_inputs():
    # Provide sample data for testing
    return {
        "target_name": "WASP-21",
        "sector": 56,
        "flux_type": "pdcsap",
        "pipeline": "SPOC",
    }


def test_tql(sample_inputs):
    ql = TessQuickLook(**sample_inputs)
    assert isinstance(ql, TessQuickLook)
    fig = ql.plot_tql()
    assert isinstance(fig, Figure)
