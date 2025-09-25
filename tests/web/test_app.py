import pytest
from app.app import app


@pytest.fixture
def client():
    app.config["TESTING"] = True
    with app.test_client() as client:
        yield client


def test_index_page(client):
    """Test if index page loads"""
    rv = client.get("/")
    assert rv.status_code == 200
    assert b"TESS QuickLook Pipeline" in rv.data


def test_post_quicklook(client, monkeypatch):
    """Test POST with sanitized TOI target"""
    monkeypatch.setattr(
        "quicklook.tql.TessQuickLook.plot_tql", lambda self, **kw: (None, "out.png", "out.h5")
    )
    rv = client.post(
        "/", data={"name": "TOI-5169", "sector": "-1", "pipeline": "spoc", "fluxtype": "pdcsap"}
    )
    assert rv.status_code == 200
    assert b"out.png" in rv.data
