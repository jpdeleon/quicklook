import pytest


@pytest.fixture
def client():
    from quicklook.app.app import app

    with app.test_client() as test_client:
        yield test_client


def test_root_page_keeps_standalone_urls(client):
    rv = client.get("/")
    assert rv.status_code == 200
    assert b'const BASE = ""' in rv.data
    assert b'href="/static/css/quicklook.css"' in rv.data


def test_forwarded_prefix_applies_to_html_assets_and_browser_calls(client):
    rv = client.get("/", headers={"X-Forwarded-Prefix": "/tess-quicklook"})
    assert rv.status_code == 200
    assert b'const BASE = "/tess-quicklook"' in rv.data
    assert b'href="/tess-quicklook/static/css/quicklook.css"' in rv.data
    assert b"fetch(BASE + '/jobs-json')" in rv.data
    assert b"location.host + BASE + '/ws/'" in rv.data
