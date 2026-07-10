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
    assert b"icon: BASE + '/static/favicon.ico'" in rv.data


@pytest.mark.parametrize(
    ("route", "active_page"),
    [
        ("/gallery", b"Gallery"),
        ("/compare", b"Compare"),
        ("/tls-summary", b"Summary"),
    ],
)
def test_secondary_pages_render_with_prefixed_navigation(client, route, active_page):
    rv = client.get(route, headers={"X-Forwarded-Prefix": "/tess-quicklook"})
    assert rv.status_code == 200
    assert b'const BASE = "/tess-quicklook"' in rv.data
    assert b'href="/tess-quicklook/"' in rv.data
    assert b'href="/tess-quicklook/gallery"' in rv.data
    assert b'href="/tess-quicklook/compare"' in rv.data
    assert b'href="/tess-quicklook/tls-summary"' in rv.data
    assert active_page in rv.data


def test_gallery_and_compare_live_job_links_keep_prefix(client):
    for route in ("/gallery", "/compare"):
        rv = client.get(route, headers={"X-Forwarded-Prefix": "/tess-quicklook"})
        assert b"href=\"' + BASE + '/?watch='" in rv.data


def test_summary_browser_calls_keep_prefix(client):
    rv = client.get("/tls-summary", headers={"X-Forwarded-Prefix": "/tess-quicklook"})
    assert b"window.location.replace(BASE + '/tls-summary?dir='" in rv.data
    assert b"const url = BASE + '/list-dirs'" in rv.data
