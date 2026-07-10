import pytest
from quicklook.app import app


@pytest.fixture
def client():
    # app.config["TESTING"] = True
    with app.test_client() as client:
        yield client


@pytest.mark.skip(reason="web app under development")
def test_index_page(client):
    """Test if index page loads"""
    rv = client.get("/")
    assert rv.status_code == 200
    assert b"TESS QuickLook Pipeline" in rv.data


@pytest.mark.skip(reason="web app under development")
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


# --- path traversal --------------------------------------------------------
#
# /submit writes LOG_DIR/<name>.log before the target is validated against the
# archive, so a name carrying `..` or a separator used to create a file outside
# the output tree. These exercise the real routes, not just the sanitizer.

TRAVERSAL_NAMES = [
    "../../../../tmp/pwned_by_quicklook",
    "/etc/passwd",
    "..",
    "subdir/TOI-1234",
    "TOI-1234/../../../etc/passwd",
]


@pytest.fixture
def submit_client():
    from quicklook.app.app import app as flask_app

    with flask_app.test_client() as c:
        yield c


@pytest.mark.parametrize("name", TRAVERSAL_NAMES)
def test_submit_rejects_traversal_names(submit_client, name):
    rv = submit_client.post("/submit", data={"name": name})
    assert rv.status_code == 400
    assert rv.get_json()["ok"] is False


@pytest.mark.parametrize("name", TRAVERSAL_NAMES)
def test_submit_writes_no_file_outside_log_dir(submit_client, name, tmp_path):
    from quicklook.app.app import LOG_DIR

    before = set(LOG_DIR.rglob("*"))
    submit_client.post("/submit", data={"name": name})
    assert set(LOG_DIR.rglob("*")) == before
    assert not (tmp_path.parent / "pwned_by_quicklook.log").exists()


def test_batch_submit_skips_bad_name_without_dropping_the_batch(submit_client):
    """A malicious entry is skipped; it must not abort the whole request."""
    rv = submit_client.post("/batch-submit", json={"targets": ["../evil"], "params": {}})
    assert rv.status_code == 200
    body = rv.get_json()
    assert body["submitted"] == []
    assert body["skipped"] == ["../evil"]


def test_debug_is_off_unless_env_var_is_set(monkeypatch):
    """The Werkzeug debug console must never be the default."""
    import quicklook.app.app as app_module

    monkeypatch.delenv("QUICKLOOK_DEBUG", raising=False)
    captured = {}
    monkeypatch.setattr(app_module.app, "run", lambda **kw: captured.update(kw))
    app_module.main()
    assert captured["debug"] is False
    assert "host" not in captured  # loopback-only bind

    monkeypatch.setenv("QUICKLOOK_DEBUG", "1")
    app_module.main()
    assert captured["debug"] is True
