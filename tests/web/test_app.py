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
    # loopback-only bind: default host must never reach beyond localhost
    assert captured.get("host", "127.0.0.1") == "127.0.0.1"

    monkeypatch.setenv("QUICKLOOK_DEBUG", "1")
    app_module.main()
    assert captured["debug"] is True


def test_run_gui_forwards_host_port_debug(monkeypatch):
    """run_gui passes host/port/debug straight through to Flask's app.run."""
    import quicklook.app.app as app_module

    captured = {}
    monkeypatch.setattr(app_module.app, "run", lambda **kw: captured.update(kw))
    app_module.run_gui(host="0.0.0.0", port=8080, debug=True)
    assert captured["host"] == "0.0.0.0"
    assert captured["port"] == 8080
    assert captured["debug"] is True


# --- matplotlib figure lifecycle -------------------------------------------
#
# pyplot keeps every figure in a global registry. The CLI exits per target, but
# the GUI worker thread outlives each job, so a figure left open accumulates
# across a batch.


def _drive_one_job(app_module, monkeypatch, tmp_path, plot_fn):
    """Run run_quicklook_background once against a fake pipeline."""
    from threading import Event

    name = "FIG-TEST"
    monkeypatch.setattr(app_module, "LOG_DIR", tmp_path)
    monkeypatch.setattr(app_module, "_save_job_history", lambda **kw: None)

    class FakeQuickLook:
        def __init__(self, **kwargs):
            pass

        plot_tql = staticmethod(plot_fn)

    monkeypatch.setattr(app_module, "TessQuickLook", FakeQuickLook)
    app_module.jobs[name] = {
        "status": "queued",
        "log_file": "",
        "error": "",
        "cancel_event": Event(),
        "submitted_at": 0.0,
        "params": {},
        "step_times": {},
    }
    app_module.run_quicklook_background(name=name, cancel_event=Event(), save=False)
    return app_module.jobs[name]


def test_worker_closes_figure_on_success(monkeypatch, tmp_path):
    import matplotlib.pyplot as pl
    import quicklook.app.app as app_module

    pl.close("all")

    def plot_tql(**kwargs):
        return pl.figure(), "out.png", "out.h5"

    info = _drive_one_job(app_module, monkeypatch, tmp_path, plot_tql)
    assert info["status"] == "done"
    assert pl.get_fignums() == []


def test_worker_closes_figure_when_the_job_raises(monkeypatch, tmp_path):
    """A figure opened before the failure must not survive the job."""
    import matplotlib.pyplot as pl
    import quicklook.app.app as app_module

    pl.close("all")

    def plot_tql(**kwargs):
        pl.figure()
        raise RuntimeError("boom")

    info = _drive_one_job(app_module, monkeypatch, tmp_path, plot_tql)
    assert info["status"] == "error"
    assert pl.get_fignums() == []


def test_worker_log_preserves_interleaved_logger_and_stdout(monkeypatch, tmp_path):
    """Logger and stdout must append without overwriting each other."""
    from loguru import logger
    import quicklook.app.app as app_module

    def plot_tql(**kwargs):
        quicklook_logger = logger.patch(lambda record: record.update(name="quicklook.tql"))
        quicklook_logger.info("logger before stdout")
        app_module._tls_stdout.write("stdout between logger messages\n")
        app_module._tls_stdout.flush()
        quicklook_logger.info("logger after stdout")

    info = _drive_one_job(app_module, monkeypatch, tmp_path, plot_tql)
    log_text = (tmp_path / "FIG-TEST.log").read_text()

    assert info["status"] == "done"
    assert "logger before stdout" in log_text
    assert "stdout between logger messages" in log_text
    assert "logger after stdout" in log_text
    assert log_text.index("logger before stdout") < log_text.index("stdout between logger messages")
    assert log_text.index("stdout between logger messages") < log_text.index("logger after stdout")
