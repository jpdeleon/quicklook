"""Search inputs in the web GUI: multi-sector runs, the custom ephemeris, result figures.

The GUI used to read Sector as an integer, so "all" silently became -1 (latest
sector) and Precession search always skipped as "single sector"; a malformed
custom ephemeris was dropped silently (TLS ran instead) or failed in the worker.
"""

import os
import time

import pytest

import quicklook.app.app as app_module
from quicklook.app.app import _job_dedup_signature, _parse_params, _validate_params


@pytest.fixture
def client(monkeypatch):
    submitted = []
    monkeypatch.setattr(
        app_module, "_submit_job", lambda name, params: submitted.append((name, params)) or name
    )
    with app_module.app.test_client() as c:
        c.submitted = submitted
        yield c


def _form(**extra):
    return {"name": "TOI-1355", "pipeline": "spoc", **extra}


# --- parsing and validation -------------------------------------------------------


@pytest.mark.parametrize(
    "sector, expected",
    [("-1", -1), ("", -1), ("57", 57), ("all", "all"), ("57,58,77-78", [57, 58, 77, 78])],
)
def test_parse_params_reads_sector_spec(sector, expected):
    assert _parse_params({"sector": sector})["sector"] == expected


def test_parse_params_rejects_bad_sector_instead_of_using_latest():
    with pytest.raises(app_module.ParamError, match="Sector"):
        _parse_params({"sector": "57;58"})


def test_multi_sector_needs_exptime():
    with pytest.raises(app_module.ParamError, match="Exposure Time"):
        _validate_params(_parse_params({"sector": "all"}))
    assert _validate_params(_parse_params({"sector": "all", "exptime": "120"}))


def test_precession_needs_multi_sector():
    with pytest.raises(app_module.ParamError, match="Precession search"):
        _validate_params(_parse_params({"sector": "57", "precession_search": "on"}))
    ok = _parse_params({"sector": "57,58", "exptime": "120", "precession_search": "on"})
    assert _validate_params(ok)["precession_search"] is True


def test_precession_rejected_for_fan_out_runs():
    params = _parse_params({"precession_search": "on"})
    with pytest.raises(app_module.ParamError, match="untick"):
        _validate_params(params, fan_out=True)
    # Fan-out runs ignore the Sector field, so no exposure time is needed.
    assert _validate_params(_parse_params({}), fan_out=True)


def test_dedup_signature_distinguishes_sector_lists():
    a = _parse_params({"sector": "57,58", "exptime": "120"})
    b = _parse_params({"sector": "57,58,77", "exptime": "120"})
    assert _job_dedup_signature(a) != _job_dedup_signature(b)
    assert hash(_job_dedup_signature(a))


# --- routes ---------------------------------------------------------------------------


def test_submit_all_reaches_the_pipeline_as_all(client):
    rv = client.post("/submit", data=_form(sector="all", exptime="120", precession_search="on"))
    assert rv.get_json()["ok"] is True
    (_, params), *_ = client.submitted
    assert params["sector"] == "all" and params["exptime"] == 120


@pytest.mark.parametrize(
    "form, reason",
    [
        (_form(sector="57;58"), "Sector"),
        (_form(sector="all"), "Exposure Time"),
        (_form(sector="-1", precession_search="on"), "Precession search"),
    ],
)
def test_submit_rejects_with_reason(client, form, reason):
    rv = client.post("/submit", data=form)
    assert rv.status_code == 400
    assert reason in rv.get_json()["reason"]
    assert client.submitted == []


def test_each_sector_submit_rejects_precession(client):
    rv = client.post(
        "/each-sector-submit",
        json={"name": "TOI-1355", "sectors": [57, 58], "params": {"precession_search": "on"}},
    )
    assert rv.status_code == 400 and client.submitted == []


def test_each_pipeline_submit_rejects_precession(client):
    rv = client.post(
        "/each-pipeline-submit",
        json={"name": "TOI-1355", "pipelines": ["spoc"], "params": {"precession_search": "on"}},
    )
    assert rv.status_code == 400 and client.submitted == []


def test_worker_passes_sector_list_to_quicklook(monkeypatch, tmp_path):
    from threading import Event

    seen = {}

    class FakeQuickLook:
        def __init__(self, **kwargs):
            seen.update(kwargs)

        def plot_tql(self, **kwargs):
            return None, "out.png", "out.h5"

    monkeypatch.setattr(app_module, "LOG_DIR", tmp_path)
    monkeypatch.setattr(app_module, "_save_job_history", lambda **kw: None)
    monkeypatch.setattr(app_module, "TessQuickLook", FakeQuickLook)
    name = "SECTOR-LIST-TEST"
    app_module.jobs[name] = {
        "status": "queued",
        "log_file": "",
        "error": "",
        "cancel_event": Event(),
        "submitted_at": 0.0,
        "params": {},
        "step_times": {},
    }
    app_module.run_quicklook_background(
        name=name, cancel_event=Event(), save=False, sector=[57, 58, 77], exptime=120
    )
    assert seen["sector"] == [57, 58, 77]


# --- result panel -----------------------------------------------------------------------


def _touch(path, age):
    path.write_bytes(b"png")
    stamp = time.time() - age
    os.utime(path, (stamp, stamp))


def test_results_json_shows_quicklook_figure_with_precession_as_extra(monkeypatch, tmp_path):
    stem = "TOI1355_sall_spoc_pdcsap_120s"
    _touch(tmp_path / f"{stem}_tls.h5", 30)
    _touch(tmp_path / f"{stem}.png", 20)
    _touch(tmp_path / f"{stem}_precession.png", 10)  # newest, but not the main figure
    monkeypatch.setattr(app_module, "OUTPUT_DIR", tmp_path)
    with app_module.app.test_client() as c:
        data = c.get("/results-json/TOI-1355").get_json()
    assert data["image"].endswith(f"{stem}.png")
    assert [e["name"] for e in data["extras"]] == [f"{stem}_precession.png"]


def test_results_json_without_h5_falls_back_to_newest(monkeypatch, tmp_path):
    _touch(tmp_path / "TOI1355_old.png", 20)
    _touch(tmp_path / "TOI1355_new.png", 10)
    monkeypatch.setattr(app_module, "OUTPUT_DIR", tmp_path)
    with app_module.app.test_client() as c:
        data = c.get("/results-json/TOI-1355").get_json()
    assert data["image"].endswith("TOI1355_new.png") and data["extras"] == []


def test_index_sector_field_accepts_text():
    html = app_module.app.test_client().get("/").get_data(as_text=True)
    assert 'type="text" name="sector"' in html
    assert "function parseSectorSpec" in html
    assert 'id="resultExtras"' in html


# --- custom ephemeris (TTV search on active stars, masking) ----------------------------

TOI6109B = "2458791.0538 0.001 5.690476 0.00001 1.75 0.2"


@pytest.mark.parametrize(
    "text", [TOI6109B, TOI6109B.replace(" ", ", "), "  " + TOI6109B.replace(" ", "  ") + " "]
)
def test_custom_ephem_accepts_spaces_or_commas(text):
    assert _parse_params({"custom_ephem": text})["custom_ephem"] == [
        2458791.0538, 0.001, 5.690476, 0.00001, 1.75, 0.2,
    ]  # fmt: skip


@pytest.mark.parametrize("text", [None, "", "   "])
def test_empty_custom_ephem_is_none(text):
    assert _parse_params({"custom_ephem": text})["custom_ephem"] is None


@pytest.mark.parametrize(
    "text",
    [
        "2458791.0538 0.001 5.690476 0.00001 1.75",  # five values: used to be dropped silently
        "2458791.0538 0.001 abc 0.00001 1.75 0.2",
        "2458791.0538 0.001 nan 0.00001 1.75 0.2",
        "2458791.0538 0.001 -5.69 0.00001 1.75 0.2",
    ],
)
def test_bad_custom_ephem_is_rejected(client, text):
    rv = client.post("/submit", data=_form(custom_ephem=text))
    assert rv.status_code == 400
    assert "Custom Ephemeris" in rv.get_json()["reason"]
    assert client.submitted == []


def test_ttv_job_with_ephem_and_sector_list(client):
    form = _form(sector="18,58,85", exptime="120", ttv_search="on", custom_ephem=TOI6109B)
    assert client.post("/submit", data=form).get_json()["ok"] is True
    (_, params), *_ = client.submitted
    assert params["sector"] == [18, 58, 85] and params["ttv_search"] is True
    assert params["custom_ephem"][2] == pytest.approx(5.690476)


def test_duo_single_event_job_with_sector_list(client):
    form = _form(sector="19,48", exptime="120", single_event_search="on")
    assert client.post("/submit", data=form).get_json()["ok"] is True
    (_, params), *_ = client.submitted
    assert params["sector"] == [19, 48] and params["single_event_search"] is True


def test_results_json_lists_ttv_companion_figures(monkeypatch, tmp_path):
    stem = "TOI6109_s18-58+85_spoc_pdcsap_120s"
    _touch(tmp_path / f"{stem}_tls.h5", 40)
    _touch(tmp_path / f"{stem}.png", 30)
    _touch(tmp_path / f"{stem}_ttv.png", 20)
    _touch(tmp_path / f"{stem}_ttv_transits.png", 10)
    monkeypatch.setattr(app_module, "OUTPUT_DIR", tmp_path)
    with app_module.app.test_client() as c:
        data = c.get("/results-json/TOI-6109").get_json()
    assert data["image"].endswith(f"{stem}.png")
    assert [e["name"] for e in data["extras"]] == [f"{stem}_ttv.png", f"{stem}_ttv_transits.png"]


# --- finding a job's outputs --------------------------------------------------------------


@pytest.mark.parametrize(
    "target, token",
    [
        ("TOI-216", "TOI0216"),
        ("TOI-216.01", "TOI0216"),
        ("TOI-1355", "TOI1355"),
        ("TIC 12345678", "TIC12345678"),
        ("V907 Sco", "V907Sco"),
        ("TOI-1355_spoc_2", "TOI1355_spoc_2"),
    ],
)
def test_output_name_token_matches_pipeline_filenames(target, token):
    assert app_module._output_name_token(target) == token


def test_results_json_finds_zero_padded_toi(monkeypatch, tmp_path):
    # TessQuickLook writes TOI-216 as TOI0216; the old *TOI216* search found nothing.
    stem = "TOI0216_s1-6_pdcsap_sc"
    _touch(tmp_path / f"{stem}_tls.h5", 20)
    _touch(tmp_path / f"{stem}.png", 15)
    _touch(tmp_path / f"{stem}_ttv.png", 10)
    monkeypatch.setattr(app_module, "OUTPUT_DIR", tmp_path)
    with app_module.app.test_client() as c:
        data = c.get("/results-json/TOI-216").get_json()
    assert data["image"].endswith(f"{stem}.png")
    assert [e["name"] for e in data["extras"]] == [f"{stem}_ttv.png"]


def test_results_json_prefers_paths_recorded_by_the_worker(monkeypatch, tmp_path):
    # A de-duplicated job name ("TOI-1355_spoc") does not appear in the filename.
    stem = "TOI1355_sall_pdcsap_sc"
    _touch(tmp_path / f"{stem}_tls.h5", 30)
    _touch(tmp_path / f"{stem}.png", 20)
    _touch(tmp_path / f"{stem}_precession.png", 10)
    _touch(tmp_path / "TOI1355_s84_pdcsap_sc.png", 1)  # newer run of another job
    monkeypatch.setattr(app_module, "OUTPUT_DIR", tmp_path)
    monkeypatch.setitem(
        app_module.jobs,
        "TOI-1355_spoc",
        {
            "output_png": str(tmp_path / f"{stem}.png"),
            "output_h5": str(tmp_path / f"{stem}_tls.h5"),
        },
    )
    with app_module.app.test_client() as c:
        data = c.get("/results-json/TOI-1355_spoc").get_json()
    assert data["image"].endswith(f"{stem}.png") and data["h5"].endswith(f"{stem}_tls.h5")
    assert [e["name"] for e in data["extras"]] == [f"{stem}_precession.png"]


def test_worker_records_output_paths(monkeypatch, tmp_path):
    from threading import Event

    class FakeQuickLook:
        def __init__(self, **kwargs):
            pass

        def plot_tql(self, **kwargs):
            return None, "/out/X.png", "/out/X_tls.h5"

    monkeypatch.setattr(app_module, "LOG_DIR", tmp_path)
    monkeypatch.setattr(app_module, "_save_job_history", lambda **kw: None)
    monkeypatch.setattr(app_module, "TessQuickLook", FakeQuickLook)
    name = "RECORD-TEST"
    app_module.jobs[name] = {
        "status": "queued", "log_file": "", "error": "", "cancel_event": Event(),
        "submitted_at": 0.0, "params": {}, "step_times": {},
    }  # fmt: skip
    app_module.run_quicklook_background(name=name, cancel_event=Event(), save=False)
    assert app_module.jobs[name]["output_png"] == "/out/X.png"
    assert app_module.jobs[name]["output_h5"] == "/out/X_tls.h5"
