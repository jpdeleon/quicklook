"""mast_query retries transient MAST timeouts (seen under parallel batch runs)."""

import pytest
import requests

from quicklook import tglc


class _Resp:
    headers = {"h": "v"}
    content = b'{"ok": 1}'


@pytest.mark.unit
def test_mast_query_retries_read_timeouts(monkeypatch):
    calls = []

    def fake_post(*a, **k):
        calls.append(k["timeout"])
        if len(calls) < 3:
            raise requests.exceptions.ReadTimeout("slow")
        return _Resp()

    monkeypatch.setattr(tglc.requests, "post", fake_post)
    monkeypatch.setattr(tglc.time, "sleep", lambda s: None)

    head, content = tglc.mast_query({"service": "x"})

    assert content == '{"ok": 1}'
    assert len(calls) == 3


@pytest.mark.unit
def test_mast_query_gives_up_after_max_attempts(monkeypatch):
    def fake_post(*a, **k):
        raise requests.exceptions.ConnectionError("down")

    monkeypatch.setattr(tglc.requests, "post", fake_post)
    monkeypatch.setattr(tglc.time, "sleep", lambda s: None)

    with pytest.raises(requests.exceptions.ConnectionError):
        tglc.mast_query({"service": "x"})
