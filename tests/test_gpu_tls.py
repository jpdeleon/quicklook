import sys
import types
from pathlib import Path

import numpy as np

from quicklook import tql


def test_gpu_extra_installs_cuda_toolkit_headers():
    pyproject = (Path(__file__).parents[1] / "pyproject.toml").read_text()

    assert 'gpu = ["gputls", "cupy-cuda12x[ctk]"]' in pyproject


def test_gpu_tls_is_selected_when_cuda_device_is_visible(monkeypatch):
    sentinel = object()
    monkeypatch.setattr(tql, "_gpu_device_count", lambda: 1)
    monkeypatch.setitem(sys.modules, "gputls", types.SimpleNamespace(gtls=sentinel))

    assert tql._get_gpu_tls() is sentinel


def test_cpu_tls_is_fallback_when_no_cuda_device_is_visible(monkeypatch):
    monkeypatch.setattr(tql, "_gpu_device_count", lambda: 0)

    assert tql._get_gpu_tls() is None


def test_cpu_tls_is_fallback_when_gtls_is_not_installed(monkeypatch):
    monkeypatch.setattr(tql, "_gpu_device_count", lambda: 1)
    monkeypatch.setitem(sys.modules, "gputls", None)

    assert tql._get_gpu_tls() is None


def test_cpu_tls_is_fallback_when_gpu_detection_fails(monkeypatch):
    def fail_detection():
        raise RuntimeError("CUDA driver unavailable")

    monkeypatch.setattr(tql, "_gpu_device_count", fail_detection)

    assert tql._get_gpu_tls() is None


def test_cpu_tls_is_fallback_when_gtls_execution_fails(monkeypatch):
    cpu_result = object()

    class BrokenGTLS:
        def __init__(self, *args, **kwargs):
            pass

        def power(self, **kwargs):
            raise RuntimeError("Failed to find CUDA headers")

    class RecordingTLS:
        calls = 0

        def __init__(self, *args, **kwargs):
            RecordingTLS.calls += 1

        def power(self, **kwargs):
            return cpu_result

    monkeypatch.setattr(tql, "_get_gpu_tls", lambda: BrokenGTLS)
    monkeypatch.setattr(tql, "tls", RecordingTLS)

    values = np.linspace(0, 1, 10)
    qlook = tql.TessQuickLook.__new__(tql.TessQuickLook)
    qlook.flat_lc = types.SimpleNamespace(
        time=types.SimpleNamespace(value=values),
        flux=types.SimpleNamespace(value=np.ones(10)),
        flux_err=types.SimpleNamespace(value=np.full(10, 0.01)),
    )
    qlook.Porb_min = 0.5
    qlook.Porb_max = 10.0
    qlook.tls_use_threads = None
    qlook.use_star_priors = False
    qlook.verbose = False

    qlook.run_tls()

    assert qlook.tls_results is cpu_result
    assert RecordingTLS.calls == 1


def test_gtls_result_is_adapted_to_tls_conventions():
    result = types.SimpleNamespace(
        depth=0.01,
        periods=np.array([1.0, 2.0, 3.0]),
        power=np.array([0.1, 1.0, 0.1]),
        period=2.0,
        duration=0.1,
        T0=0.5,
        SDE=8.0,
    )

    class Model:
        def showFit(self):
            return [1.0, 0.99], [0.0, 0.5], [1.0, 0.99]

    adapted = tql._adapt_gtls_result(result, Model())

    assert adapted.depth == 0.01
    assert adapted["rp_rs"] == 0.1
    assert adapted.period_uncertainty == 1.0
    np.testing.assert_array_equal(adapted.model_folded_phase, [0.0, 0.5])
