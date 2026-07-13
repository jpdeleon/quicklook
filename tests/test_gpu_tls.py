import sys
import types

import numpy as np

from quicklook import tql


def test_gpu_tls_is_selected_when_cuda_device_is_visible(monkeypatch):
    sentinel = object()
    monkeypatch.setattr(tql, "_gpu_device_count", lambda: 1)
    monkeypatch.setitem(sys.modules, "gputls", types.SimpleNamespace(gtls=sentinel))

    assert tql._get_gpu_tls() is sentinel


def test_cpu_tls_is_fallback_when_no_cuda_device_is_visible(monkeypatch):
    monkeypatch.setattr(tql, "_gpu_device_count", lambda: 0)

    assert tql._get_gpu_tls() is None


def test_cpu_tls_is_fallback_when_gpu_detection_fails(monkeypatch):
    def fail_detection():
        raise RuntimeError("CUDA driver unavailable")

    monkeypatch.setattr(tql, "_gpu_device_count", fail_detection)

    assert tql._get_gpu_tls() is None


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

    assert adapted.depth == 0.99
    assert adapted["rp_rs"] == 0.1
    assert adapted.period_uncertainty == 1.0
    np.testing.assert_array_equal(adapted.model_folded_phase, [0.0, 0.5])
