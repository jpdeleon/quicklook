"""Exit status of scripts/benchmark_tglc.py.

Batch runners (the tql-gallery benchmark suite) need to tell "the comparison ran
and the light curves disagree" apart from "the script crashed", so a
disagreement exits with EXIT_DISAGREE (3) rather than 1.
"""

import importlib.util
from pathlib import Path

import pytest

from quicklook.benchmark_tglc import BenchmarkResult

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "benchmark_tglc.py"


@pytest.fixture
def script(monkeypatch, tmp_path):
    spec = importlib.util.spec_from_file_location("benchmark_tglc_script", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.chdir(tmp_path)
    return module


def _fake_run(passed: bool):
    def run_benchmark(target, sector, **kwargs):
        result = BenchmarkResult(
            target=target,
            sector=sector,
            n_local=100,
            n_mast=100,
            n_matched=100,
            overlap_frac=1.0,
            passed=passed,
            failures=() if passed else ("cal_aper_flux: pearson r 0.5 < 0.9",),
        )
        return result, None, None

    return run_benchmark


def _argv(tmp_path):
    return ["TIC 1", "--sector", "3", "--no-plot", "--outdir", str(tmp_path / "out")]


def test_agreement_exits_zero(script, monkeypatch, tmp_path):
    monkeypatch.setattr(script, "run_benchmark", _fake_run(passed=True))

    assert script.main(_argv(tmp_path)) == 0


def test_disagreement_exits_with_distinct_code(script, monkeypatch, tmp_path):
    monkeypatch.setattr(script, "run_benchmark", _fake_run(passed=False))

    code = script.main(_argv(tmp_path))

    assert code == script.EXIT_DISAGREE == 3


def test_disagreement_still_writes_json(script, monkeypatch, tmp_path):
    monkeypatch.setattr(script, "run_benchmark", _fake_run(passed=False))

    script.main(_argv(tmp_path))

    assert (tmp_path / "out" / "TIC1_s3_benchmark.json").exists()


def test_crash_propagates_instead_of_masquerading_as_disagreement(script, monkeypatch, tmp_path):
    def boom(*args, **kwargs):
        raise RuntimeError("TESScut down")

    monkeypatch.setattr(script, "run_benchmark", boom)

    with pytest.raises(RuntimeError):
        script.main(_argv(tmp_path))
