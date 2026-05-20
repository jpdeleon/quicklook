## Development & Testing

- **Install dev tools**:

```bash
pip install -e ".[dev,notebooks]"
```

- **Initialize pre-commit**:
```bash
pre-commit install
```

- **Run tests**:

```bash
tox
```

- **Run only notebooks tests**:

```bash
pytest -v -m notebook
```

- **Run linting, typing, and formatting**:

```bash
tox -e lint
tox -e type
tox -e black
```

- **Benchmark performance** (locally, not in CI):

```bash
pytest --benchmark-only
```


```python
import os
import pytest
import pytest_benchmark

@pytest.mark.skipif(os.getenv("CI") == "true", reason="Skip performance tests in CI")
@pytest.mark.benchmark
def test_tql_runtime(benchmark, planet_inputs):
    def run_ql():
        ql = TessQuickLook(**planet_inputs)
        fig = ql.plot_tql()
    benchmark(run_ql)
```

- **Re-install tox**:

```bash
tox -r
```

- **Check minimum Python required by dependencies**:

`check_pyproject_deps.py` reads `pyproject.toml`, looks up each dependency's
`Requires-Python` in the currently installed environment, and prints the
overall minimum Python version needed.

```bash
python check_pyproject_deps.py
```
