## Development & Testing

- **Install dev tools**:
```bash
pip install -e ".[dev,notebooks]"
```

- **Run tests**:

```bash
tox
```

- **Run only notebooks tests**:

```bash
pytest -v -m notebook
```

- **Run linting and formatting**:

```bash
tox -e lint
tox -e black
```

- **Benchmark performance** (locally, not in CI):

```python
import os
import pytest
import pytest_benchmark

@pytest.mark.skipif(os.getenv("CI") == "true", reason="Skip performance tests in CI")
@pytest.mark.benchmark
def test_tql_runtime(benchmark):
    def run_ql():
        ql = TessQuickLook(**planet_inputs)
        fig = ql.plot_tql()
    benchmark(run_ql)
```
