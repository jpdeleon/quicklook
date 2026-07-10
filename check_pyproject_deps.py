import sys
import os
from packaging.version import parse as parse_version
from packaging.specifiers import SpecifierSet

# For Python ≥3.11
try:
    import tomllib
except ImportError:
    import tomli as tomllib  # pip install tomli

try:
    from importlib.metadata import distributions
except ImportError:
    from importlib_metadata import distributions  # pip install importlib-metadata


def read_pyproject(file_path="pyproject.toml"):
    with open(file_path, "rb") as f:
        data = tomllib.load(f)
    # Poetry dependencies
    poetry_deps = data.get("tool", {}).get("poetry", {}).get("dependencies", {})
    # PEP 621 dependencies
    pep621_deps = data.get("project", {}).get("dependencies", {})
    # Merge, remove python itself
    deps = {}
    if poetry_deps:
        deps.update({k: v for k, v in poetry_deps.items() if k.lower() != "python"})
    if pep621_deps:
        # pep621_deps can be list like ["requests>=2.0"]
        for item in pep621_deps:
            name = item.split()[0].split(">=")[0].split("==")[0].split("<=")[0]
            deps[name] = item
    return list(deps.keys())


def get_requires_python(pkg_name):
    try:
        dist = next(d for d in distributions() if d.metadata["Name"].lower() == pkg_name.lower())
        requires_python = dist.metadata.get("Requires-Python")
        return requires_python or "Any"
    except StopIteration:
        return "Not installed / metadata missing"


def min_python_from_specifier(spec):
    """Return the lowest Python version that satisfies a specifier string."""
    if spec in (None, "Any"):
        return None
    try:
        spec_set = SpecifierSet(spec)
        # Find the lowest version manually (best-effort)
        for major in range(2, 4):
            for minor in range(0, 20):
                v = f"{major}.{minor}"
                if parse_version(v) in spec_set:
                    return v
    except Exception:
        return None
    return None


def main(pyproject_file="pyproject.toml"):
    packages = read_pyproject(pyproject_file)
    results = []

    for pkg in packages:
        req_python = get_requires_python(pkg)
        min_py = min_python_from_specifier(req_python)
        results.append((pkg, req_python, min_py))

    print(f"{'Package':25} | {'Requires-Python':20} | {'Min Python'}")
    print("-" * 70)
    min_versions = []
    for pkg, spec, min_py in results:
        print(f"{pkg:25} | {spec:20} | {min_py or 'Unknown'}")
        if min_py:
            min_versions.append(parse_version(min_py))

    if min_versions:
        overall_min = str(min(min_versions))
        print("\nOverall minimum Python version to satisfy all dependencies:", overall_min)
    else:
        print("\nCould not determine overall minimum Python version (check metadata).")


if __name__ == "__main__":
    main()
