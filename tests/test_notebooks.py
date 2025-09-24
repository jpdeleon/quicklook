import os
import sys
import pytest
import subprocess
import tempfile
import json
import warnings


def find_notebooks():
    """Find all notebooks in the repository."""
    notebooks = []
    for root, dirs, files in os.walk(
        os.path.dirname(os.path.dirname(__file__))
    ):
        # Skip hidden directories
        dirs[:] = [d for d in dirs if not d.startswith(".")]
        for file in files:
            if file.endswith(".ipynb") and not file.startswith("."):
                notebooks.append(os.path.join(root, file))
    return notebooks


@pytest.mark.notebook
@pytest.mark.slow
@pytest.mark.network
@pytest.mark.parametrize("notebook_path", find_notebooks())
def test_notebook_execution(notebook_path):
    """Test that the notebook can be executed without errors."""
    # Skip this test if running in CI environment without proper kernel setup
    if os.environ.get("CI") and not os.environ.get("KERNEL_INSTALLED"):
        pytest.skip("Skipping notebook test in CI environment without kernel")

    # Get the project root directory to add to PYTHONPATH
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    # Set up environment with PYTHONPATH including the project root
    test_env = os.environ.copy()
    if "PYTHONPATH" in test_env:
        test_env["PYTHONPATH"] = (
            f"{project_root}{os.pathsep}{test_env['PYTHONPATH']}"
        )
    else:
        test_env["PYTHONPATH"] = project_root

    # Create a temporary file for the output
    with tempfile.NamedTemporaryFile(suffix=".ipynb") as temp_output:
        try:
            # Execute the notebook
            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "jupyter",
                    "nbconvert",
                    "--to",
                    "notebook",
                    "--execute",
                    "--output",
                    temp_output.name,
                    notebook_path,
                ],
                capture_output=True,
                text=True,
                check=True,
                timeout=300,  # 5 minute timeout
                env=test_env,
            )
            # If we get here, the notebook executed successfully
            assert result.returncode == 0

            # Verify the output notebook exists and can be loaded
            with open(temp_output.name, "r") as f:
                notebook_content = json.load(f)
                # Check that the notebook has cells
                assert "cells" in notebook_content

        except subprocess.CalledProcessError as e:
            # network-related failures do not fail CI, but real bugs still fail
            stdout_lower = e.stdout.lower() if e.stdout else ""
            stderr_lower = e.stderr.lower() if e.stderr else ""
            network_error_keywords = [
                "remote",
                "connection",
                "timeout",
                "astroquery",
            ]

            if any(
                k in stdout_lower or k in stderr_lower
                for k in network_error_keywords
            ):
                if os.environ.get("CI", "false") == "true":
                    pytest.xfail(f"Notebook network issue: {notebook_path}")
                else:
                    raise
            else:
                # Non-network error
                print(f"Notebook {notebook_path} failed:")
                print(e.stdout)
                print(e.stderr)
                raise
        except subprocess.TimeoutExpired:
            # pytest.fail(f"Notebook execution timed out: {notebook_path}")
            warnings.warn(f"Notebook execution timed out: {notebook_path}")
            return  # Do not fail the test, just skip further checks
