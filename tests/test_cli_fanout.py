"""--each-sector / --each-pipeline must forward every option of a single `ql run`.

Regression: the fan-out helper used to hard-code custom_ephem=None,
window_length=None, sigma_clip_flat=None, sigma_clip_raw, gp_kernel_size,
use_star_priors and archival_survey, so those flags were silently ignored.
"""

from concurrent.futures import ThreadPoolExecutor
from unittest.mock import patch

import pytest
from typer.testing import CliRunner

from quicklook.cli.app import app

OPTIONS = [
    "--custom-ephem", "2459705.80", "0.001", "2.8626", "0.00001", "2.36", "0.01",
    "--window-length", "0.5",
    "--window-t14-multiplier", "5",
    "--sigma-clip-flat", "10", "3",
    "--sigma-clip-raw", "8", "4",
    "--gp-kernel-size", "7",
    "--use-priors",
    "--survey", "poss2ukstu_red",
]  # fmt: skip

EXPECTED = {
    "custom_ephem": (2459705.80, 0.001, 2.8626, 0.00001, 2.36, 0.01),
    "window_length": 0.5,
    "window_t14_multiplier": 5.0,
    "sigma_clip_flat": (10.0, 3.0),
    "sigma_clip_raw": (8.0, 4.0),
    "gp_kernel_size": 7,
    "use_star_priors": True,
    "archival_survey": "poss2ukstu_red",
}


def _invoke(fanout_args):
    # Threads instead of processes so the TessQuickLook mock records the calls.
    with (
        patch("concurrent.futures.ProcessPoolExecutor", ThreadPoolExecutor),
        patch("quicklook.utils.get_available_sectors", return_value=[10, 11]),
        patch("quicklook.utils.get_available_pipelines", return_value=["spoc", "qlp"]),
        patch("quicklook.tql.TessQuickLook") as ql,
    ):
        result = CliRunner().invoke(
            app, ["run", "--name", "TOI-1234", "--save"] + fanout_args + OPTIONS
        )
    assert result.exit_code == 0, result.output
    return ql.call_args_list


@pytest.mark.parametrize(
    "fanout_args, key, values",
    [
        (["--each-sector", "--pipeline", "tglc"], "sector", {10, 11}),
        (["--each-pipeline"], "pipeline", {"spoc", "qlp"}),
    ],
)
def test_fanout_forwards_all_options(fanout_args, key, values):
    calls = _invoke(fanout_args)
    assert {c.kwargs[key] for c in calls} == values
    for call in calls:
        for name, expected in EXPECTED.items():
            assert call.kwargs[name] == expected, name


def test_single_run_and_fanout_pass_identical_options():
    with patch("quicklook.tql.TessQuickLook") as ql:
        CliRunner().invoke(app, ["run", "--name", "TOI-1234", "--sector", "10", "--save"] + OPTIONS)
    single = dict(ql.call_args.kwargs)
    fanned = dict(_invoke(["--each-sector"])[0].kwargs)
    # sector varies by construction; tls_use_threads is deliberately split
    # across --jobs for --each-sector, so it differs from a single run.
    for kw in (single, fanned):
        kw.pop("sector")
        kw.pop("tls_use_threads")
    assert fanned == single
