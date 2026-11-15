"""Cross-target check for single events at the same epoch (issue #42, step 6)."""

import numpy as np
import pandas as pd
import pytest
from typer.testing import CliRunner

from quicklook import h5io
from quicklook.cli.app import app
from quicklook.shared_epochs import (
    flag_shared_events,
    load_events,
    shared_epochs_report,
    vet_targets,
)


def _events(rows):
    return pd.DataFrame(rows, columns=["ticid", "sector", "event_class", "t0"])


@pytest.mark.unit
def test_events_from_three_targets_at_one_epoch_are_shared():
    ev = _events(
        [
            (1, "14", "mono", 1695.90),
            (2, "14", "mono", 1696.00),
            (3, "14", "duo", 1696.10),
            (4, "14", "mono", 1700.00),
        ]
    )

    out = flag_shared_events(ev, tol=0.25, min_targets=3)

    assert out["shared"].tolist() == [True, True, True, False]
    assert out["n_targets_nearby"].tolist() == [3, 3, 3, 1]
    assert "shared" not in ev.columns  # input untouched


@pytest.mark.unit
def test_several_events_of_one_target_do_not_count_as_shared():
    ev = _events([(1, "14", "mono", 1696.0), (1, "14", "mono", 1696.1), (2, "14", "mono", 1696.2)])

    assert not flag_shared_events(ev, tol=0.25, min_targets=3)["shared"].any()


@pytest.mark.unit
def test_tolerance_is_symmetric_and_not_binned():
    # A fixed 0.5 d grid would split 1695.98 and 1696.02 across bins.
    ev = _events(
        [(1, "14", "mono", 1695.98), (2, "14", "mono", 1696.02), (3, "14", "mono", 1696.2)]
    )

    assert flag_shared_events(ev, tol=0.25)["shared"].all()


@pytest.mark.unit
def test_vet_targets_relabels_only_fully_shared_single_event_targets():
    ev = _events(
        [
            (1, "14", "mono", 1696.0),
            (2, "14", "mono", 1696.1),
            (3, "14", "periodic", 1696.05),
            (4, "14", "duo", 1696.0),
            (4, "14", "duo", 1705.0),  # one real-looking event left
        ]
    )

    targets = vet_targets(flag_shared_events(ev)).set_index("ticid")

    assert targets.loc[1, "event_class_vetted"] == "systematic"
    assert targets.loc[2, "event_class_vetted"] == "systematic"
    assert targets.loc[3, "event_class_vetted"] == "periodic"
    assert targets.loc[4, "event_class_vetted"] == "duo"
    assert (targets.loc[4, "n_events"], targets.loc[4, "n_shared"]) == (2, 1)


@pytest.mark.unit
def test_empty_inputs():
    flagged = flag_shared_events(_events([]))

    assert flagged.empty and "shared" in flagged.columns
    assert vet_targets(flagged).empty


def _write_h5(path, ticid, t0s, event_class="mono", sector=14):
    h5io.save(
        path,
        {
            "ticid": ticid,
            "sector": sector,
            "event_class": event_class,
            "event_t0s": np.asarray(t0s, dtype=float),
            "event_snrs": np.full(len(t0s), 9.0),
            "event_depths": np.full(len(t0s), 1e-3),
            "event_durations": np.full(len(t0s), 0.1),
            "flux_flat": np.ones(10),  # large arrays are not read
        },
    )


@pytest.mark.unit
def test_load_events_reads_h5_and_skips_targets_without_events(tmp_path):
    _write_h5(tmp_path / "a_tls.h5", 1, [1696.0, 1700.0])
    _write_h5(tmp_path / "b_tls.h5", 2, [], event_class="none")

    ev = load_events(sorted(tmp_path.glob("*_tls.h5")))

    assert ev["ticid"].tolist() == [1, 1]
    assert ev["t0"].tolist() == [1696.0, 1700.0]
    assert ev["sector"].tolist() == ["14", "14"]
    assert np.allclose(ev["snr"], 9.0)


@pytest.mark.unit
def test_report_raises_without_h5(tmp_path):
    with pytest.raises(FileNotFoundError):
        shared_epochs_report(tmp_path)


@pytest.mark.unit
def test_cli_writes_event_and_target_tables(tmp_path):
    batch = tmp_path / "batch"
    batch.mkdir()
    for tic, t0 in ((1, 1696.0), (2, 1696.1), (3, 1696.2), (4, 1702.0)):
        _write_h5(batch / f"TIC{tic}_s14_tglc_psf_lc_tls.h5", tic, [t0])

    result = CliRunner().invoke(app, ["shared-epochs", str(batch)])

    assert result.exit_code == 0, result.output
    assert "3/4 events shared; 3 mono/duo target(s) re-labelled systematic" in result.output
    targets = pd.read_csv(str(batch) + "_shared_targets.csv").set_index("ticid")
    assert targets.loc[4, "event_class_vetted"] == "mono"
    assert pd.read_csv(str(batch) + "_shared_events.csv")["shared"].sum() == 3


@pytest.mark.unit
def test_cli_rejects_directory_without_h5(tmp_path):
    result = CliRunner().invoke(app, ["shared-epochs", str(tmp_path)])

    assert result.exit_code != 0
