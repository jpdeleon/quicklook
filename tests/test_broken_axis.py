"""Broken time axis between TESS sectors (ported from allesfitter2 plot_utils)."""

import os

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as pl
import numpy as np
import pytest

from quicklook.broken_axis import (
    SECTOR_GAP_DAYS,
    broken_time_axes,
    detect_time_gaps,
    finish_broken_axes,
)

CADENCE = 120 / 86400


def _two_sectors():
    return np.concatenate([np.arange(0.0, 27.0, CADENCE), np.arange(700.0, 727.0, CADENCE)])


def test_detect_time_gaps_splits_only_on_large_gaps():
    t = np.array([0.0, 1.0, 2.0, 4.5, 5.0, 40.0, 41.0])
    assert detect_time_gaps(t, 2.0) == [(0.0, 2.0), (4.5, 5.0), (40.0, 41.0)]
    assert detect_time_gaps(t, 10.0) == [(0.0, 5.0), (40.0, 41.0)]


def test_detect_time_gaps_handles_unsorted_and_tiny_inputs():
    assert detect_time_gaps(np.array([40.0, 0.0, 1.0]), 10.0) == [(0.0, 1.0), (40.0, 40.0)]
    assert detect_time_gaps(np.array([3.0]), 10.0) == [(3.0, 3.0)]
    assert detect_time_gaps(np.array([]), 10.0) == [(0.0, 0.0)]


def test_single_sector_keeps_the_original_axes():
    fig, axes = pl.subplots(3, 3)
    t = np.arange(0.0, 27.0, CADENCE)
    t = t[(t < 12.5) | (t > 13.5)]  # orbit gap stays below the threshold
    subs = broken_time_axes(axes[1, 1], t)
    assert subs == [axes[1, 1]]
    pl.close(fig)


def test_two_sectors_replace_the_panel_with_one_axes_per_sector():
    # Minor and right-side ticks on, as in quicklook's "science" style.
    style = {"ytick.minor.visible": True, "ytick.right": True}
    with pl.rc_context(style):
        fig, axes = pl.subplots(3, 3)
        slot = axes[1, 1].get_subplotspec()
        subs = broken_time_axes(axes[1, 1], _two_sectors())
        subs[0].figure.canvas.draw()
    assert len(subs) == 2
    assert axes[1, 1] not in fig.axes  # replaced, not overlaid
    lo0, hi0 = subs[0].get_xlim()
    lo1, hi1 = subs[1].get_xlim()
    assert lo0 < 0.5 and 26.5 < hi0 < 30
    assert 697 < lo1 < 700.5 and hi1 > 726.5
    assert subs[1].get_shared_y_axes().joined(subs[0], subs[1])
    assert not subs[0].spines["right"].get_visible()
    assert not subs[1].spines["left"].get_visible()
    # Minor ticks too (the science style turns them on), or they line the break.
    assert not any(t.tick1line.get_visible() for t in subs[1].yaxis.get_minor_ticks())
    assert not any(t.tick2line.get_visible() for t in subs[0].yaxis.get_minor_ticks())
    for ax in subs:
        assert ax.get_subplotspec().get_topmost_subplotspec() == slot
    pl.close(fig)


def test_finish_moves_xlabel_to_the_slot_and_keeps_one_legend():
    fig, axes = pl.subplots(3, 3)
    subs = broken_time_axes(axes[0, 0], _two_sectors())
    for ax in subs:
        ax.plot([0, 727], [1, 1], label="flux")
        ax.legend()
        ax.set_xlabel("Time [BTJD]")
        ax.set_ylabel("Flux")
    host = finish_broken_axes(subs)
    assert host is not None and host.get_xlabel() == "Time [BTJD]"
    assert all(ax.get_xlabel() == "" for ax in subs)
    assert subs[0].get_ylabel() == "Flux" and subs[1].get_ylabel() == ""
    assert subs[0].get_legend() is not None and subs[1].get_legend() is None
    assert subs[0].get_legend()._loc == 0  # "best": a narrow sub-axes has no safe corner
    pl.close(fig)


def test_finish_is_a_no_op_for_a_single_axes():
    fig, ax = pl.subplots()
    ax.set_xlabel("Time")
    assert finish_broken_axes([ax]) is None
    assert ax.get_xlabel() == "Time"
    pl.close(fig)


def test_threshold_keeps_consecutive_sectors_together():
    assert SECTOR_GAP_DAYS == pytest.approx(10.0)
    t = np.concatenate([np.arange(0.0, 27.0, CADENCE), np.arange(28.0, 55.0, CADENCE)])
    assert len(detect_time_gaps(t, SECTOR_GAP_DAYS)) == 1
