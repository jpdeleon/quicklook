"""Broken time axis for light curves that span separate TESS sectors.

Ported from ``allesfitter/plot_utils.py`` in jpdeleon/allesfitter2
(``detect_time_gaps``, ``broken_xaxis_subplots``, ``_draw_break_marks``), and
adapted to replace an existing panel of a ``pl.subplots`` grid in place:

* ``broken_time_axes(ax, time)`` swaps ``ax`` for one sub-axes per data
  segment (gaps > ``SECTOR_GAP_DAYS``) in the same grid slot, sharing y, with
  inner spines hidden and diagonal break marks. Without such a gap it returns
  ``[ax]`` unchanged, so callers can loop either way.
* ``finish_broken_axes(axes)`` tidies labels after plotting: one y label, one
  legend, and a single x label centred under the whole slot.
"""

from __future__ import annotations

import numpy as np
from matplotlib import gridspec as _gridspec

# Larger than TESS orbit gaps (~1 d) and the gap between consecutive sectors
# (~1-2 d), so only separate sectors are broken apart.
SECTOR_GAP_DAYS = 10.0
BREAK_MARK_SIZE = 0.015
EDGE_PAD_FRAC = 0.01
SUBPLOT_WSPACE = 0.05
LEGEND_FONTSIZE = "x-small"
_HOST_LABELPAD = 18  # points; keeps the shared x label below the tick labels

__all__ = ["SECTOR_GAP_DAYS", "broken_time_axes", "detect_time_gaps", "finish_broken_axes"]


def detect_time_gaps(time, gap_threshold_days):
    """Return ``(start, end)`` segments separated by gaps > ``gap_threshold_days``.

    Robust to unsorted input. Empty or single-point inputs return one segment.
    """
    t = np.asarray(time, dtype=float)
    t = t[np.isfinite(t)]
    if t.size == 0:
        return [(0.0, 0.0)]
    if t.size == 1:
        return [(float(t[0]), float(t[0]))]
    t = np.sort(t)
    gap_idx = np.where(np.diff(t) > gap_threshold_days)[0]
    if gap_idx.size == 0:
        return [(float(t[0]), float(t[-1]))]
    starts = np.concatenate(([0], gap_idx + 1))
    ends = np.concatenate((gap_idx, [t.size - 1]))
    return [(float(t[s]), float(t[e])) for s, e in zip(starts, ends)]


def _draw_break_marks(axes, size):
    """Diagonal "axis-break" marks at each interior break (axes-relative ``size``)."""
    kw = dict(color="k", clip_on=False, lw=1.0)
    n = len(axes)
    for i, ax in enumerate(axes):
        tr = ax.transAxes
        if i < n - 1:
            ax.plot((1 - size, 1 + size), (-size, +size), transform=tr, **kw)
            ax.plot((1 - size, 1 + size), (1 - size, 1 + size), transform=tr, **kw)
        if i > 0:
            ax.plot((-size, +size), (1 - size, 1 + size), transform=tr, **kw)
            ax.plot((-size, +size), (-size, +size), transform=tr, **kw)


def broken_time_axes(ax, time, gap_threshold_days: float = SECTOR_GAP_DAYS) -> list:
    """Replace ``ax`` by one sub-axes per data segment in the same grid slot.

    Sub-axes get equal widths (each sector stays readable however long the gap
    between them), share y, and have per-segment x limits with a small pad.
    Returns ``[ax]`` untouched when no gap exceeds the threshold.
    """
    segments = detect_time_gaps(time, gap_threshold_days)
    if len(segments) == 1:
        return [ax]
    fig, spec = ax.figure, ax.get_subplotspec()
    ax.remove()
    sub_gs = _gridspec.GridSpecFromSubplotSpec(
        1, len(segments), subplot_spec=spec, wspace=SUBPLOT_WSPACE
    )
    first = fig.add_subplot(sub_gs[0])
    axes = [first] + [fig.add_subplot(sub_gs[i], sharey=first) for i in range(1, len(segments))]
    for sub, (lo, hi) in zip(axes, segments):
        pad = max(hi - lo, np.finfo(float).eps) * EDGE_PAD_FRAC
        sub.set_xlim(lo - pad, hi + pad)
    n = len(axes)
    for i, sub in enumerate(axes):
        if i < n - 1:
            sub.spines["right"].set_visible(False)
            sub.tick_params(which="both", right=False, labelright=False)
        if i > 0:
            sub.spines["left"].set_visible(False)
            sub.tick_params(which="both", left=False, labelleft=False)
    _draw_break_marks(axes, BREAK_MARK_SIZE)
    # Plotting full-range data later would reset the limits; keep them fixed.
    for sub in axes:
        sub.set_autoscalex_on(False)
    return axes


def _shrink_legend(ax):
    """Redraw a legend smaller and placed by matplotlib: a sub-axes is narrow,
    so a full-size legend at a fixed corner spills out or hides data."""
    legend = ax.get_legend()
    if legend is None:
        return
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, loc="best", fontsize=LEGEND_FONTSIZE, framealpha=0.8)


def finish_broken_axes(axes):
    """After plotting: one y label and legend, and one x label under the whole slot.

    Returns the label-carrying host axes, or ``None`` for a single axes.
    """
    if len(axes) < 2:
        return None
    xlabel = next((a.get_xlabel() for a in axes if a.get_xlabel()), "")
    for i, sub in enumerate(axes):
        sub.set_xlabel("")
        if i > 0:
            sub.set_ylabel("")
            legend = sub.get_legend()
            if legend is not None:
                legend.remove()
    _shrink_legend(axes[0])
    host = axes[0].figure.add_subplot(axes[0].get_subplotspec().get_topmost_subplotspec())
    host.set_frame_on(False)
    host.tick_params(
        left=False, bottom=False, top=False, right=False, labelleft=False, labelbottom=False
    )
    host.patch.set_alpha(0)
    host.set_zorder(-1)
    host.set_xlabel(xlabel, labelpad=_HOST_LABELPAD)
    return host
