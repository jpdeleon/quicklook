"""Panels for single-event (mono-transit) results in the quicklook figure."""

from __future__ import annotations

import numpy as np

from quicklook.single_event import (
    SingleEvent,
    SingleEventResult,
    local_models,
    min_period,
    period_is_allowed,
)

EVENT_COLOR = "C3"
EXPLAINED_COLOR = "0.5"
_ZOOM_BIN_FRACTION = 0.1  # zoom-panel bin width, in transit durations
_MIN_BIN = 5 / 1440  # days


def plot_single_event_statistic(
    result: SingleEventResult,
    ax,
    unexplained: tuple[SingleEvent, ...] = (),
):
    """Panel 5 (mono layout): detection SNR vs trial mid-time, threshold and events."""
    ax.plot(result.grid_t0, result.grid_snr, lw=0.7, color="k")
    ax.axhline(
        result.snr_min, ls="--", color="r", lw=1, label=f"threshold (SNR={result.snr_min:g})"
    )
    keep = {e.t0 for e in unexplained}
    for ev in result.events:
        color = EVENT_COLOR if ev.t0 in keep else EXPLAINED_COLOR
        ax.plot(ev.t0, ev.snr, "v", color=color, ms=9, mec="k", zorder=5)
    if result.events:
        ax.plot([], [], "v", color=EVENT_COLOR, mec="k", label="single event")
    ax.set_xlabel("Time [BTJD]")
    ax.set_ylabel("Single-event SNR")
    ax.legend(loc="upper right", fontsize="small")
    return ax


def _binned(t, y, width):
    edges = np.arange(t.min(), t.max() + width, width)
    idx = np.digitize(t, edges)
    centers, means = [], []
    for i in np.unique(idx):
        sel = idx == i
        if sel.sum() >= 2:
            centers.append(t[sel].mean())
            means.append(y[sel].mean())
    return np.array(centers), np.array(means)


def plot_event_zoom(time, flux, event: SingleEvent, ax, color: str = EVENT_COLOR, label=None):
    """Panel 7 (mono layout): the event with its transit and best systematics model."""
    m = local_models(time, flux, event.t0, event.duration)
    hours = (m.time - event.t0) * 24
    ax.plot(hours, m.flux, ".", color="0.75", ms=2, zorder=1)
    width = max(event.duration * _ZOOM_BIN_FRACTION, _MIN_BIN) * 24
    bx, by = _binned(hours, m.flux, width)
    ax.plot(bx, by, "o", color="k", ms=3, zorder=2)
    dbic = m.bic[m.systematics_kind] - m.bic["transit"]
    ax.plot(hours, m.transit, color=color, lw=2, label=label or "transit", zorder=3)
    ax.plot(
        hours,
        m.systematics,
        color="C0",
        ls="--",
        lw=1.5,
        label=rf"{m.systematics_kind} ($\Delta$BIC={dbic:.0f})",
        zorder=3,
    )
    ax.set_xlabel(f"Hours from T0={event.t0:.3f} BTJD")
    ax.set_ylabel("Normalized flux")
    ax.legend(loc="lower left", fontsize="small")
    return ax


def plot_period_constraint(time, event: SingleEvent, ax):
    """Panel 8 (mono layout): which periods the data allow for a single event."""
    time = np.sort(np.asarray(time, dtype=float))
    span = time[-1] - time[0]
    periods = np.arange(event.duration, 2 * span, event.duration / 4)
    allowed = np.array([period_is_allowed(time, event.t0, event.duration, p) for p in periods])
    ax.fill_between(periods, 0, allowed, step="mid", color="C2", alpha=0.5, label="allowed")
    p_min = min_period(time, event.t0, event.duration)
    ax.axvline(p_min, color="k", ls="--", label=f"P min = {p_min:.1f} d")
    ax.set_ylim(0, 1.15)
    ax.set_yticks([])
    ax.set_xlim(0, periods[-1])
    ax.set_xlabel("Period [d]")
    ax.set_title("Periods not excluded by the data", fontsize="medium")
    ax.legend(loc="upper left", fontsize="small")
    return ax


def mark_events(ax, events, color: str = EVENT_COLOR):
    """Shade each event's duration on a light-curve panel."""
    for ev in events:
        ax.axvspan(ev.t0 - ev.duration / 2, ev.t0 + ev.duration / 2, color=color, alpha=0.25)
    return ax
