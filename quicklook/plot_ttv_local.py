"""Byproduct figures of the local-baseline TTV diagnostics (``*_ttv.png``, ``*_ttv_transits.png``)."""

from __future__ import annotations

import matplotlib.pyplot as pl
import numpy as np

from quicklook.ttv import trapezoid
from quicklook.ttv_local import BASELINE_ORDER, HALF_WINDOW, LocalTTVResult

MIN = 1440.0
GALLERY_COLS = 6
_SOURCE_LABEL = {"tls": "TLS", "chain": "chained single events", "toi": "TOI", "custom": "custom"}


def plot_oc(result: LocalTTVResult, ax) -> None:
    """Single-transit (dots; x = rejected) and visit-mean (open squares) offsets."""
    for s in result.singles:
        ok = s.fit.usable
        ax.errorbar(
            s.epoch,
            s.fit.offset * MIN,
            s.err * MIN,
            fmt="o" if ok else "x",
            color="C0" if ok else "0.7",
            ms=4,
            alpha=0.85 if ok else 0.5,
            capsize=2,
        )
    for v in result.usable_visits():
        ax.errorbar(
            np.mean(v.epochs),
            v.fit.offset * MIN,
            v.err * MIN,
            fmt="s",
            mfc="none",
            mec="k",
            ecolor="k",
            ms=10,
            capsize=4,
        )
    ax.axhline(0, color="k", ls=":", lw=0.8)
    title = f"O$-$C vs {_SOURCE_LABEL.get(result.source, result.source)} ephemeris (P={result.period:.5f} d)"
    if np.isfinite(result.visit_pvalue):
        title += f"\nvisit means: $\\chi^2$={result.visit_chi2:.1f}/{result.visit_dof} dof, p={result.visit_pvalue:.1e}"
    ax.set_title(title, fontsize="medium")
    ax.set_xlabel("Epoch")
    ax.set_ylabel("O$-$C [min]")


def plot_duration_offset_map(result: LocalTTVResult, ax) -> None:
    im = ax.pcolormesh(
        result.grid_offsets * MIN,
        result.grid_durations * 24,
        result.grid_snr,
        shading="auto",
        cmap="magma",
    )
    pl.colorbar(im, ax=ax, pad=0.01).set_label(r"$\sqrt{\Delta\chi^2}$")
    i, j = np.unravel_index(np.nanargmax(result.grid_snr), result.grid_snr.shape)
    ax.plot(result.grid_offsets[j] * MIN, result.grid_durations[i] * 24, "c+", ms=14, mew=2)
    ax.set_xlabel("Common offset [min]")
    ax.set_ylabel("Duration [h]")
    ax.set_title(f"All transits stacked: best {result.duration * 24:.2f} h", fontsize="medium")


def plot_visit_scans(result: LocalTTVResult, ax) -> None:
    for k, v in enumerate(result.visits):
        if v.fit is None:
            continue
        lbl = f"visit {v.visit + 1} ({len(v.epochs)} tr): {v.fit.offset * MIN:+.0f}$\\pm${v.err * MIN:.0f} min"
        ax.plot(v.fit.offsets * MIN, v.fit.snr_curve, color=f"C{k % 10}", lw=1.5, label=lbl)
    ax.axvline(0, color="0.5", lw=0.8)
    ax.set_xlabel("Common offset per visit [min]")
    ax.set_ylabel("SNR")
    ax.set_title("Visit-mean offset scans", fontsize="medium")
    ax.legend(fontsize="x-small")


def plot_ttv_diagnostics(result: LocalTTVResult, title: str = "") -> pl.Figure:
    fig = pl.figure(figsize=(14, 9), constrained_layout=True)
    gs = fig.add_gridspec(2, 2)
    plot_oc(result, fig.add_subplot(gs[0, :]))
    plot_duration_offset_map(result, fig.add_subplot(gs[1, 0]))
    plot_visit_scans(result, fig.add_subplot(gs[1, 1]))
    cal = [c for c in result.calibration.values() if np.isfinite(c.single_sigma)]
    sub = "raw flux, cubic baseline per transit"
    if cal:
        sub += (
            "; injection-calibrated single-transit σ per visit: "
            + ", ".join(f"{c.single_sigma * MIN:.0f}" for c in cal)
            + " min"
        )
    fig.suptitle(f"{title} TTV diagnostics\n{sub}", fontsize="large")
    return fig


def plot_transit_gallery(result: LocalTTVResult, time, flux, title: str = "") -> pl.Figure:
    """Each transit after removing its own cubic baseline, with the fitted trapezoid."""
    t, f = np.asarray(time, dtype=float), np.asarray(flux, dtype=float)
    n = max(len(result.singles), 1)
    rows = int(np.ceil(n / GALLERY_COLS))
    fig, axes = pl.subplots(
        rows, GALLERY_COLS, figsize=(2.4 * GALLERY_COLS, 2.1 * rows), squeeze=False
    )
    for ax in axes.flat:
        ax.axis("off")
    for ax, s in zip(axes.flat, result.singles):
        ax.axis("on")
        tc0 = result.t0 + s.epoch * result.period
        sel = np.abs(t - tc0) < HALF_WINDOW
        if sel.sum() < BASELINE_ORDER + 3:
            continue
        x = (t[sel] - tc0) / HALF_WINDOW
        model = -trapezoid(t[sel], tc0 + s.fit.offset, result.duration)
        A = np.column_stack([np.vander(x, BASELINE_ORDER + 1), model])
        coef, *_ = np.linalg.lstsq(A, f[sel], rcond=None)
        detr = f[sel] - A[:, :-1] @ coef[:-1] + 1
        hours = (t[sel] - tc0) * 24
        ax.plot(hours, detr, ".", color="0.6", ms=1)
        ax.plot(hours, 1 + model * coef[-1], color="C0" if s.fit.usable else "0.4", lw=1.2)
        ax.axvline(0, color="0.5", lw=0.6)
        ax.set_title(
            f"E{s.epoch}: {s.fit.offset * MIN:+.0f} min, SNR {s.fit.snr:.1f}",
            fontsize=7,
            color="k" if s.fit.usable else "0.6",
        )
        ax.tick_params(labelsize=6)
    fig.suptitle(f"{title}  each transit on its own cubic baseline (grey = not used)", fontsize=10)
    fig.tight_layout()
    return fig


__all__ = [
    "plot_duration_offset_map",
    "plot_oc",
    "plot_transit_gallery",
    "plot_ttv_diagnostics",
    "plot_visit_scans",
]
