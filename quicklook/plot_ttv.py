"""Figure panels for the TTV search (step 5 of #42)."""

from __future__ import annotations

import matplotlib.pyplot as pl
import numpy as np

HOURS_PER_DAY = 24.0
MINUTES_PER_DAY = 1440.0
ALIGNED_HALF_WINDOW = 1.5  # durations shown on each side of the fitted mid-time
_N_BINS = 40


def _binned(x, y, n_bins):
    edges = np.linspace(x.min(), x.max(), n_bins + 1)
    idx = np.clip(np.digitize(x, edges) - 1, 0, n_bins - 1)
    centres = (edges[:-1] + edges[1:]) / 2
    med = np.array(
        [np.nanmedian(y[idx == i]) if np.any(idx == i) else np.nan for i in range(n_bins)]
    )
    return centres, med


def plot_aligned_transits(time, flux, times, ax=None):
    """Every fitted transit on its own mid-time, coloured by epoch, plus the binned median."""
    ax = pl.gca() if ax is None else ax
    t, f = np.asarray(time, dtype=float), np.asarray(flux, dtype=float)
    if not times:
        ax.text(0.5, 0.5, "No transits fitted", ha="center", va="center", transform=ax.transAxes)
        return ax
    half = ALIGNED_HALF_WINDOW * max(tt.duration for tt in times)
    epochs = np.array([tt.epoch for tt in times])
    cmap = pl.get_cmap("viridis")
    norm = pl.Normalize(epochs.min(), max(epochs.max(), epochs.min() + 1))
    xs, ys = [], []
    for tt in times:
        sel = np.abs(t - tt.tc) < half
        x = (t[sel] - tt.tc) * HOURS_PER_DAY
        ax.plot(x, f[sel], ".", ms=2, color=cmap(norm(tt.epoch)), alpha=0.5)
        xs.append(x)
        ys.append(f[sel])
    x_all, y_all = np.concatenate(xs), np.concatenate(ys)
    ok = np.isfinite(y_all)
    if ok.sum() > _N_BINS:
        xc, yc = _binned(x_all[ok], y_all[ok], _N_BINS)
        ax.plot(xc, yc, "o-", color="r", ms=4, lw=1.5, label="binned")
    sm = pl.cm.ScalarMappable(cmap=cmap, norm=norm)
    pl.colorbar(sm, ax=ax, pad=0.01).set_label("Epoch")
    ax.set_xlabel("Hours from fitted mid-transit")
    ax.set_ylabel("Normalized Flux")
    ax.set_title(f"{len(times)} transits aligned on fitted times")
    ax.legend(fontsize="small", loc="lower left")
    return ax


def plot_oc(stats, ax=None, n_alternatives: int = 0):
    """O-C diagram (minutes) against the best linear ephemeris."""
    ax = pl.gca() if ax is None else ax
    if stats is None or stats.ephemeris is None:
        ax.text(
            0.5,
            0.5,
            "Fewer than 3 transits timed",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        return ax
    oc, err = stats.oc * MINUTES_PER_DAY, stats.oc_err * MINUTES_PER_DAY
    color = "C3" if stats.detected else "k"
    ax.errorbar(stats.epochs, oc, yerr=err, fmt="o", color=color, ms=4, capsize=2)
    ax.axhline(0, color="0.5", ls="--", lw=1)
    verdict = "TTVs detected" if stats.detected else "linear ephemeris OK"
    title = (
        f"O$-$C: amp={stats.amplitude * MINUTES_PER_DAY:.1f} min, p={stats.pvalue:.1e} ({verdict})"
    )
    ax.set_title(title, fontsize="medium")
    if n_alternatives:
        ax.text(
            0.02,
            0.95,
            f"{n_alternatives} other cycle count(s) also fit",
            transform=ax.transAxes,
            va="top",
            fontsize="small",
            color="C1",
        )
    ax.set_xlabel("Epoch")
    ax.set_ylabel("O$-$C [min]")
    return ax


__all__ = ["plot_aligned_transits", "plot_oc"]
