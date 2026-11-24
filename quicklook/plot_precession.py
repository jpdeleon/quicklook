"""Figure for the precession search: b, T14 and depth per sector against time."""

from __future__ import annotations

import matplotlib.pyplot as pl
import numpy as np

from quicklook.precession import DAYS_PER_YEAR, PrecessionResult

MINUTES_PER_DAY = 1440.0
_MODEL_STYLES = (("static", "--", "C7"), ("linear", "-", "C3"), ("quadratic", ":", "C2"))


def _years(t_mid, t_ref):
    return (np.asarray(t_mid, dtype=float) - t_ref) / DAYS_PER_YEAR


def _draw_models(ax, result: PrecessionResult, t_yr):
    models = result.models
    grid = np.linspace(t_yr.min(), t_yr.max(), 200)
    for name, ls, color in _MODEL_STYLES:
        trend = getattr(models, name)
        if trend is None:
            continue
        curve = np.polynomial.polynomial.polyval(grid, trend.coeffs)
        ax.plot(grid, curve, ls=ls, color=color, label=f"{name} (BIC={trend.bic:.1f})")


def plot_precession(result: PrecessionResult, title: str = "") -> pl.Figure:
    """Three stacked panels; the b panel carries the static/linear/quadratic fits."""
    fits = result.sector_fits
    fig, axes = pl.subplots(3, 1, figsize=(7, 8), sharex=True)
    if not fits or result.models is None:
        axes[0].text(0.5, 0.5, result.reason or "no fit", ha="center", transform=axes[0].transAxes)
        return fig
    t_yr = _years([f.t_mid for f in fits], result.models.t_ref)
    panels = (
        ("b", "b_err", 1.0, "impact parameter b"),
        ("t14", "t14_err", 24.0, "T14 [h]"),
        ("depth", "depth_err", 1e6, "(Rp/Rs)$^2$ [ppm]"),
    )
    for ax, (attr, err_attr, scale, label) in zip(axes, panels):
        y = np.array([getattr(f, attr) for f in fits]) * scale
        yerr = np.array([getattr(f, err_attr) for f in fits]) * scale
        ax.errorbar(t_yr, y, yerr=yerr, fmt="o", color="k", capsize=3)
        for x, yy, f in zip(t_yr, y, fits):
            ax.annotate(f"S{f.sector}", (x, yy), textcoords="offset points", xytext=(4, 4))
        ax.set_ylabel(label)
    _draw_models(axes[0], result, t_yr)
    axes[0].legend(fontsize=8, loc="best")
    verdict = "candidate" if result.is_candidate else "no candidate"
    axes[0].set_title(
        f"{title}  db/dt = {result.db_dt:.4f} ± {result.db_dt_err:.4f} /yr "
        f"({result.significance:.1f}σ, ΔBIC={result.delta_bic:.1f}) → {verdict}",
        fontsize=9,
    )
    if result.flags:
        axes[-1].text(
            0.01,
            0.02,
            "flags: " + ", ".join(result.flags),
            transform=axes[-1].transAxes,
            fontsize=8,
            color="C3",
        )
    axes[-1].set_xlabel(f"years since BTJD {result.models.t_ref:.1f}")
    fig.tight_layout()
    return fig
