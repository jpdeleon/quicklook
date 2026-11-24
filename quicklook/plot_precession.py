"""Figures for the precession search: per-sector shape parameters against time.

Planet branch: b, T14, (Rp/Rs)^2. EB branch: inclination, T14, eclipse depth.
"""

from __future__ import annotations

import matplotlib.pyplot as pl
import numpy as np

from quicklook.precession import DAYS_PER_YEAR, PrecessionResult
from quicklook.precession_eb import EBPrecessionResult

MINUTES_PER_DAY = 1440.0
_MODEL_STYLES = (("static", "--", "C7"), ("linear", "-", "C3"), ("quadratic", ":", "C2"))


def _years(t_mid, t_ref):
    return (np.asarray(t_mid, dtype=float) - t_ref) / DAYS_PER_YEAR


def _draw_models(ax, result, t_yr):
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
    panels = (
        ("b", "b_err", 1.0, "impact parameter b"),
        ("t14", "t14_err", 24.0, "T14 [h]"),
        ("depth", "depth_err", 1e6, "(Rp/Rs)$^2$ [ppm]"),
    )
    headline = (
        f"db/dt = {result.db_dt:.4f} ± {result.db_dt_err:.4f} /yr "
        f"({result.significance:.1f}σ, ΔBIC={result.delta_bic:.1f})"
    )
    return _plot_panels(result, panels, f"{title}  {headline}")


def plot_eb_precession(result: EBPrecessionResult, title: str = "") -> pl.Figure:
    """Inclination (with the trend fits), primary-eclipse T14 and depth per sector."""
    panels = (
        ("inc", "inc_err", 1.0, "inclination i [deg]"),
        ("t14", "t14_err", 24.0, "T14 [h]"),
        ("depth_data", "depth_data_err", 1e3, "primary depth [ppt]"),
    )
    headline = (
        f"di/dt = {result.di_dt:.3f} ± {result.di_dt_err:.3f} deg/yr "
        f"({result.significance:.1f}σ, ΔBIC={result.delta_bic:.1f})"
    )
    fig = _plot_panels(result, panels, f"{title}  {headline}")
    if result.models is not None:
        fits = result.sector_fits
        t_yr = _years([f.t_mid for f in fits], result.models.t_ref)
        model_depth = np.array([f.depth for f in fits]) * 1e3
        fig.axes[2].plot(t_yr, model_depth, "s", mfc="none", color="C0", label="model")
        fig.axes[2].legend(fontsize=8, loc="best")
    return fig


def _plot_panels(result, panels, headline: str) -> pl.Figure:
    fits = result.sector_fits
    fig, axes = pl.subplots(3, 1, figsize=(7, 8), sharex=True)
    if not fits or result.models is None:
        axes[0].text(0.5, 0.5, result.reason or "no fit", ha="center", transform=axes[0].transAxes)
        return fig
    t_yr = _years([f.t_mid for f in fits], result.models.t_ref)
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
    axes[0].set_title(f"{headline.strip()} → {verdict}", fontsize=9)
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
