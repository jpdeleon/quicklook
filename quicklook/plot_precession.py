"""Figures for the precession search: per-sector shape parameters against time.

Planet branch: b, T14, (Rp/Rs)^2. EB branch: inclination, T14, eclipse depth.
``plot_sector_folds`` shows each sector's phase-folded transit (planet branch).
"""

from __future__ import annotations

import matplotlib.pyplot as pl
import numpy as np

from quicklook.precession import (
    DAYS_PER_YEAR,
    DEFAULT_LIMB_DARKENING,
    HALF_WINDOW_DURATIONS,
    MIN_HALF_WINDOW,
    PrecessionResult,
    _light_curve,
    extract_sector_transits,
    import_batman,
    transit_model,
)
from quicklook.precession_eb import EBPrecessionResult

MINUTES_PER_DAY = 1440.0
MAX_PLOTTED_B_ERR = 1.0  # sectors with a larger b uncertainty are not drawn
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
    """Three stacked panels; the b panel carries the static/linear/quadratic fits.

    Sectors with ``b_err > MAX_PLOTTED_B_ERR`` are left out of every panel (and
    named in the figure). The trend fits and the verdict still use all sectors.
    """
    panels = (
        ("b", "b_err", 1.0, "impact parameter b"),
        ("t14", "t14_err", 24.0, "T14 [h]"),
        ("depth", "depth_err", 1e6, "(Rp/Rs)$^2$ [ppm]"),
    )
    headline = (
        f"db/dt = {result.db_dt:.4f} ± {result.db_dt_err:.4f} /yr "
        f"({result.significance:.1f}σ, ΔBIC={result.delta_bic:.1f})"
    )
    fits = result.sector_fits
    shown = tuple(f for f in fits if f.b_err <= MAX_PLOTTED_B_ERR)
    hidden = tuple(f for f in fits if f not in shown)
    return _plot_panels(result, panels, f"{title}  {headline}", shown, hidden)


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
    fig = _plot_panels(result, panels, f"{title}  {headline}", result.sector_fits)
    if result.models is not None:
        fits = result.sector_fits
        t_yr = _years([f.t_mid for f in fits], result.models.t_ref)
        model_depth = np.array([f.depth for f in fits]) * 1e3
        fig.axes[2].plot(t_yr, model_depth, "s", mfc="none", color="C0", label="model")
        fig.axes[2].legend(fontsize=8, loc="best")
    return fig


def _plot_panels(result, panels, headline: str, fits, hidden=()) -> pl.Figure:
    fig, axes = pl.subplots(3, 1, figsize=(7, 8), sharex=True)
    if not result.sector_fits or result.models is None:
        axes[0].text(0.5, 0.5, result.reason or "no fit", ha="center", transform=axes[0].transAxes)
        return fig
    t_yr = _years([f.t_mid for f in fits], result.models.t_ref)
    for ax, (attr, err_attr, scale, label) in zip(axes, panels):
        ax.set_ylabel(label)
        if not fits:
            continue
        y = np.array([getattr(f, attr) for f in fits]) * scale
        yerr = np.array([getattr(f, err_attr) for f in fits]) * scale
        ax.errorbar(t_yr, y, yerr=yerr, fmt="o", color="k", capsize=3)
        for x, yy, f in zip(t_yr, y, fits):
            ax.annotate(f"S{f.sector}", (x, yy), textcoords="offset points", xytext=(4, 4))
    if fits:
        _draw_models(axes[0], result, t_yr)
        axes[0].legend(fontsize=8, loc="best")
    if hidden:
        names = ", ".join(f"S{f.sector}" for f in hidden)
        axes[-1].text(
            0.99,
            0.97,
            f"not shown (b uncertainty > {MAX_PLOTTED_B_ERR:g}): {names}",
            transform=axes[-1].transAxes,
            ha="right",
            va="top",
            fontsize=8,
            color="C7",
        )
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


FOLD_MAX_COLS = 3
FOLD_BIN_POINTS = 8  # raw points per binned marker
FOLD_MODEL_POINTS = 400


def _binned(x, y, n_bins: int):
    """Mean of ``y`` in ``n_bins`` equal-width bins of ``x`` (empty bins dropped)."""
    edges = np.linspace(x.min(), x.max(), n_bins + 1)
    idx = np.clip(np.digitize(x, edges) - 1, 0, n_bins - 1)
    counts = np.bincount(idx, minlength=n_bins)
    keep = counts > 0
    xs = np.bincount(idx, weights=x, minlength=n_bins)[keep] / counts[keep]
    ys = np.bincount(idx, weights=y, minlength=n_bins)[keep] / counts[keep]
    return xs, ys


def _fold_model(fit, a_rs, period, t0_eff, half, ld):
    """Best-fit transit curve of one sector on a fine phase grid (days from mid-transit)."""
    batman = import_batman()
    phase = np.linspace(-half, half, FOLD_MODEL_POINTS)
    params, model = transit_model(batman, t0_eff + phase, fit.cadence, period, t0_eff, ld)
    return phase, _light_curve(params, model, t0_eff, a_rs, 0.0, fit.rp_rs, fit.b)


def plot_sector_folds(
    result: PrecessionResult,
    time,
    flux,
    flux_err,
    spans,
    period: float,
    t0: float,
    duration: float,
    title: str = "",
    limb_darkening=DEFAULT_LIMB_DARKENING,
) -> pl.Figure:
    """Phase-folded transit of each sector with its best-fit model and ``b``.

    Folded on the ephemeris shifted by the sector's fitted mid-time offset, so
    a changing ``b`` shows up as a changing transit shape rather than a shift.
    The last panel overlays every sector's model (coloured by time) to compare them.
    """
    sectors = {
        s.span.sector: s
        for s in extract_sector_transits(time, flux, flux_err, spans, period, t0, duration)
    }
    fits = [f for f in result.sector_fits if f.sector in sectors]
    if not fits:
        fig, ax = pl.subplots(figsize=(5, 3))
        ax.text(0.5, 0.5, "no sector with data", ha="center", transform=ax.transAxes)
        return fig
    n_panels = len(fits) + 1
    ncols = min(FOLD_MAX_COLS, n_panels)
    nrows = int(np.ceil(n_panels / ncols))
    fig, axes = pl.subplots(
        nrows, ncols, figsize=(4.2 * ncols, 3.2 * nrows), sharex=True, sharey=True, squeeze=False
    )
    flat_axes = axes.ravel()
    half = min(max(HALF_WINDOW_DURATIONS * duration, MIN_HALF_WINDOW), 0.25 * period)
    colors = pl.cm.viridis(np.linspace(0, 0.9, len(fits)))
    overlay = flat_axes[len(fits)]
    for ax, fit, color in zip(flat_axes, fits, colors):
        data = sectors[fit.sector]
        t0_eff = t0 + fit.dt0
        ph = ((data.time - t0_eff + 0.5 * period) % period - 0.5 * period) * 24.0
        norm = data.flux / fit.baseline
        ax.plot(ph, norm, ".", ms=2, color="0.75", rasterized=True)
        n_bins = max(int(np.sqrt(ph.size)), FOLD_BIN_POINTS)
        bx, by = _binned(ph, norm, n_bins)
        ax.plot(bx, by, "o", ms=3, color="k")
        phase, curve = _fold_model(fit, result.a_rs, period, t0_eff, half, limb_darkening)
        ax.plot(phase * 24.0, curve, "-", color="C3", lw=1.5)
        overlay.plot(phase * 24.0, curve, "-", color=color, lw=1.5, label=f"S{fit.sector}")
        ax.set_title(
            f"S{fit.sector}: b = {fit.b:.2f} ± {fit.b_err:.2f}, T14 = {fit.t14 * 24:.2f} h",
            fontsize=9,
        )
    overlay.set_title("best-fit models", fontsize=9)
    overlay.legend(fontsize=7, loc="lower right", ncol=2)
    for ax in flat_axes[n_panels:]:
        ax.set_visible(False)
    for ax in axes[-1]:
        ax.set_xlabel("hours from mid-transit")
    for ax in axes[:, 0]:
        ax.set_ylabel("normalised flux")
    verdict = "candidate" if result.is_candidate else "no candidate"
    fig.suptitle(
        f"{title}  phase-folded per sector, P = {period:.5f} d → {verdict}".strip(), fontsize=10
    )
    fig.tight_layout()
    return fig
