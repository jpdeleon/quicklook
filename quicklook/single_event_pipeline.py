"""Single-event search glue for ``TessQuickLook`` (kept out of tql.py for size).

``SingleEventMixin`` expects the host to provide ``raw_lc``, ``flat_lc``,
``tls_results``, ``exofop_data``, ``window_length``, ``verbose`` and
``single_event_search``. The search detrends ``raw_lc`` itself (see
``quicklook.single_event_detrend``); ``flat_lc`` stays TLS's.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

import astropy.units as u
import matplotlib.pyplot as pl
import numpy as np
from loguru import logger

from quicklook.broken_axis import broken_time_axes, finish_broken_axes
from quicklook.duo import DuoPair, find_duo
from quicklook.plot_single import (
    mark_events,
    plot_alias_weights,
    plot_duo_zoom,
    plot_event_zoom,
    plot_period_constraint,
    plot_single_event_statistic,
)
from quicklook.single_event import (
    SingleEvent,
    SingleEventResult,
    _cadence,
    classify_events,
    duration_grid,
    longest_segment,
    min_period,
    single_event_search,
)
from quicklook.single_event_detrend import (
    WINDOW_TO_DURATION,
    choose_window,
    detrend,
    event_mask,
    is_fast_rotator,
)
from quicklook.utils import get_params_from_exofop

# TLS signals weaker than this do not make a target "periodic" (SDE 7 is
# roughly a 1% false-alarm level in Hippke & Heller 2019). Two transits are not
# enough either: TLS's period is then just one of the dT/n aliases, which the
# duotransit pairing lists instead.
PERIODIC_SDE_MIN = 7.0
PERIODIC_MIN_TRANSITS = 3
MONO_FIGURE_SUFFIX = "_mono.png"
MAX_REDETREND_ITER = 3  # mask-and-re-detrend passes after the first search

_EVENT_ARRAY_FIELDS = {
    "event_t0s": "t0",
    "event_durations": "duration",
    "event_depths": "depth",
    "event_snrs": "snr",
}
_PRIMARY_FIELDS = (
    "t0",
    "duration",
    "depth",
    "snr",
    "delta_bic_null",
    "delta_bic_systematics",
    "dist_to_gap",
)


def _density(params) -> float | None:
    """M/R^3 in solar units from an ExoFOP stellar-parameter dict, or None."""
    try:
        radius, mass = float(params["srad"]), float(params["mass"])
    except (TypeError, KeyError, ValueError):
        return None
    if not (np.isfinite(radius) and np.isfinite(mass) and radius > 0 and mass > 0):
        return None
    return mass / radius**3


class SingleEventMixin:
    """Runs the single-event search and draws its panels."""

    # Provided by the host class (TessQuickLook).
    raw_lc: Any
    flat_lc: Any
    window_length: float
    tls_results: Any
    exofop_data: Any
    verbose: bool

    single_event_result: Optional[SingleEventResult] = None
    unexplained_events: tuple[SingleEvent, ...] = ()
    primary_event: Optional[SingleEvent] = None
    event_class: Optional[str] = None
    p_min: float = np.nan
    duo: Optional[DuoPair] = None
    _panel0_axes: Optional[list] = None
    # The search's own detrended light curve and settings (PR 3 of #42).
    se_time: Optional[np.ndarray] = None
    se_flux: Optional[np.ndarray] = None
    se_trend: Optional[np.ndarray] = None
    se_window: float = np.nan
    se_window_reason: str = ""
    se_max_duration: float = np.nan
    se_redetrended: bool = False
    _broken_panels: dict = {}

    # --- search -------------------------------------------------------------

    def _stellar_density(self) -> float:
        """Stellar density in solar units from ExoFOP, or 1 if unknown.

        Uses the TIC row (index 1, as the summary panel does) when it has both
        radius and mass, otherwise the newest entry that does. The newest entry
        alone can lack a mass (e.g. TOI-2180), which would silently give solar
        density.
        """
        entries = list((self.exofop_data or {}).get("stellar_parameters") or [])
        tic_row = entries[1:2]
        latest = get_params_from_exofop(self.exofop_data or {}, name="stellar_parameters")
        for params in tic_row + ([latest] if latest else []) + entries[::-1]:
            density = _density(params)
            if density is not None:
                return density
        return 1.0

    def _tls_is_periodic(self) -> bool:
        res = self.tls_results
        if res is None:
            return False
        sde = res.get("SDE")
        n_transits = res.get("distinct_transit_count", PERIODIC_MIN_TRANSITS)
        return bool(
            sde is not None
            and np.isfinite(sde)
            and sde >= PERIODIC_SDE_MIN
            and n_transits >= PERIODIC_MIN_TRANSITS
        )

    def _explained_by_tls(self, event) -> bool:
        res = self.tls_results
        period, epoch, dur = res.get("period"), res.get("T0"), res.get("duration")
        if not all(v is not None and np.isfinite(v) for v in (period, epoch, dur)):
            return False
        phase = ((event.t0 - epoch + period / 2) % period) - period / 2
        return abs(phase) < max(dur, event.duration)

    def run_single_event_search(self):
        """Search the flattened light curve for individual dips and classify the target.

        Adds ``event_class``, ``n_events``, ``n_unexplained_events``, the
        strongest unexplained event (``single_event_*``), ``p_min`` and the
        per-event arrays to ``tls_results`` (and hence the HDF5). Does nothing
        unless ``single_event_search`` is enabled.
        """
        if not getattr(self, "single_event_search", False):
            return None
        rho_star = self._stellar_density()
        time, flux, res = self._detrend_and_search(rho_star)
        periodic = self._tls_is_periodic()
        unexplained = tuple(e for e in res.events if not (periodic and self._explained_by_tls(e)))
        self.single_event_result = res
        self.unexplained_events = unexplained
        self.primary_event = unexplained[0] if unexplained else None
        self.duo = None if periodic else find_duo(time, unexplained, rho_star)
        self.event_class = classify_events(
            len(res.events), periodic, len(unexplained), has_duo=self.duo is not None
        )
        ev = self.primary_event
        if self.duo is not None:
            self.p_min = self.duo.p_min
        else:
            self.p_min = min_period(time, ev.t0, ev.duration) if ev is not None else np.nan
        self._store_single_event_results()
        if self.verbose:
            self._log_single_event_results()
        return res

    def _raw_normalised(self):
        t = np.asarray(self.raw_lc.time.value, dtype=float)
        f = np.asarray(self.raw_lc.flux.value, dtype=float)
        return t, f / np.nanmedian(f)

    def _search_durations(self, grid) -> np.ndarray:
        """The stellar duration grid, capped at window / 3 (with a warning)."""
        cap = self.se_window / WINDOW_TO_DURATION
        self.se_max_duration = float(min(grid[-1], cap))
        if grid[-1] <= cap:
            return grid
        logger.warning(
            f"Single-event search: detrending window {self.se_window:.2f} d preserves "
            f"dips up to ~{cap * 24:.1f} h, so the first pass skips longer trial "
            f"durations (up to {grid[-1] * 24:.1f} h); detected events are then "
            "re-measured on a light curve detrended with them masked."
        )
        return np.geomspace(grid[0], max(cap, 2 * grid[0]), grid.size)

    def _detrend_and_search(self, rho_star):
        """Detrend ``raw_lc`` for the search, search, then repeat with events masked.

        The first pass caps trial durations at window / 3, since longer dips are
        partly removed by the detrending. With the events masked the trend no
        longer depends on the window inside them, so the second pass uses the
        full stellar grid and measures long transits at their real duration.
        """
        time, raw = self._raw_normalised()
        # LOCoR also disables the rotation-aware window but leaves it at 0 d;
        # only a positive window counts as one the user chose.
        user_set = not getattr(self, "_use_rotation_aware_window", True) and (
            float(self.window_length) > 0
        )
        fast = is_fast_rotator(time, raw)
        self.se_window, self.se_window_reason = choose_window(
            float(self.window_length), rho_star, user_set, fast
        )
        full_grid = duration_grid(longest_segment(time), _cadence(time), rho_star)
        durations = self._search_durations(full_grid)
        edge = float(getattr(self, "edge_cutoff", 0.0) or 0.0)
        flux, trend = detrend(time, raw, self.se_window, edge_cutoff=edge)
        res = single_event_search(time, flux, rho_star=rho_star, durations=durations)
        self.se_redetrended = bool(res.events)
        mask = np.zeros(time.size, dtype=bool)
        for _ in range(MAX_REDETREND_ITER):
            grown = mask | event_mask(time, res.events)
            if not res.events or np.array_equal(grown, mask):
                break
            # The mask only grows: a first-pass event capped below its real
            # duration masks too little, and the next pass measures it longer.
            mask = grown
            flux, trend = detrend(time, raw, self.se_window, mask=mask, edge_cutoff=edge)
            res = single_event_search(time, flux, rho_star=rho_star, durations=full_grid)
        self.se_time, self.se_flux, self.se_trend = time, flux, trend
        return time, flux, res

    def _store_single_event_results(self):
        res, ev = self.single_event_result, self.primary_event
        assert res is not None
        values = {
            "single_event_search": True,
            "event_class": self.event_class,
            "n_events": len(res.events),
            "n_unexplained_events": len(self.unexplained_events),
            "p_min": float(self.p_min),
            "single_event_window": float(self.se_window),
            "single_event_window_reason": self.se_window_reason,
            "single_event_max_duration": float(self.se_max_duration),
            "single_event_redetrended": bool(self.se_redetrended),
        }
        for name in _PRIMARY_FIELDS:
            values[f"single_event_{name}"] = float(getattr(ev, name)) if ev else np.nan
        for key, attr in _EVENT_ARRAY_FIELDS.items():
            values[key] = np.array([getattr(e, attr) for e in res.events], dtype=float)
        values.update(self._duo_values())
        for key, val in values.items():
            self.tls_results[key] = val

    def _duo_values(self) -> dict:
        duo = self.duo
        if duo is None:
            return {
                "duo_t0s": np.array([], dtype=float),
                "duo_delta_t": np.nan,
                "duo_depth_ratio": np.nan,
                "n_aliases": 0,
                "p_aliases": np.array([], dtype=float),
                "p_alias_weights": np.array([], dtype=float),
                "p_best": np.nan,
            }
        return {
            "duo_t0s": np.array([duo.first.t0, duo.second.t0]),
            "duo_delta_t": float(duo.delta_t),
            "duo_depth_ratio": float(duo.depth_ratio),
            "n_aliases": len(duo.aliases),
            "p_aliases": np.array([a.period for a in duo.aliases]),
            "p_alias_weights": np.array([a.weight for a in duo.aliases]),
            "p_best": float(duo.best.period),
        }

    def _log_single_event_results(self):
        ev = self.primary_event
        msg = (
            f"Single-event search: {len(self.single_event_result.events)} event(s), "
            f"{len(self.unexplained_events)} not explained by TLS -> class={self.event_class} "
            f"(window {self.se_window:.2f} d [{self.se_window_reason}], durations <= "
            f"{self.se_max_duration * 24:.1f} h"
            f"{', re-detrended with events masked' if self.se_redetrended else ''})"
        )
        if self.duo is not None:
            duo = self.duo
            msg += (
                f"; duotransit T0s={duo.first.t0:.3f},{duo.second.t0:.3f} BTJD "
                f"(dT={duo.delta_t:.2f} d), {len(duo.aliases)} allowed aliases, "
                f"most likely P={duo.best.period:.2f} d, P_min={self.p_min:.2f} d"
            )
        elif ev is not None:
            msg += (
                f"; strongest at T0={ev.t0:.3f} BTJD, depth={ev.depth * 1e3:.2f} ppt, "
                f"duration={ev.duration * 24:.1f} h, SNR={ev.snr:.1f}, P_min={self.p_min:.1f} d"
            )
        logger.info(msg)

    # --- summary panel --------------------------------------------------------

    def _single_event_candidate_rows(self, r_star: float) -> list:
        """Candidate-column rows for the mono layout of the summary panel."""
        n_events = len(self.single_event_result.events) if self.single_event_result else 0
        pm_rows = [("Class", f"{self.event_class}"), ("Num. events", f"{n_events}")]
        if self.duo is not None:
            return pm_rows + self._duo_candidate_rows(r_star)
        ev = self.primary_event
        if ev is None:
            return pm_rows
        rows = pm_rows + [
            ("T0 (event)", f"{ev.t0:.4f} BTJD"),
            ("Dur. (event)", f"{ev.duration * 24:.2f} hr"),
            ("Depth (event)", f"{ev.depth * 1e3:.2f} ppt"),
            ("SNR (event)", f"{ev.snr:.1f}"),
            (r"$\Delta$BIC (vs flat)", f"{ev.delta_bic_null:.0f}"),
            (r"$\Delta$BIC (vs sys.)", f"{ev.delta_bic_systematics:.0f}"),
            ("P min", f"{self.p_min:.1f} d"),
            ("Gap distance", f"{ev.dist_to_gap:.2f} d"),
        ]
        if np.isfinite(r_star) and ev.depth > 0:
            rp = np.sqrt(ev.depth) * r_star * u.Rsun.to(u.Rearth)
            rows.append(("Rp (event)", f"{rp:.2f} " + r"R$_{\oplus}$"))
        return rows

    def _duo_candidate_rows(self, r_star: float) -> list:
        duo = self.duo
        assert duo is not None
        depth = (duo.first.depth + duo.second.depth) / 2
        rows = [
            ("T0 (1st)", f"{duo.first.t0:.4f} BTJD"),
            ("T0 (2nd)", f"{duo.second.t0:.4f} BTJD"),
            (r"$\Delta$T", f"{duo.delta_t:.3f} d"),
            ("Dur. (mean)", f"{(duo.first.duration + duo.second.duration) * 12:.2f} hr"),
            ("Depth (mean)", f"{depth * 1e3:.2f} ppt"),
            ("Depth ratio", f"{duo.depth_ratio:.2f}"),
            ("SNR (1st, 2nd)", f"{duo.first.snr:.1f}, {duo.second.snr:.1f}"),
            ("Aliases", f"{len(duo.aliases)}"),
            ("P (most likely)", f"{duo.best.period:.2f} d ({duo.best.weight:.0%})"),
            ("P min", f"{self.p_min:.2f} d"),
        ]
        if np.isfinite(r_star) and depth > 0:
            rp = np.sqrt(depth) * r_star * u.Rsun.to(u.Rearth)
            rows.append(("Rp (mean)", f"{rp:.2f} " + r"R$_{\oplus}$"))
        return rows

    # --- figure ---------------------------------------------------------------

    def _time_panel_axes(self, ax) -> list:
        """``[ax]``, or one axes per sector for mono/duo layouts with sector gaps."""
        if self.event_class not in ("mono", "duo"):
            return [ax]
        return broken_time_axes(ax, self.flat_lc.time.value)

    def _se_light_curve(self):
        """The search's detrended light curve (falls back to ``flat_lc``)."""
        if self.se_time is not None and self.se_flux is not None:
            return self.se_time, self.se_flux
        return self.flat_lc.time.value, self.flat_lc.flux.value

    def _plot_se_trend(self, ax):
        """Overlay the single-event search's own trend on the raw light-curve panel."""
        if self.se_time is None or self.se_trend is None:
            return
        scale = float(np.nanmedian(self.raw_lc.flux.value))
        masked = " (events masked)" if self.se_redetrended else ""
        ax.plot(
            self.se_time,
            self.se_trend * scale,
            color="C2",
            ls="--",
            lw=1.5,
            label=f"single-event trend\nwindow={self.se_window:.2f} d{masked}",
        )
        ax.legend(fontsize="x-small", loc="best")

    def _plot_single_event_panels(self, axes):
        """Fill panels 1, 4, 5, 7 and 8 (indices 0, 3, 4, 6, 7) for the mono/duo layout."""
        assert self.single_event_result is not None
        time, flux = self._se_light_curve()
        events = self.unexplained_events
        panel0 = self._panel0_axes or [axes[0]]
        for ax in panel0:
            mark_events(ax, events)
            self._plot_se_trend(ax)
        panel3 = self._time_panel_axes(axes[3])
        for ax in panel3:
            ax.plot(time, flux, ".", color="k", ms=1)
            mark_events(ax, events)
            ax.set_xlabel("Time [BTJD]")
            ax.set_ylabel("Normalized flux")
        finish_broken_axes(panel3)
        panel4 = self._time_panel_axes(axes[4])
        for ax in panel4:
            plot_single_event_statistic(self.single_event_result, ax, events)
        finish_broken_axes(panel4)
        self._broken_panels = {0: panel0, 3: panel3, 4: panel4}
        ev = self.primary_event
        if ev is None:
            for ax in (axes[6], axes[7]):
                ax.text(0.5, 0.5, "No single event", ha="center", transform=ax.transAxes)
                ax.axis("off")
            return
        if self.duo is not None:
            plot_duo_zoom(time, flux, self.duo, axes[6])
            plot_alias_weights(self.duo, axes[7])
            return
        plot_event_zoom(time, flux, ev, axes[6])
        plot_period_constraint(time, ev, axes[7])

    def _save_single_event_figure(self, fp) -> Path | None:
        """Write ``<fp>_mono.png`` when a periodic target also has unexplained events."""
        ev = self.primary_event
        res = self.single_event_result
        if self.event_class != "periodic" or ev is None or res is None:
            return None
        time, flux = self._se_light_curve()
        fig, axes = pl.subplots(2, 2, figsize=(12, 8), tight_layout=True)
        axes = axes.flatten()
        plot_single_event_statistic(res, axes[0], self.unexplained_events)
        plot_event_zoom(time, flux, ev, axes[1])
        plot_period_constraint(time, ev, axes[2])
        axes[3].axis("off")
        lines = [f"{label}: {value}" for label, value in self._single_event_candidate_rows(np.nan)]
        axes[3].text(0, 1, "\n".join(lines), va="top", family="monospace", fontsize="small")
        fig.suptitle("Single event not explained by the periodic TLS signal")
        out = Path(str(fp) + MONO_FIGURE_SUFFIX)
        fig.savefig(out, dpi=120, bbox_inches="tight")
        pl.close(fig)
        logger.info(f"Saved: {out}")
        return out
