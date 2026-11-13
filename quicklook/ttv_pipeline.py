"""TTV search glue for ``TessQuickLook`` (step 5 of #42; kept out of tql.py for size).

``TTVMixin`` expects the host to provide ``flat_lc``, ``tls_results``,
``verbose`` and ``ttv_search``, plus ``SingleEventMixin``'s state for chaining
single events. Two entry points:

* ``run_ttv_search`` right after TLS, from the TLS ephemeris, so the iterative
  re-search can mask transits at their fitted times (``primary_transit_mask``).
* ``run_ttv_from_chain`` after the single-event search, when TLS found nothing
  periodic but three or more single events line up within a tolerance.
"""

from __future__ import annotations

from typing import Any, Optional

import numpy as np
from loguru import logger

from quicklook.plot_ttv import plot_aligned_transits, plot_oc
from quicklook.ttv import (
    TTV_MIN_TRANSITS,
    CycleSolution,
    TransitTime,
    TTVStats,
    chain_events,
    cycle_count_solutions,
    default_delta,
    fit_transit_times,
    ttv_mask,
    ttv_statistics,
)

MINUTES_PER_DAY = 1440.0
_TTV_ARRAY_FIELDS = {
    "ttv_epochs": "epoch",
    "ttv_tcs": "tc",
    "ttv_tc_errs": "tc_err",
    "ttv_durations": "duration",
    "ttv_depths": "depth",
    "ttv_snrs": "snr",
}


class TTVMixin:
    """Per-transit timing fits, O-C statistics and the TTV figure panels."""

    # Provided by the host class (TessQuickLook / SingleEventMixin).
    flat_lc: Any
    tls_results: Any
    verbose: bool

    ttv_times: tuple[TransitTime, ...] = ()
    ttv_stats: Optional[TTVStats] = None
    ttv_source: str = ""
    ttv_delta: float = np.nan
    ttv_cycle_solutions: list = []
    _ttv_lc: Optional[tuple] = None  # (time, flux) the times were fitted on

    # --- search -------------------------------------------------------------

    def run_ttv_search(self):
        """Fit every transit of a periodic TLS signal; no-op otherwise.

        Needs ``--ttv-search`` and a TLS detection that passes the same
        periodic test as the single-event search (SDE >= 7, >= 3 transits).
        """
        if not getattr(self, "ttv_search", False):
            return None
        if not self._tls_is_periodic():
            if self.verbose:
                logger.info("TTV search: TLS signal is not periodic; nothing to time.")
            self._store_ttv_results()
            return None
        res = self.tls_results
        time, flux = self.flat_lc.time.value, self.flat_lc.flux.value
        return self._fit_and_store(time, flux, res["period"], res["T0"], res["duration"], "tls")

    def run_ttv_from_chain(self):
        """Chain three or more single events into a TTV-tolerant ephemeris.

        Runs only when the TLS-based search found nothing to time. A chain
        that yields at least three fitted transits makes the target
        ``periodic`` (the period is the chain's, not TLS's).
        """
        if not getattr(self, "ttv_search", False) or self.ttv_stats is not None:
            return None
        if not getattr(self, "single_event_search", False) or self.single_event_result is None:
            return None
        chain = chain_events(self.unexplained_events)
        if chain is None:
            return None
        period, t0, members = chain
        duration = float(np.median([e.duration for e in members]))
        time, flux = self._se_light_curve()
        stats = self._fit_and_store(time, flux, period, t0, duration, "chain")
        if stats is not None and stats.n >= TTV_MIN_TRANSITS:
            self.event_class = "periodic"
            self.tls_results["event_class"] = "periodic"
        return stats

    def _fit_and_store(self, time, flux, period, t0, duration, source):
        period, t0, duration = float(period), float(t0), float(duration)
        self.ttv_delta = default_delta(period, duration)
        self.ttv_times = fit_transit_times(time, flux, period, t0, duration, self.ttv_delta)
        self.ttv_stats = ttv_statistics(self.ttv_times)
        self.ttv_source = source
        self._ttv_lc = (np.asarray(time), np.asarray(flux))
        self.ttv_cycle_solutions = (
            cycle_count_solutions(self.ttv_times) if self.ttv_stats.n >= 2 else []
        )
        self._store_ttv_results()
        if self.verbose:
            self._log_ttv_results()
        return self.ttv_stats

    # --- outputs ------------------------------------------------------------

    def _store_ttv_results(self):
        stats, eph = self.ttv_stats, getattr(self.ttv_stats, "ephemeris", None)
        alternatives = [s for s in self.ttv_cycle_solutions if any(s.shifts)]
        values = {
            "ttv_search": True,
            "ttv_source": self.ttv_source,
            "ttv_detected": bool(stats.detected) if stats else False,
            "ttv_amplitude": float(stats.amplitude) if stats else np.nan,
            "ttv_rms": float(stats.rms) if stats else np.nan,
            "ttv_pvalue": float(stats.pvalue) if stats else np.nan,
            "ttv_pvalue_robust": float(stats.pvalue_robust) if stats else np.nan,
            "n_ttv_transits": len(self.ttv_times),
            "ttv_period": eph.period if eph else np.nan,
            "ttv_t0": eph.t0 if eph else np.nan,
            "ttv_delta": float(self.ttv_delta),
            "ttv_oc": np.asarray(stats.oc if stats else [], dtype=float),
            "ttv_n_cycle_solutions": len(self.ttv_cycle_solutions),
            "ttv_alt_periods": np.array([s.period for s in alternatives], dtype=float),
        }
        for key, attr in _TTV_ARRAY_FIELDS.items():
            values[key] = np.array([getattr(tt, attr) for tt in self.ttv_times], dtype=float)
        for key, val in values.items():
            self.tls_results[key] = val

    def _log_ttv_results(self):
        stats = self.ttv_stats
        msg = f"TTV search ({self.ttv_source}): {len(self.ttv_times)} transit(s) fitted"
        if stats is not None and stats.ephemeris is not None:
            msg += (
                f", P={stats.ephemeris.period:.5f} d, O-C amplitude "
                f"{stats.amplitude * MINUTES_PER_DAY:.1f} min, p={stats.pvalue:.2g} "
                f"(robust {stats.pvalue_robust:.2g}) -> "
                f"{'TTVs detected' if stats.detected else 'consistent with a linear ephemeris'}"
            )
        if len(self.ttv_cycle_solutions) > 1:
            periods = ", ".join(f"{s.period:.5f}" for s in self.ttv_cycle_solutions)
            msg += f"; cycle count ambiguous across gaps, P in {{{periods}}} d"
        logger.info(msg)

    # --- masking for the iterative search ------------------------------------

    def primary_transit_mask(self, lc, margin: float = 1.2) -> np.ndarray:
        """Mask the primary signal at its fitted times where available.

        Epochs without a fitted time keep the linear-ephemeris mask, so a
        transit too weak to fit on its own is still removed.
        """
        linear = self._transit_mask_for_lc(lc, self.tls_results, margin=margin)
        if not self.ttv_times or self.ttv_source != "tls":
            return linear
        time = np.asarray(lc.time.value, dtype=float)
        res = self.tls_results
        epochs = np.round((time - res["T0"]) / res["period"]).astype(int)
        fitted = {tt.epoch for tt in self.ttv_times}
        unfitted = ~np.isin(epochs, list(fitted))
        return ttv_mask(time, self.ttv_times, margin=margin) | (linear & unfitted)

    # --- figure ---------------------------------------------------------------

    @property
    def has_ttv_layout(self) -> bool:
        """TTV panels replace odd/even and secondary when >= 3 transits were fitted."""
        return (
            self.ttv_stats is not None
            and self.ttv_stats.n >= TTV_MIN_TRANSITS
            and getattr(self, "event_class", None) not in ("mono", "duo")
        )

    def ttv_transit_mask(self, lc, margin: float = 1.0) -> np.ndarray:
        """Cadences of ``lc`` inside the fitted transits (for panel 4)."""
        return ttv_mask(np.asarray(lc.time.value, dtype=float), self.ttv_times, margin=margin)

    def _plot_ttv_panels(self, axes):
        """Panel 7: transits aligned on their fitted times. Panel 8: O-C."""
        time, flux = self._ttv_lc
        plot_aligned_transits(time, flux, self.ttv_times, ax=axes[6])
        alternatives: list[CycleSolution] = [s for s in self.ttv_cycle_solutions if any(s.shifts)]
        plot_oc(self.ttv_stats, ax=axes[7], n_alternatives=len(alternatives))

    def _ttv_candidate_rows(self) -> list:
        stats = self.ttv_stats
        eph = stats.ephemeris
        rows = [
            ("TTV source", "TLS" if self.ttv_source == "tls" else "chained single events"),
            ("Mean period", f"{eph.period:.5f}" + r"$\pm$" + f"{eph.period_err:.5f} d"),
            ("Transits timed", f"{stats.n}"),
            ("TTV amplitude", f"{stats.amplitude * MINUTES_PER_DAY:.1f} min"),
            ("TTV p-value", f"{stats.pvalue:.2g} (robust {stats.pvalue_robust:.2g})"),
            ("TTVs", "detected" if stats.detected else "not significant"),
        ]
        if len(self.ttv_cycle_solutions) > 1:
            rows.append(("Cycle solutions", f"{len(self.ttv_cycle_solutions)}"))
        return rows


__all__ = ["TTVMixin"]
