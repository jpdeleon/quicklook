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

from pathlib import Path
from typing import Any, Optional

import matplotlib.pyplot as pl
import numpy as np
from loguru import logger

from quicklook.plot_ttv import plot_aligned_transits, plot_oc
from quicklook.plot_ttv_local import plot_transit_gallery, plot_ttv_diagnostics
from quicklook.ttv import (
    TTV_MIN_TRANSITS,
    CycleSolution,
    TransitTime,
    TTVStats,
    chain_events,
    cycle_count_solutions,
    default_delta,
    fit_transit_times,
    inflate_errors,
    red_noise_beta,
    ttv_mask,
    ttv_statistics,
)
from quicklook.ttv_local import LocalTTVResult, run_local_ttv

MINUTES_PER_DAY = 1440.0
TTV_FIGURE_SUFFIX = "_ttv.png"
TTV_GALLERY_SUFFIX = "_ttv_transits.png"
MIN_EPHEM_PERIOD = 0.3  # days; TFOP fills missing values with 0.1
TTV_PERIOD_MATCH_TOL = 0.01  # timed period vs catalogue period (relative)
ROTATION_ALIAS_TOL = 0.03  # same tolerance as tql's rotation_alias_risk
ROTATION_ALIAS_HARMONICS = (1, 2, 3, 1 / 2, 1 / 3)
_TTV_SOURCE_LABELS = {
    "tls": "TLS",
    "chain": "chained single events",
    "custom": "custom ephemeris",
}
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
    ttv_beta: float = 1.0
    ttv_period_match: Optional[bool] = None  # None: no catalogue period to compare
    ttv_rotation_alias: Optional[bool] = None  # None: no rotation period
    ttv_cycle_solutions: list = []
    _ttv_lc: Optional[tuple] = None  # (time, flux) the times were fitted on
    ttv_local: Optional[LocalTTVResult] = None
    ttv_injections: int = 10

    # --- search -------------------------------------------------------------

    def run_ttv_search(self):
        """Fit every transit of the user's ephemeris, else of a periodic TLS signal.

        Needs ``--ttv-search``. A ``--custom-ephem`` (not masked out with
        ``--mask-ephem``) is timed whatever TLS found: on young, spotted stars TLS
        locks onto the rotation period and would time starspots instead. Without
        one, the TLS detection must pass the same periodic test as the
        single-event search (SDE >= 7, >= 3 transits).
        """
        if not getattr(self, "ttv_search", False):
            return None
        custom = self._custom_ephemeris()
        if custom is not None:
            period, t0, duration = custom
            if self.verbose:
                logger.info(
                    f"TTV search: timing the custom ephemeris (P={period:.5f} d) "
                    f"instead of the TLS signal (P={float(self.tls_results['period']):.5f} d)."
                )
            time, flux = self.flat_lc.time.value, self.flat_lc.flux.value
            return self._fit_and_store(time, flux, period, t0, duration, "custom")
        if not self._tls_is_periodic():
            if self.verbose:
                logger.info("TTV search: TLS signal is not periodic; nothing to time.")
            self._store_ttv_results()
            return None
        res = self.tls_results
        time, flux = self.flat_lc.time.value, self.flat_lc.flux.value
        return self._fit_and_store(time, flux, res["period"], res["T0"], res["duration"], "tls")

    def _custom_ephemeris(self):
        """(period, t0, duration) in days from ``--custom-ephem``, else None.

        None when the ephemeris came from the catalogue, when ``--mask-ephem``
        removed those transits, or when the period or duration is unusable.
        """
        if getattr(self, "ephem_source", None) != "custom" or getattr(self, "mask_ephem", False):
            return None
        period = getattr(self, "toi_period", None)
        epoch = getattr(self, "toi_epoch", None)
        dur = getattr(self, "toi_dur", None)
        if period is None or epoch is None or dur is None:
            return None
        period, t0, duration = float(period[0]), float(epoch[0]), float(dur[0])
        if not (np.isfinite(period) and period >= MIN_EPHEM_PERIOD and np.isfinite(t0)):
            return None
        if not (np.isfinite(duration) and duration > 0):
            return None
        return period, t0, duration

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
        times = fit_transit_times(time, flux, period, t0, duration, self.ttv_delta)
        # White-noise errors understate the O-C scatter of active stars.
        self.ttv_beta = red_noise_beta(time, flux, period, t0, duration)
        self.ttv_times = inflate_errors(times, self.ttv_beta)
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

    # --- vetting ---------------------------------------------------------------

    def _catalogue_periods(self) -> list:
        """Known periods (d): the loaded ephemeris plus every ExoFOP planet entry."""
        periods = []
        toi = getattr(self, "toi_period", None)
        if toi is not None:
            periods.append(toi[0])
        exofop = getattr(self, "exofop_data", None) or {}
        for entry in exofop.get("planet_parameters") or []:
            periods.append(entry.get("per") if isinstance(entry, dict) else None)
        out = []
        for p in periods:
            try:
                p = float(p)
            except (TypeError, ValueError):
                continue
            if np.isfinite(p) and p >= MIN_EPHEM_PERIOD:
                out.append(p)
        return out

    def _rotation_period(self):
        res = self.tls_results
        gls = res.get("Prot_gls") if hasattr(res, "get") else None
        prot = gls[0] if gls is not None else getattr(self, "Prot_ls", None)
        try:
            prot = float(prot)
        except (TypeError, ValueError):
            return None
        return prot if np.isfinite(prot) and prot > 0 else None

    def vet_ttv_detection(self) -> bool:
        """Final TTV verdict: raw detection, on a catalogue period, not rotation.

        Runs after the rotation period is known (from ``append_tls_results``).
        ``ttv_detected_raw`` keeps the statistical verdict; ``ttv_detected`` is
        False when the timed period matches no catalogue period (when there is
        one) or aliases the stellar rotation period.
        """
        stats = self.ttv_stats
        eph = getattr(stats, "ephemeris", None)
        raw = bool(stats.detected) if stats is not None else False
        period = eph.period if eph is not None else np.nan
        self.ttv_period_match = self.ttv_rotation_alias = None
        if np.isfinite(period):
            known = self._catalogue_periods()
            if known:
                self.ttv_period_match = any(
                    abs(period - p) <= TTV_PERIOD_MATCH_TOL * p for p in known
                )
            prot = self._rotation_period()
            if prot is not None:
                self.ttv_rotation_alias = any(
                    abs(period - n * prot) <= ROTATION_ALIAS_TOL * n * prot
                    for n in ROTATION_ALIAS_HARMONICS
                )
        vetted = raw and self.ttv_period_match is not False and not self.ttv_rotation_alias
        for key, val in (
            ("ttv_detected_raw", raw),
            ("ttv_detected", vetted),
            ("ttv_period_match", self.ttv_period_match),
            ("ttv_rotation_alias", self.ttv_rotation_alias),
        ):
            self.tls_results[key] = val
        if raw and not vetted:
            reasons = []
            if self.ttv_period_match is False:
                reasons.append("timed period matches no catalogue period")
            if self.ttv_rotation_alias:
                reasons.append("timed period aliases the rotation period")
            logger.info(f"TTV vetting: detection rejected ({'; '.join(reasons)}).")
        return vetted

    @property
    def ttv_detected(self) -> bool:
        """Vetted verdict when vetting ran, else the statistical one."""
        res = self.tls_results
        if hasattr(res, "get") and res.get("ttv_detected_raw") is not None:
            return bool(res.get("ttv_detected"))
        return bool(self.ttv_stats is not None and self.ttv_stats.detected)

    # --- local-baseline diagnostics (byproduct figures) ---------------------

    def _diagnostic_ephemeris(self):
        """(period, t0, duration, source): timed TLS/chain ephemeris, else TOI/custom."""
        stats = self.ttv_stats
        if (
            stats is not None
            and stats.ephemeris is not None
            and stats.n >= TTV_MIN_TRANSITS
            and np.isfinite(stats.ephemeris.period)
            and np.isfinite(stats.ephemeris.t0)
        ):
            dur = float(np.median([tt.duration for tt in self.ttv_times]))
            return stats.ephemeris.period, stats.ephemeris.t0, dur, self.ttv_source
        period = getattr(self, "toi_period", None)
        epoch = getattr(self, "toi_epoch", None)
        dur = getattr(self, "toi_dur", None)
        if period is None or epoch is None or float(period[0]) < MIN_EPHEM_PERIOD:
            return None
        if getattr(self, "mask_ephem", False):  # those transits were cut from raw_lc
            return None
        duration = float(dur[0]) if dur is not None and float(dur[0]) * 24 > 0.5 else None
        source = "custom" if getattr(self, "ephem_source", None) == "custom" else "toi"
        return float(period[0]), float(epoch[0]), duration, source

    def run_ttv_diagnostics(self):
        """Local-baseline timing of every transit on the raw flux (``*_ttv.png``).

        Runs with ``--ttv-search``. Uses the timed TLS/chain ephemeris when the TTV
        search fitted >= 3 transits, otherwise the TOI (or ``--ephem``) ephemeris,
        so it still works when TLS locks onto stellar variability.
        """
        if not getattr(self, "ttv_search", False):
            return None
        ephem = self._diagnostic_ephemeris()
        if ephem is None:
            if self.verbose:
                logger.info("TTV diagnostics: no timed or TOI ephemeris; skipped.")
            return None
        period, t0, duration, source = ephem
        lc = self.raw_lc
        time = np.asarray(lc.time.value, dtype=float)
        flux = np.asarray(lc.flux.value, dtype=float)
        flux = flux / np.nanmedian(flux)
        self._ttv_raw = (time, flux)
        self.ttv_local = run_local_ttv(
            time, flux, period, t0, duration, source=source, n_injections=int(self.ttv_injections)
        )
        self._store_ttv_local()
        if self.verbose and self.ttv_local is not None:
            r = self.ttv_local
            visits = ", ".join(
                f"{v.fit.offset * MINUTES_PER_DAY:+.0f}+-{v.err * MINUTES_PER_DAY:.0f}"
                for v in r.usable_visits()
            )
            logger.info(
                f"TTV diagnostics ({source} ephemeris, P={period:.5f} d): duration "
                f"{r.duration * 24:.2f} h, {len(r.usable_singles())}/{len(r.singles)} transits "
                f"timed, visit-mean O-C [min]: {visits or 'none'}"
                + (f", p_linear={r.visit_pvalue:.2g}" if np.isfinite(r.visit_pvalue) else "")
            )
        return self.ttv_local

    def _store_ttv_local(self):
        r = self.ttv_local
        vals = {"ttv_local_source": "" if r is None else r.source}
        if r is not None:
            vs, ss = r.usable_visits(), r.singles
            vals.update(
                ttv_local_period=r.period,
                ttv_local_t0=r.t0,
                ttv_local_duration=r.duration,
                ttv_local_best_offset=r.best_common_offset,
                ttv_local_visit_epochs=np.array([np.mean(v.epochs) for v in vs], dtype=float),
                ttv_local_visit_oc=np.array([v.fit.offset for v in vs], dtype=float),
                ttv_local_visit_err=np.array([v.err for v in vs], dtype=float),
                ttv_local_visit_snr=np.array([v.fit.snr for v in vs], dtype=float),
                ttv_local_visit_chi2=r.visit_chi2,
                ttv_local_visit_pvalue=r.visit_pvalue,
                ttv_local_epochs=np.array([s.epoch for s in ss], dtype=float),
                ttv_local_oc=np.array([s.fit.offset for s in ss], dtype=float),
                ttv_local_err=np.array([s.err for s in ss], dtype=float),
                ttv_local_snr=np.array([s.fit.snr for s in ss], dtype=float),
                ttv_local_usable=np.array([s.fit.usable for s in ss], dtype=bool),
            )
        for key, val in vals.items():
            self.tls_results[key] = val

    def _save_ttv_figures(self, fp) -> list:
        """Write ``<fp>_ttv.png`` and ``<fp>_ttv_transits.png`` when diagnostics ran."""
        r = self.ttv_local
        if r is None:
            return []
        title = getattr(self, "target_name", "")
        out = []
        for suffix, make in (
            (TTV_FIGURE_SUFFIX, lambda: plot_ttv_diagnostics(r, title)),
            (TTV_GALLERY_SUFFIX, lambda: plot_transit_gallery(r, *self._ttv_raw, title)),
        ):
            fig = make()
            path = Path(str(fp) + suffix)
            fig.savefig(path, dpi=120, bbox_inches="tight")
            pl.close(fig)
            logger.info(f"Saved: {path}")
            out.append(path)
        return out

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
            "ttv_beta": float(self.ttv_beta),
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
                f"{stats.amplitude * MINUTES_PER_DAY:.1f} min, beta={self.ttv_beta:.2f}, "
                f"p={stats.pvalue:.2g} "
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
        verdict = {
            "detected": "TTVs detected",
            "rejected (vetting)": "TTVs rejected by vetting",
        }.get(self._ttv_verdict_label(), "linear ephemeris OK")
        plot_oc(self.ttv_stats, ax=axes[7], n_alternatives=len(alternatives), verdict=verdict)

    def _ttv_verdict_label(self) -> str:
        if self.ttv_detected:
            return "detected"
        if self.ttv_stats is not None and self.ttv_stats.detected:
            return "rejected (vetting)"
        return "not significant"

    def _ttv_candidate_rows(self) -> list:
        stats = self.ttv_stats
        eph = stats.ephemeris
        rows = [
            ("TTV source", _TTV_SOURCE_LABELS.get(self.ttv_source, self.ttv_source)),
            ("Mean period", f"{eph.period:.5f}" + r"$\pm$" + f"{eph.period_err:.5f} d"),
            ("Transits timed", f"{stats.n}"),
            ("TTV amplitude", f"{stats.amplitude * MINUTES_PER_DAY:.1f} min"),
            ("TTV p-value", f"{stats.pvalue:.2g} (robust {stats.pvalue_robust:.2g})"),
            ("TTVs", self._ttv_verdict_label()),
        ]
        if len(self.ttv_cycle_solutions) > 1:
            rows.append(("Cycle solutions", f"{len(self.ttv_cycle_solutions)}"))
        return rows


__all__ = ["TTVMixin"]
