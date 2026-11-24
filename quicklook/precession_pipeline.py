"""Precession search glue for ``TessQuickLook`` (#63; kept out of tql.py for size).

``PrecessionMixin`` expects the host to provide ``flat_lc``, ``tls_results``,
``verbose``, ``precession_search`` and ``_sector_spans`` (recorded when several
sectors are stitched), plus ``TTVMixin._custom_ephemeris``,
``SingleEventMixin._tls_is_periodic`` and ``_stellar_prior_kwargs``.

``precession_mode`` picks the branch: ``planet`` (per-sector ``b``),
``eb`` (per-sector inner-binary inclination, ``quicklook.precession_eb``) or
``auto`` (``eb`` when SIMBAD classifies the target as an eclipsing binary).
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Optional, Union

import matplotlib.pyplot as pl
import numpy as np
from loguru import logger

from quicklook.plot_precession import plot_eb_precession, plot_precession
from quicklook.precession import (
    MIN_SECTORS,
    PrecessionResult,
    a_rs_from_star,
    search_precession,
)
from quicklook.precession_eb import EBPrecessionResult, search_eb_precession

PRECESSION_MODES = ("auto", "planet", "eb")
_EB_OTYPE = re.compile(r"eclbin|eclipsing binary|^eb\*", re.IGNORECASE)
AnyPrecessionResult = Union[PrecessionResult, EBPrecessionResult]

PRECESSION_FIGURE_SUFFIX = "_precession.png"
_SECTOR_ARRAY_FIELDS = {
    "precession_sectors": "sector",
    "precession_t_mid": "t_mid",
    "precession_b": "b",
    "precession_b_err": "b_err",
    "precession_rp_rs": "rp_rs",
    "precession_rp_rs_err": "rp_rs_err",
    "precession_depth": "depth",
    "precession_depth_err": "depth_err",
    "precession_t14": "t14",
    "precession_t14_err": "t14_err",
    "precession_dt0": "dt0",
    "precession_cadence": "cadence",
    "precession_crowdsap": "crowdsap",
    "precession_flfrcsap": "flfrcsap",
    "precession_corr_b_rp": "corr_b_rp",
}
_SCALAR_FIELDS = (
    "a_rs",
    "a_rs_err",
    "db_dt",
    "db_dt_err",
    "significance",
    "delta_bic",
    "dt14_dt",
    "dt14_dt_err",
    "drp_dt",
    "drp_dt_err",
    "baseline_yr",
    "di_dt",
    "years_to_exit",
)
_EB_SECTOR_ARRAY_FIELDS = {
    "precession_sectors": "sector",
    "precession_t_mid": "t_mid",
    "precession_inc": "inc",
    "precession_inc_err": "inc_err",
    "precession_b": "b",
    "precession_b_err": "b_err",
    "precession_depth": "depth",
    "precession_depth_err": "depth_err",
    "precession_depth_data": "depth_data",
    "precession_depth_data_err": "depth_data_err",
    "precession_t14": "t14",
    "precession_t14_err": "t14_err",
    "precession_dt0": "dt0",
    "precession_cadence": "cadence",
    "precession_crowdsap": "crowdsap",
    "precession_flfrcsap": "flfrcsap",
    "precession_has_secondary": "has_secondary",
}
_EB_SCALAR_FIELDS = (
    "di_dt",
    "di_dt_err",
    "db_dt",
    "db_dt_err",
    "significance",
    "delta_bic",
    "ddepth_dt",
    "ddepth_dt_err",
    "baseline_yr",
    "max_nodal_period_yr",
)
_EB_BINARY_FIELDS = (
    "r_sum",
    "r_sum_err",
    "k",
    "k_err",
    "sb_ratio",
    "sb_ratio_err",
    "dt_sec",
    "dt_sec_err",
    "third_light",
)


def is_eclipsing_binary_otype(otype: Optional[str]) -> bool:
    """True when a SIMBAD object type (code or description) names an eclipsing binary."""
    return bool(otype) and _EB_OTYPE.search(str(otype)) is not None


class PrecessionMixin:
    """Per-sector transit-shape fits and the b(t) trend test."""

    flat_lc: Any
    tls_results: Any
    verbose: bool

    precession_search: bool = False
    precession_mode: str = "auto"
    precession_result: Optional[AnyPrecessionResult] = None
    precession_a_rs_prior: Optional[tuple] = None
    _sector_spans: tuple = ()

    def run_precession_search(self) -> Optional[AnyPrecessionResult]:
        """Fit each sector's transit (or eclipse) shape and test for a secular drift.

        Needs ``--precession-search`` and a multi-sector light curve. Uses the
        custom ephemeris when given, else a periodic TLS signal. The EB branch
        needs the binary period with ``T0`` at the primary eclipse.
        """
        if not getattr(self, "precession_search", False):
            return None
        spans = tuple(getattr(self, "_sector_spans", ()) or ())
        if len(spans) < MIN_SECTORS:
            logger.warning(
                "Precession search needs >= 2 stitched sectors "
                "(use --sector all or a sector list); skipping."
            )
            return self._finish_precession(PrecessionResult((), None, reason="single sector"))
        ephem = self._precession_ephemeris()
        if ephem is None:
            if self.verbose:
                logger.info("Precession search: no periodic signal or custom ephemeris to fit.")
            return self._finish_precession(PrecessionResult((), None, reason="no periodic signal"))
        period, t0, duration = ephem
        lc = self.flat_lc
        flux_err = None if lc.flux_err is None else np.asarray(lc.flux_err.value, dtype=float)
        data = (
            np.asarray(lc.time.value, dtype=float),
            np.asarray(lc.flux.value, dtype=float),
            flux_err,
            spans,
            period,
            t0,
            duration,
        )
        if self.resolve_precession_mode() == "eb":
            logger.info(f"Precession search: eclipsing-binary branch (P = {period:.6f} d)")
            return self._finish_precession(search_eb_precession(*data))
        self.precession_a_rs_prior = self._precession_a_rs_prior(period)
        result = search_precession(*data, a_rs_prior=self.precession_a_rs_prior)
        return self._finish_precession(result)

    def resolve_precession_mode(self) -> str:
        """``planet`` or ``eb`` from ``precession_mode`` (``auto`` asks SIMBAD)."""
        mode = str(getattr(self, "precession_mode", "auto") or "auto").lower()
        if mode not in PRECESSION_MODES:
            raise ValueError(f"precession_mode must be one of {PRECESSION_MODES}, got {mode!r}")
        if mode != "auto":
            return mode
        otype = getattr(self, "simbad_obj_type", None)
        return "eb" if is_eclipsing_binary_otype(otype) else "planet"

    def _precession_ephemeris(self) -> Optional[tuple]:
        custom = self._custom_ephemeris() if hasattr(self, "_custom_ephemeris") else None
        if custom is not None:
            return tuple(float(v) for v in custom)
        if not (hasattr(self, "_tls_is_periodic") and self._tls_is_periodic()):
            return None
        res = self.tls_results
        return float(res["period"]), float(res["T0"]), float(res["duration"])

    def _precession_a_rs_prior(self, period: float) -> Optional[tuple]:
        """``(a/Rs, sigma)`` from the ExoFOP R*, M*, else None (a/Rs left free)."""
        if not hasattr(self, "_stellar_prior_kwargs"):
            return None
        try:
            kw = self._stellar_prior_kwargs(purpose="the precession a/Rs prior")
        except (AttributeError, KeyError, TypeError, ValueError) as exc:
            logger.warning(f"Precession search: no stellar parameters ({exc}); a/Rs is free.")
            return None
        r, m = kw.get("R_star", np.nan), kw.get("M_star", np.nan)
        r_err = kw.get("R_star_max", np.nan) - r
        m_err = kw.get("M_star_max", np.nan) - m
        return a_rs_from_star(period, r, m, r_err, m_err)

    def _finish_precession(self, result: AnyPrecessionResult) -> AnyPrecessionResult:
        self.precession_result = result
        if self.tls_results is not None:
            self._store_precession_results()
        self._log_precession_results()
        return result

    def _store_precession_results(self):
        r = self.precession_result
        is_eb = isinstance(r, EBPrecessionResult)
        values = {
            "precession_search": True,
            "precession_mode": "eb" if is_eb else "planet",
            "precession_candidate": bool(r.is_candidate),
            "precession_flags": ",".join(r.flags),
            "precession_reason": r.reason,
        }
        if is_eb:
            scalars, arrays = _EB_SCALAR_FIELDS, _EB_SECTOR_ARRAY_FIELDS
            for name in _EB_BINARY_FIELDS:
                value = getattr(r.binary, name) if r.binary is not None else np.nan
                values[f"precession_eb_{name}"] = float(value)
        else:
            scalars, arrays = _SCALAR_FIELDS, _SECTOR_ARRAY_FIELDS
            values["precession_a_rs_prior"] = bool(self.precession_a_rs_prior is not None)
        for name in scalars:
            values[f"precession_{name}"] = float(getattr(r, name))
        for key, attr in arrays.items():
            values[key] = np.array([getattr(f, attr) for f in r.sector_fits], dtype=float)
        for key, val in values.items():
            self.tls_results[key] = val

    def _log_precession_results(self):
        r = self.precession_result
        msg = f"Precession search: {len(r.sector_fits)} sector(s) fitted"
        if r.models is not None and isinstance(r, EBPrecessionResult):
            msg += (
                f" over {r.baseline_yr:.2f} yr, di/dt={r.di_dt:.3f}±{r.di_dt_err:.3f} deg/yr "
                f"({r.significance:.1f} sigma, dBIC={r.delta_bic:.1f})"
            )
            if np.isfinite(r.max_nodal_period_yr):
                msg += f", P_nodal <= {r.max_nodal_period_yr:.0f} yr"
        elif r.models is not None:
            msg += (
                f" over {r.baseline_yr:.2f} yr, db/dt={r.db_dt:.4f}±{r.db_dt_err:.4f}/yr "
                f"({r.significance:.1f} sigma, dBIC={r.delta_bic:.1f})"
            )
        msg += f" -> {'CANDIDATE' if r.is_candidate else r.reason}"
        if r.flags:
            msg += f" [flags: {', '.join(r.flags)}]"
        logger.info(msg)

    def _save_precession_figure(self, fp) -> list:
        """Write ``<fp>_precession.png`` when the search produced per-sector fits."""
        r = self.precession_result
        if r is None or not r.sector_fits:
            return []
        plot = plot_eb_precession if isinstance(r, EBPrecessionResult) else plot_precession
        fig = plot(r, getattr(self, "target_name", ""))
        path = Path(str(fp) + PRECESSION_FIGURE_SUFFIX)
        fig.savefig(path, dpi=120, bbox_inches="tight")
        pl.close(fig)
        logger.info(f"Saved: {path}")
        return [path]
