"""
TODO:
* Add momentum dumps as in TESSLatte:
https://github.com/noraeisner/LATTE/blob/7ac35c8a51949345bc076fd30a456e74fce70c51/LATTE/LATTEutils.py#L3501C13-L3501C63
"""

import math
import re
import sys
import traceback
import textwrap
import warnings
from pathlib import Path
from time import time as timer
from loguru import logger
from importlib.resources import files
from quicklook.exceptions import NoDataError, InvalidInputError
import matplotlib.pyplot as pl
import numpy as np
import pandas as pd
from transitleastsquares import transitleastsquares as tls
from wotan import flatten
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clip
from astropy.wcs import WCS
import astropy.units as u
from astroquery.simbad import Simbad
from astroquery.mast import Catalogs
import lightkurve as lk
from quicklook import h5io
from quicklook.notch_locor import notch_flatten, locor_flatten
from quicklook.long_period import default_period_ceiling, epoch_check, needs_standard_fallback
from quicklook.broken_axis import finish_broken_axes
from quicklook.single_event_pipeline import SingleEventMixin
from quicklook.ttv_pipeline import TTVMixin
from quicklook.utils import (
    get_exofop_json,
    get_params_from_exofop,
    TESS_TIME_OFFSET,
    TESS_pix_scale,
)
from quicklook.gls import Gls
from quicklook.plot import (
    _get_frac_depth,
    use_style,
    get_dss_data,
    plot_gaia_sources_on_survey,
    plot_gaia_sources_on_tpf,
    plot_odd_even_transit,
    plot_secondary_eclipse,
    plot_tls,
    plot_periodogram,
    plot_gls_periodogram,
)
from quicklook.inject import InjectionParams, run_grid, select_best_window
from quicklook.tglc import get_tglc_lc, tglc_quality_mask

# FITSFixedWarning: 'datfix' made the change 'Invalid time in DATE-OBS
warnings.filterwarnings("ignore", category=Warning, message=".*datfix.*")
warnings.filterwarnings("ignore", category=Warning, message=".*obsfix.*")

# Pipeline registry is centralized in quicklook.pipelines; these
# re-exports preserve the old import path for callers (and notebooks).
from quicklook.pipelines import (  # noqa: E402
    ALL_TESS_PIPELINES,
    FULL_FRAME_TESS_PIPELINES,
)


DATA_PATH = files("quicklook").joinpath("data")
simbad_obj_list_file = Path(DATA_PATH, "simbad_obj_types.csv")
use_style("science")

# Flatten methods handled by quicklook.notch_locor instead of wotan.
NOTCH_LOCOR_METHODS = ("notch", "locor")

__all__ = ["TessQuickLook"]


class _TLSResult(dict):
    """TLS-compatible mapping used to normalize GTLS results."""

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError as exc:
            raise AttributeError(name) from exc


def _gpu_device_count():
    """Return the number of CUDA devices visible to CuPy."""
    import cupy

    return cupy.cuda.runtime.getDeviceCount()


def _get_gpu_tls():
    """Return GTLS when it is installed and a CUDA GPU is visible."""
    try:
        device_count = _gpu_device_count()
    except Exception as exc:
        logger.debug(f"GPU TLS unavailable ({exc}); using CPU TLS")
        return None

    if device_count < 1:
        logger.info("No CUDA GPU detected; using CPU TLS")
        return None

    try:
        from gputls import gtls
    except Exception as exc:
        logger.warning(f"CUDA GPU detected but GTLS could not be loaded ({exc}); using CPU TLS")
        return None

    return gtls


def _adapt_gtls_result(result, model):
    """Expose a GTLS result through the mapping/attribute API QuickLook uses."""
    values = vars(result).copy()
    fractional_depth = float(values["depth"])
    values["depth"] = fractional_depth
    values["rp_rs"] = np.sqrt(max(0.0, fractional_depth))
    values["odd_even_mismatch"] = np.nan

    periods = np.asarray(values["periods"])
    power = np.asarray(values["power"])
    try:
        peak = int(np.nanargmax(power))
        half_max = power[peak] / 2
        lower = np.flatnonzero(power[:peak] <= half_max)
        upper = np.flatnonzero(power[peak + 1 :] <= half_max)
        if len(lower) and len(upper):
            lo = lower[-1]
            hi = peak + 1 + upper[0]
            values["period_uncertainty"] = 0.5 * (periods[hi] - periods[lo])
        else:
            values["period_uncertainty"] = np.inf
    except (ValueError, IndexError):
        values["period_uncertainty"] = np.inf

    model_flux, model_phase, _ = model.showFit()
    values["model_folded_phase"] = np.asarray(model_phase)
    values["model_folded_model"] = np.asarray(model_flux)
    return _TLSResult(values)


def format_sector_ranges(sectors):
    """Compact sorted sector list: [66, 67, 68, 93] -> "66-68+93"."""
    runs = []
    for s in sorted({int(s) for s in sectors}):
        if runs and s == runs[-1][1] + 1:
            runs[-1][1] = s
        else:
            runs.append([s, s])
    return "+".join(str(a) if a == b else f"{a}-{b}" for a, b in runs)


def _download_bitmask(author: str, quality_bitmask):
    """Bitmask to hand lightkurve when downloading a product from ``author``.

    lightkurve's TGLC reader overflows on ``'hardest'`` (int16 ``tess_flags``;
    see LIGHTKURVE_TGLC_BUG_REPORT.txt), so TGLC HLSP files are read unmasked
    and :func:`_apply_tglc_hlsp_mask` applies the requested bitmask instead.
    """
    return "none" if str(author).upper() == "TGLC" else quality_bitmask


def _apply_tglc_hlsp_mask(lc, quality_bitmask):
    """Drop cadences of a MAST TGLC light curve flagged under ``quality_bitmask``."""
    good = tglc_quality_mask(lc, quality_bitmask)
    if good.all():
        return lc
    logger.info(
        f"Masked {int((~good).sum())}/{len(lc)} TGLC cadences "
        f"(quality_bitmask={quality_bitmask!r} + TGLC background flag)."
    )
    return lc[good]


def _count_valid_transits(counts, depths) -> int:
    """Transits with in-transit cadences and a finite depth.

    TLS returns scalars (NaN) instead of arrays for degenerate fits, so both
    inputs are promoted to 1-D first.
    """
    counts = np.atleast_1d(np.asarray(counts, dtype=float))
    depths = np.atleast_1d(np.asarray(depths, dtype=float))
    n = min(len(counts), len(depths))
    return int(np.sum((counts[:n] > 0) & np.isfinite(depths[:n])))


_TGLC_FLUX_COLUMNS = {"aperture": "cal_aper_flux", "psf": "cal_psf_flux"}


def _has_finite_flux(lc) -> bool:
    return bool(np.isfinite(np.asarray(lc.flux.value, dtype=float)).any())


def _select_tglc_flux(lc, flux_type):
    """Pick the TGLC HLSP flux column; fall back to aperture if it is all NaN.

    Returns ``(lc, used)`` where ``used`` is "psf", "aperture" or the original
    ``flux_type`` when the default column was kept. TGLC leaves PSF flux NaN
    for some stars (e.g. HATS-67 s10), which would otherwise yield an empty
    light curve.
    """
    used = flux_type
    col = _TGLC_FLUX_COLUMNS.get(flux_type)
    if col is not None and col in lc.colnames:
        lc = lc.select_flux(col)
    if _has_finite_flux(lc) or "cal_aper_flux" not in lc.colnames:
        return lc, used
    logger.warning(
        f"TGLC {flux_type!r} flux is all NaN for this target; using aperture photometry."
    )
    return lc.select_flux("cal_aper_flux"), "aperture"


def _aperture_if_flux_all_nan(lc):
    """Local ePSF light curve: swap all-NaN PSF flux for aperture flux.

    The local extraction derives ``flux_err`` from the selected flux, so the
    error is recomputed (MAD scatter) for the aperture flux as well. Seen for
    stars near the cutout edge (Kepler-1892, Kepler-501 s14).
    """
    if _has_finite_flux(lc) or "cal_aper_flux" not in lc.colnames:
        return lc
    aper = np.asarray(lc["cal_aper_flux"], dtype=float)
    if not np.isfinite(aper).any():
        return lc
    logger.warning("TGLC PSF flux is all NaN for this target; using aperture photometry.")
    sector = getattr(lc, "sector", None)
    out = lc.copy()
    out["flux"] = aper * out.flux.unit
    mad = 1.4826 * np.nanmedian(np.abs(aper - np.nanmedian(aper)))
    out["flux_err"] = np.full_like(aper, mad) * out.flux.unit
    out.meta["FLUX_ORIGIN"] = "cal_aper_flux"
    out.sector = sector
    return out


def _require_usable_flux(lc):
    """Raise a clear error when the local ePSF extraction gave no photometry.

    Stars on the CCD edge (Kepler-1892, Kepler-501 s14) come back with all-NaN
    PSF flux and a constant aperture flux; TLS would otherwise fail later with
    an opaque ``max() arg is an empty sequence``.
    """
    flux = np.asarray(lc.flux.value, dtype=float)
    finite = flux[np.isfinite(flux)]
    if finite.size < 3 or np.ptp(finite) == 0:
        raise ValueError(
            f"Local TGLC extraction for sector {getattr(lc, 'sector', '?')} produced no "
            "usable photometry (PSF flux all NaN, aperture flux constant); the target "
            "is likely on the CCD edge. Try another sector or pipeline."
        )


class TessQuickLook(SingleEventMixin, TTVMixin):
    def __init__(
        self,
        target_name: str,
        sector=-1,
        pipeline: str = "SPOC",
        flux_type="pdcsap",
        exptime: float = None,
        pg_method: str = "gls",
        flatten_method: str = "biweight",
        gp_kernel: str = "periodic_auto",  # works if flatten_method=='gp'
        gp_kernel_size: float = 5,  # works if flatten_method=='gp'
        # bin_size : float = None, # useful for dense lc (exp<=120s)
        window_length: float = None,
        edge_cutoff: float = 0.1,
        sigma_clip_raw: tuple = None,
        sigma_clip_flat: tuple = None,
        custom_ephem: list = None,
        mask_ephem: bool = False,
        Porb_limits: tuple = None,
        extend_period_ceiling: bool = True,
        archival_survey="dss1",
        show_simbad: bool = False,
        show_plot: bool = True,
        verbose: bool = True,
        savefig: bool = False,
        savetls: bool = False,
        overwrite: bool = False,
        suffix: str = None,
        quality_bitmask: str = "hard",
        outdir: str = ".",
        tls_use_threads: int = None,
        use_star_priors: bool = False,
        phase_xlim: float = None,
        iterative_search: bool = False,
        min_sde_iterative: float = 5.0,
        max_planets: int = None,
        cancel_check=None,
        neighbor_check: bool = True,
        single_event_search: bool = False,
        ttv_search: bool = False,
    ):
        # start timer
        self.timer_start = timer()
        self.target_name = target_name
        logger.info(f"Generating quicklook for {self.target_name}...")
        self.verbose = verbose
        self.show_plot = show_plot
        self.exofop_data = get_exofop_json(target_name)
        self.parse_exofop_info()
        self.custom_ephem = custom_ephem
        self.parse_custom_ephem()
        self.simbad_obj_type = self.get_simbad_obj_type()
        self.flux_type = flux_type
        self.exptime = exptime
        self.pg_method = pg_method
        self.sigma_clip_raw = sigma_clip_raw
        self.sigma_clip_flat = sigma_clip_flat
        self.quality_bitmask = quality_bitmask
        self.flatten_method = flatten_method
        self.tls_use_threads = tls_use_threads
        self.use_star_priors = use_star_priors
        self.phase_xlim = phase_xlim
        self.iterative_search = iterative_search
        self.min_sde_iterative = min_sde_iterative
        self.max_planets = max_planets
        self.single_event_search = single_event_search
        self.ttv_search = ttv_search
        # Optional zero-arg callable polled inside the long TGLC ePSF loop;
        # returning truthy raises InterruptedError so the worker thread
        # actually stops when the GUI Cancel button is clicked.
        self.cancel_check = cancel_check
        # False skips the TESS TPF/Gaia/DSS neighbor panel (e.g. non-TESS light curves).
        self.neighbor_check = neighbor_check
        self._check_cancel("before light-curve download")
        self.raw_lc = self.get_lc(
            author=pipeline,
            sector=sector,
            exptime=self.exptime,  # cadence=cadence
        )
        self._check_cancel("after light-curve download")
        self.overwrite = overwrite
        self.outdir = outdir
        self.suffix = suffix
        self.mask_ephem = mask_ephem
        _ = self.check_output_file_exists()
        self.gp_kernel = gp_kernel  # squared_exp, matern, periodic, periodic_auto
        self.gp_kernel_size = gp_kernel_size
        self.edge_cutoff = edge_cutoff

        if window_length is None and flatten_method == "locor":
            # window_length doubles as LOCoR's rotation period, which the
            # toi_dur/0.5 d heuristic below doesn't make sense for; defer to
            # the GLS-based period estimate computed in flatten_raw_lc().
            self._use_rotation_aware_window = False
            self.window_length = 0.0
        elif window_length is None:
            self._use_rotation_aware_window = True
            self.window_length = (
                self.toi_dur[0] * 3
                if (self.toi_dur is not None) and (self.toi_dur[0] * 3 >= 0.1)
                else 0.5
            )
        else:
            self._use_rotation_aware_window = False
            self.window_length = window_length

        self.tmask = self.get_transit_mask()
        # A known ephemeris that predicts no transit within the downloaded
        # data is expected for long-period planets whose transit simply does
        # not fall in this ~27-day sector (e.g. TOI-2074, P=177.6 d). This is
        # not a failure: flattening, the TLS search, and the plots still run.
        # Warn (pointing at the nearest predicted transit) and continue rather
        # than aborting the whole quicklook.
        if self.toi_epoch is not None and self.tmask.sum() == 0:
            self._warn_no_transits_in_sector()
        if self.mask_ephem and self.tmask.sum() > 0:
            if self.verbose:
                logger.info(
                    f"Masking transits in raw lightcurve using {self.ephem_source} ephem..."
                )
            self.raw_lc = self.raw_lc[~self.tmask]
            # update tmask
            self.tmask = self.get_transit_mask()
        self._check_cancel("before flattening light curve")
        # self.flat_lc, self.trend_lc = self.raw_lc.flatten(return_trend=True)
        self.flat_lc, self.trend_lc = self.flatten_raw_lc()
        self._check_cancel("after flattening light curve")
        self.Porb_min = 0.1 if Porb_limits is None else Porb_limits[0]
        flat_time = self.flat_lc.time.value
        self._flat_time_min = flat_time.min()
        self._flat_time_max = flat_time.max()
        self.extend_period_ceiling = extend_period_ceiling
        self.Porb_max = (
            default_period_ceiling(self._flat_time_max - self._flat_time_min, extend_period_ceiling)
            if Porb_limits is None
            else Porb_limits[1]
        )
        self._check_cancel("before TLS transit search")
        self.run_tls()
        self.run_epoch_check()
        self._check_cancel("after TLS transit search")
        self.fold_lc = self.flat_lc.fold(
            period=self.tls_results.period,
            epoch_time=self.tls_results.T0,
            normalize_phase=False,
            wrap_phase=self.tls_results.period / 2,
        )
        # compute_advanced_vetting_metrics() already ran once inside run_tls(),
        # but self.fold_lc didn't exist yet at that point, so secondary_depth
        # and secondary_sde (which need it) were always computed against
        # fold_lc=None and left as NaN on the primary candidate regardless of
        # whether the data could support them. Re-run now that fold_lc is
        # available; depth_variance_ratio/duration_consistency_ratio (which
        # don't need fold_lc) just recompute to the same values.
        self.compute_advanced_vetting_metrics()
        # Before the iterative search, so it can mask transits at their fitted times.
        self.run_ttv_search()
        self._check_cancel("before iterative TLS search")
        if self.iterative_search:
            if self.verbose:
                logger.info("Running iterative multi-planet TLS search...")
            self.run_iterative_tls()
        else:
            self.iterative_tls_results = [self.tls_results] if self.tls_results is not None else []
        self._check_cancel("after iterative TLS search")
        self.run_single_event_search()
        self.run_ttv_from_chain()
        self.run_ttv_diagnostics()
        # if self.fold_lc.primary_key is None:
        #     self.fold_lc.primary_key = ('time',)
        self.savefig = savefig
        self.savetls = savetls
        self.archival_survey = archival_survey
        self.show_simbad = show_simbad

    def _check_cancel(self, phase=""):
        """Raise InterruptedError if the GUI's cancel button was pressed.

        Called at phase boundaries so a cancel issued mid-pipeline (e.g.
        during a slow MAST download or a multi-minute TLS run) is honoured
        as soon as the current blocking call returns, instead of running
        the rest of the pipeline to completion.
        """
        if self.cancel_check is not None and self.cancel_check():
            raise InterruptedError(f"Job cancelled ({phase})" if phase else "Job cancelled")

    def __repr__(self):
        """Override to print a readable string representation of class.

        This is mainly used for debugging and logging purposes.
        """
        included_args = [
            # ===target attributes===
            "name",
            "search_radius",
            "sector",
            "exptime",
            "mission",
            "campaign",
            # "all_sectors",
            # ===tpf attributes===
            "sap_mask",
            "quality_bitmask",
            # 'aper_radius', 'threshold_sigma', 'percentile'
            # cutout_size #for FFI
            "pipeline",
        ]
        # Get the values of the included args
        args = []
        for key in self.__dict__.keys():
            val = self.__dict__.get(key)
            if key in included_args:
                if key == "target_coord":
                    # Format the coordinate string
                    coord = self.target_coord.to_string("decimal")
                    args.append(f"{key}=({coord.replace(' ', ',')})")
                elif val is not None:
                    args.append(f"{key}={val}")
        # Join the args with commas
        args = ", ".join(args)
        return f"{type(self).__name__}({args})"

    def parse_exofop_info(self):
        """
        Parse the TFOP info to get the star names,
        Gaia name, Gaia ID, and target coordinates.
        """
        self.star_names = np.array(self.exofop_data.get("basic_info")["star_names"].split(", "))
        if self.verbose:
            names_block = "\n".join(f"\t{n}" for n in self.star_names)
            logger.info(f"Catalog names:\n{names_block}")
        self.gaia_name = self.star_names[
            np.array([i[:4].lower() == "gaia" for i in self.star_names])
        ][0]
        self.gaiaid = int(self.gaia_name.split()[-1])
        ra, dec = (
            self.exofop_data.get("coordinates")["ra"],
            self.exofop_data.get("coordinates")["dec"],
        )
        self.target_coord = SkyCoord(ra=ra, dec=dec, unit="degree")

        if self.target_name.lower()[:3] == "toi":
            self.toiid = int(float(self.target_name.split("-")[-1]))
        else:
            idx = [i[:3].lower() == "toi" for i in self.star_names]
            if sum(idx) > 0:
                toiid = int(self.star_names[idx][0].split("-")[-1])
            else:
                toiid = None
            self.toiid = toiid
        self.ticid = int(self.exofop_data.get("basic_info")["tic_id"])
        if self.ticid is not None:
            self.query_name = f"TIC{self.ticid}"
        else:
            self.query_name = self.target_name.replace("-", " ")

    def parse_custom_ephem(self):
        """
        Parse the custom ephem to get the ephemeris source and parameters.
        """
        if self.custom_ephem:
            # Copy: typer passes a tuple, and the unit conversions below work in place.
            self.custom_ephem = [float(x) for x in self.custom_ephem]
            self.ephem_source = "custom"
            err_msg = "Custom ephem must be a tuple: (t0,t0err,P,Perr,t14,t14err)"
            if len(self.custom_ephem) != 6:
                raise InvalidInputError(err_msg)
            if self.verbose:
                msg = f"Using ephemeris mask:\nt0={self.custom_ephem[0]}\n"
                msg += f"P={self.custom_ephem[2]}d\nt14={self.custom_ephem[4]}"
                logger.info(msg)
            # TODO: using tfop in variable name is misleading
            if self.custom_ephem[0] > TESS_TIME_OFFSET:
                if self.verbose:
                    msg = "Custom transit epoch given in JD. "
                    msg += f"Converting to BTJD = JD-{TESS_TIME_OFFSET:,}."
                    logger.info(msg)
                self.custom_ephem[0] -= TESS_TIME_OFFSET
            self.toi_epoch = (self.custom_ephem[0], self.custom_ephem[1])
            self.toi_period = (self.custom_ephem[2], self.custom_ephem[3])
            if self.custom_ephem[4] > 1:
                if self.verbose:
                    logger.info("Custom transit duration given in hours. Converting to days.")
                self.custom_ephem[4] /= 24
                self.custom_ephem[5] /= 24
            self.toi_dur = (self.custom_ephem[4], self.custom_ephem[5])
            self.toi_depth = None
        else:
            # use tfop ephem if available
            (
                self.toi_epoch,
                self.toi_period,
                self.toi_dur,
                self.toi_depth,
            ) = (
                self.get_toi_ephem()
                if len(self.exofop_data.get("planet_parameters")) != 0
                else (None, None, None, None)
            )
            self.ephem_source = "TFOP" if self.toi_epoch is not None else None

    def _is_multi_sector(self):
        """True when the light curve stitches several sectors (``sector="all"``)."""
        return isinstance(self.sector, (list, tuple, np.ndarray))

    def _sector_label(self):
        """Filename token for the sector.

        Zero-padded number for one sector, "all" for a full stitch, or the
        compact range (e.g. "66-68+93") when explicit sectors were requested.
        """
        if not self._is_multi_sector():
            return str(self.sector).zfill(2)
        requested = getattr(self, "_sectors_requested", None)
        return format_sector_ranges(requested) if requested else "all"

    def _sector_display(self):
        """Human-readable sector text for plot titles."""
        if self._is_multi_sector():
            requested = getattr(self, "_sectors_requested", None)
            if requested:
                return f"{format_sector_ranges(requested).replace('+', ', ')} ({len(self.sector)} sectors)"
            return f"{min(self.sector)}-{max(self.sector)} ({len(self.sector)} sectors)"
        return str(self.sector)

    @property
    def _tpf_sector(self):
        """Single sector for the TPF/neighbor panel (latest one if stitched)."""
        return max(self.sector) if self._is_multi_sector() else self.sector

    def check_output_file_exists(self):
        name = self.target_name.replace(" ", "")
        if name.lower()[:3] == "toi":
            name = f"TOI{str(self.toiid).zfill(4)}"
        # Filename "lctype" token. For SPOC it's the flux column (pdcsap/sap);
        # for TGLC with an explicit aperture/PSF choice it's "tglc_aper" or
        # "tglc_psf" so the two products of the same target/sector don't
        # overwrite each other on disk. Otherwise it's just the pipeline name.
        if self.pipeline == "spoc":
            lctype = self.flux_type
        elif self.pipeline == "tglc" and self.flux_type in ("aperture", "psf"):
            lctype = "tglc_aper" if self.flux_type == "aperture" else "tglc_psf"
        else:
            lctype = self.pipeline
        fp = Path(
            self.outdir,
            f"{name}_s{self._sector_label()}_{lctype}_{self.cadence[0]}c",
        )
        if self.mask_ephem:
            fp = fp.with_stem(fp.stem + "_mask_ephem")
        if self.suffix is not None:
            fp = fp.with_stem(fp.stem + f"_{self.suffix}")
        png_file = fp.with_suffix(".png")
        if png_file.exists() and not self.overwrite:
            raise FileExistsError(f"{png_file} already exists! Set overwrite=True.")
        return fp

    def query_simbad(self):
        """
        Query Simbad to get the object type of the target star.

        Returns
        -------
        res : SimbadResult
            The result of the query, if the target is resolved.
            Otherwise, None.
        """
        # The object type is optional metadata, so an unreachable or failing
        # Simbad service (e.g. a TAP outage) must not abort the whole run.
        try:
            # See also: https://simbad.cds.unistra.fr/guide/otypes.htx
            Simbad.add_votable_fields("otype")

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")

                # Try resolving the target star by name
                for name in self.star_names:
                    r = Simbad.query_object(name)
                    if r is not None:
                        return r
                    if self.verbose:
                        logger.info(f"Simbad cannot resolve {name}.")
        except Exception as e:
            logger.warning(
                f"Simbad query failed ({type(e).__name__}: {e}); continuing without object type."
            )
        return None

    def get_simbad_obj_type(self):
        """
        Retrieves the object type of the target star from Simbad.

        Returns
        -------
        str or None
            The description of the object type if found, otherwise None.
        """
        # Query Simbad for the target star
        r = self.query_simbad()

        if r:
            # Extract the object type category
            category = r.to_pandas().squeeze()["otype"]

            if len(category) >= 4:
                return category

            # Load Simbad object type descriptions
            df = pd.read_csv(simbad_obj_list_file)
            dd = df.query("Id==@category")

            if len(dd) > 0:
                # Retrieve the description and id
                desc = dd["Description"].squeeze()
                oid = dd["Id"].squeeze()

                # Check if the description contains 'binary' and print appropriate message
                if dd["Description"].str.contains("(?i)binary").any():
                    logger.info("***" * 15)
                    logger.info(f"Simbad classifies {self.target_name} as {oid}={desc}!")
                    logger.info("***" * 15)
                else:
                    logger.info(f"Simbad classifies {self.target_name} as {oid}={desc}!")

                return desc
        # Return None if no object type is found
        return None

    def get_lc(self, **kwargs: dict) -> lk.TessLightCurve:
        """
        Retrieves a light curve for the specified target.

        Parameters
        ----------
        kwargs : dict
            Additional parameters such as radius, sector, author, cadence, exptime.

        Returns
        -------
        lk.TessLightCurve
            The retrieved light curve object.
        """
        # Determine the sector of interest
        requested_sectors = None
        self._sectors_requested = None
        if kwargs.get("sector") is None:
            sector_orig = None
        else:
            sector_orig = kwargs.pop("sector")
            # An explicit list of sectors is a partial "all": stitch only those.
            if isinstance(sector_orig, (list, tuple, set)):
                requested_sectors = sorted({int(s) for s in sector_orig})
                sector_orig = "all"
            # sector = None if sector_orig in ["all", -1] else int(sector_orig)
            sector = (
                None
                if str(sector_orig).lower() == "all" or int(sector_orig) == -1
                else int(sector_orig)
            )
            kwargs["sector"] = sector

        if sector_orig == "all":
            # assure data comes from single pipeline
            if kwargs.get("exptime") is None:
                raise InvalidInputError("Supply exptime when using sector='all' or a sector list.")

        # Single MAST search, then filter locally
        search_result_all_lcs = lk.search_lightcurve(self.query_name)
        err_msg = f"Search using '{self.query_name}' did not yield any lightcurve results."
        if len(search_result_all_lcs) == 0:
            if kwargs.get("author", "").upper() == "TGLC":
                return self._get_tglc_lc_fallback(kwargs.get("sector"))
            raise NoDataError(err_msg)

        # Extract all available sectors
        self.all_sectors = self.get_unique_sectors(search_result_all_lcs)
        # Display available light curves
        cols = ["author", "mission", "t_exptime"]
        if self.verbose:
            lc_table = search_result_all_lcs.table.to_pandas()[cols]
            logger.info(f"All available lightcurves:\n{lc_table.to_string()}")

        # Validate the requested sector. If the user asked for a specific
        # sector that the MAST search did not return (e.g. a stale value
        # left over in the GUI from a different target), fall back to the
        # latest available sector instead of raising — much friendlier
        # than killing the job, and the warning still surfaces in the log.
        if kwargs.get("sector") is None:
            if self.verbose:
                logger.info(f"Available sectors: {self.all_sectors}")
        else:
            if kwargs.get("sector") not in self.all_sectors:
                logger.warning(
                    f"Requested sector={kwargs.get('sector')} is not in the "
                    f"available sectors {self.all_sectors} for {self.query_name}; "
                    "falling back to the latest available sector."
                )
                sector_orig = -1
                kwargs["sector"] = None

        # Validate the requested author
        if kwargs.get("author") is None:
            kwargs["author"] = "SPOC"
        else:
            all_authors = set(search_result_all_lcs.table["provenance_name"].tolist())
            if kwargs["author"].upper() not in all_authors:
                if kwargs["author"].upper() == "TGLC":
                    return self._get_tglc_lc_fallback(kwargs.get("sector"))
                err_msg = f"author={kwargs.get('author')} not in {all_authors}"
                raise NoDataError(err_msg)
        self.pipeline = kwargs["author"].lower()
        self.all_pipelines = all_authors

        # Filter locally instead of a second MAST query
        mask = np.ones(len(search_result_all_lcs), dtype=bool)
        tbl = search_result_all_lcs.table
        if kwargs.get("author"):
            mask &= np.array(
                [a.upper() == kwargs["author"].upper() for a in tbl["provenance_name"]]
            )
        if kwargs.get("sector") is not None:
            mask &= np.array(
                [
                    len(m.split()) == 3 and int(m.split()[-1]) == kwargs["sector"]
                    for m in tbl["mission"]
                ]
            )
        # Falsy exptime (None, "", 0) means "no filter" — otherwise the
        # comparison `tbl["t_exptime"] == ""` evaluates to False on every
        # row and zeros out the search, which then falsely triggered the
        # TGLC ePSF fallback even when MAST had products.
        if kwargs.get("exptime"):
            mask &= np.array(tbl["t_exptime"]) == kwargs["exptime"]
        if requested_sectors:
            available = self.get_unique_sectors(search_result_all_lcs[mask])
            missing = [s for s in requested_sectors if s not in available]
            if missing:
                raise NoDataError(
                    f"Sectors {missing} not available for {self.query_name} with "
                    f"author={kwargs.get('author')} exptime={kwargs.get('exptime')}; "
                    f"available at this cadence: {available}."
                )
            mask &= np.array(
                [
                    len(m.split()) == 3 and int(m.split()[-1]) in requested_sectors
                    for m in tbl["mission"]
                ]
            )
            self._sectors_requested = requested_sectors
        search_result = search_result_all_lcs[mask]
        if len(search_result) == 0:
            if kwargs.get("author", "").upper() == "TGLC":
                return self._get_tglc_lc_fallback(kwargs.get("sector"))
            err_msg = f"Search using '{self.query_name}' "
            err_msg += f"{kwargs} did not yield any lightcurve results."
            raise NoDataError(err_msg)

        # Download and return light curve
        if sector_orig == "all":
            if self.verbose:
                filtered_table = search_result.table.to_pandas()[cols]
                logger.info(
                    f"Filtered lightcurves based on query ({kwargs}):\n{filtered_table.to_string()}"
                )
            msg = f"Downloading all {kwargs.get('author')} lcs..."
            if self.verbose:
                logger.info(msg)
            filtered_sectors = self.get_unique_sectors(search_result)
            if len(filtered_sectors) <= 1:
                raise NoDataError(f"Only {len(filtered_sectors)} sector is available.")
            lc = self._download_with_retry(
                lambda: self._stitch_masked(search_result, kwargs["author"])
            )
            # only the sectors actually stitched (same pipeline and exptime)
            self.sector = filtered_sectors

            if self.pipeline in ["spoc"]:
                # exptime = int(lc.meta["EXPOSURE"] / 10) * 10
                exptime = int(lc._meta["FRAMETIM"] * lc._meta["NUM_FRM"])
            else:
                # estimate exp time
                exptime = int(np.median(np.diff(lc.time.jd)) * 24 * 60 * 60)
            msg = f"Downloaded all {kwargs.get('author')} (exp={exptime} s) lcs "
            msg += f"in sectors {', '.join([str(s) for s in filtered_sectors])}."
            if self.verbose:
                logger.info(msg)
        else:
            # Download the light curve for the specified sector
            msg = f"Downloading {kwargs.get('author').upper()} lc..."
            if self.verbose:
                logger.info(msg)
            idx = sector_orig if sector_orig == -1 else 0
            lc = self._download_with_retry(
                lambda: self._download_masked(search_result[idx], kwargs["author"])
            )

            if self.pipeline in ["spoc"]:
                # exptime = int(lc.meta["EXPOSURE"] / 10) * 10
                exptime = int(lc._meta["FRAMETIM"] * lc._meta["NUM_FRM"])
            else:
                # estimate exp time
                exptime = int(np.median(np.diff(lc.time.jd)) * 24 * 60 * 60)
            msg = f"Downloaded {lc.meta['AUTHOR'].upper()} "
            msg += f"(exp={exptime} s) lc in sector {lc.sector}."
            if self.verbose:
                logger.info(msg)
            self.sector = lc.sector

        # Select flux type / photometry for pipelines that expose a choice.
        author = lc.meta["AUTHOR"].lower()
        if author == "spoc":
            lc = lc.select_flux(self.flux_type + "_flux")
        elif author == "tglc":
            # MAST TGLC HLSP files carry both calibrated flux columns.
            lc, used = _select_tglc_flux(lc, self.flux_type)
            if used != self.flux_type:
                self.flux_type = used

        # Set exposure time and cadence
        if self.exptime is None:
            self.exptime = exptime
        if (self.exptime <= 120) and (self.flatten_method == "gp"):
            raise InvalidInputError(
                "Using flatten_method='GP' for dense data with exp<=120s is not recommended."
            )
        # assert self.exptime == search_result.exptime[idx].value
        self.cadence = "short" if self.exptime < 1800 else "long"

        if self.pipeline in ["cdips"]:
            # TODO: improve mag err estimate
            lc.flux_err = np.full_like(lc.flux, 0.01 * u.mag)  # magnitude error
            err_msg = "CDIPS pipeline has no flux error column.\n"
            err_msg += "Assuming err=0.1 mag"
            logger.error(err_msg)  # show in red

        # Apply sigma clipping if specified
        if self.sigma_clip_raw is not None:
            if self.verbose:
                logger.info("Applying sigma clip on raw lc with ")
                logger.info(f"(lower,upper)={self.sigma_clip_raw}...")
            return lc.normalize().remove_outliers(
                sigma_lower=self.sigma_clip_raw[0],
                sigma_upper=self.sigma_clip_raw[1],
            )
        else:
            return lc.normalize()

    def get_unique_sectors(self, search_result):
        missions = search_result.table["mission"].tolist()
        all_sectors = [int(x.split()[-1]) for x in missions if len(x.split()) == 3]
        return sorted(set(all_sectors))

    def _download_masked(self, search_row, author):
        """Download one light curve, applying the quality bitmask safely for TGLC."""
        lc = search_row.download(quality_bitmask=_download_bitmask(author, self.quality_bitmask))
        if str(author).upper() == "TGLC":
            lc = _apply_tglc_hlsp_mask(lc, self.quality_bitmask)
        return lc

    def _stitch_masked(self, search_result, author):
        """Download and stitch several light curves, masking TGLC per sector."""
        lcs = search_result.download_all(
            quality_bitmask=_download_bitmask(author, self.quality_bitmask)
        )
        if str(author).upper() == "TGLC":
            lcs = lk.LightCurveCollection(
                [_apply_tglc_hlsp_mask(lc, self.quality_bitmask) for lc in lcs]
            )
        return lcs.stitch()

    def _download_with_retry(self, download_callable, max_retries=1):
        """Run a lightkurve download; if it fails with a "may be corrupt"
        error, diagnose the failure (FITS truncation vs. lightkurve
        adapter error), quarantine the suspect file, and retry once if
        retrying could plausibly help.

        Lightkurve's ``io/read.py`` re-labels *any* exception from its
        format-specific reader as "this file may be corrupt — remove
        it." That message is misleading: a file that passes
        ``astropy.io.fits.verify('exception')`` is not corrupt — the
        adapter just failed, and retrying with the same versions will
        fail the same way. Distinguish the two cases:

        * **truncated** — astropy verify also fails. Genuine partial
          download. Quarantine and retry.
        * **adapter_error** — astropy verify passes. Lightkurve's
          adapter is at fault (schema drift, version skew). Quarantine
          for post-mortem and raise immediately; no retry.
        * **missing** — file already gone (e.g. previous quarantine
          step ran). Treat like truncated for retry purposes.

        Quarantined files are renamed to
        ``<original>.bad-<unix_ts>`` so the on-disk evidence survives
        the next run.

        Parameters
        ----------
        download_callable : Callable
            Zero-arg callable that performs the lightkurve download.
        max_retries : int
            How many times to retry after a quarantine. Defaults to 1.
        """
        from lightkurve.utils import LightkurveError
        from astropy.io import fits
        from time import time as _now

        def _diagnose(path):
            if not path.exists():
                return "missing"
            try:
                with fits.open(path) as hdul:
                    for hdu in hdul:
                        hdu.verify("exception")
                return "adapter_error"
            except Exception:
                return "truncated"

        def _quarantine(path):
            ts = int(_now())
            dest = path.with_suffix(path.suffix + f".bad-{ts}")
            try:
                path.rename(dest)
                return dest
            except OSError as rename_err:
                logger.warning(f"Could not quarantine {path}: {rename_err}. Deleting instead.")
                try:
                    path.unlink()
                except OSError:
                    pass
                return None

        for attempt in range(max_retries + 1):
            try:
                return download_callable()
            except LightkurveError as e:
                msg = str(e)
                # Lightkurve's exact text:
                #   "Error in reading Data product <PATH> of type <TYPE> .
                #    This file may be corrupt due to an interrupted download.
                #    Please remove it from your disk and try again."
                m = re.search(r"Error in reading Data product\s+(\S+)", msg)
                if not m:
                    # Unrecognised error shape — nothing to clean up.
                    raise

                bad_path = Path(m.group(1))
                verdict = _diagnose(bad_path)

                if verdict == "adapter_error":
                    # File is structurally valid. Lightkurve's TGLC/SPOC
                    # adapter raised something the read wrapper turned
                    # into a "corrupt" message. Retrying won't help.
                    quarantined = _quarantine(bad_path)
                    where = (
                        f"Quarantined the suspect file at: {quarantined}"
                        if quarantined is not None
                        else "Quarantine failed; file removed from cache"
                    )
                    raise LightkurveError(
                        f"Lightkurve reported '{bad_path.name}' as corrupt, but "
                        f"astropy.io.fits.verify('exception') passes — the FITS "
                        f"file is structurally valid. This is an adapter-level "
                        f"failure (likely a lightkurve/astropy version skew or an "
                        f"upstream HLSP schema change), not an interrupted "
                        f"download. Retrying with the same package versions will "
                        f"fail the same way. {where}. "
                        f"Next steps: (1) inspect the quarantined file with "
                        f"astropy.io.fits.info() to confirm the schema, "
                        f"(2) check for a lightkurve update, or "
                        f"(3) for TGLC, use the local ePSF fallback "
                        f"(_get_tglc_lc_fallback) which bypasses the HLSP adapter."
                    ) from e

                # verdict in ("truncated", "missing") — quarantine and retry.
                if verdict == "truncated":
                    quarantined = _quarantine(bad_path)
                    if quarantined is not None:
                        logger.warning(
                            f"FITS verify failed (truncated/corrupt). Quarantined to {quarantined}."
                        )
                else:
                    logger.info(f"Suspect cache file already missing on disk: {bad_path}")

                if attempt >= max_retries:
                    attempts = max_retries + 1
                    raise LightkurveError(
                        f"MAST returned a truncated/corrupt {bad_path.name} on all "
                        f"{attempts} attempt(s). Suspect copies have been "
                        f"quarantined alongside the cache as "
                        f"<name>.bad-<timestamp> for inspection. "
                        f"This usually means the upstream HLSP product is "
                        f"temporarily unavailable or the network dropped mid-download. "
                        f"Try: (1) re-run the same command in a few minutes, "
                        f"(2) pick a different sector with --sector, or "
                        f"(3) for TGLC, the local ePSF fallback "
                        f"(_get_tglc_lc_fallback) can recompute the light curve "
                        f"directly from the FFI cutout."
                    ) from e
                logger.info(f"Retrying download (attempt {attempt + 2}/{max_retries + 1})...")

    def _get_tglc_lc_fallback(self, sector):
        """Run local TGLC ePSF extraction when MAST has no TGLC products."""
        logger.info("No TGLC products on MAST; running local ePSF extraction...")
        self.pipeline = "tglc"
        self.all_pipelines = {"TGLC"}
        # Pass the request through to get_tglc_lc unchanged: None -> first
        # available, -1 -> latest available, positive int -> that sector.
        sector_arg = None if sector is None else int(sector)
        # flux_type carries the GUI's aperture/psf choice for TGLC; anything
        # else (e.g. a stale "pdcsap") falls back to automatic selection.
        photometry = self.flux_type if self.flux_type in ("aperture", "psf", "auto") else "auto"
        lc = get_tglc_lc(
            self.query_name,
            sector=sector_arg,
            photometry=photometry,
            verbose=self.verbose,
            cancel_check=getattr(self, "cancel_check", None),
        )
        lc = _aperture_if_flux_all_nan(lc)
        _require_usable_flux(lc)
        good = tglc_quality_mask(lc, self.quality_bitmask)
        if not good.all():
            logger.info(
                f"Masked {int((~good).sum())}/{len(lc)} TGLC cadences "
                f"(quality_bitmask={self.quality_bitmask!r} + TGLC background flag)."
            )
        sector = lc.sector
        lc = lc[good]
        lc.sector = sector
        self.sector = lc.sector
        self.all_sectors = [lc.sector]
        self.exptime = lc.meta.get("EXPOSURE", 1800)
        self.cadence = "short" if self.exptime < 1800 else "long"
        if self.sigma_clip_raw is not None:
            return lc.normalize().remove_outliers(
                sigma_lower=self.sigma_clip_raw[0],
                sigma_upper=self.sigma_clip_raw[1],
            )
        return lc.normalize()

    def get_tpf(self, **kwargs: dict) -> lk.targetpixelfile.TargetPixelFile:
        """
        Search for, download, and return a TPF file.

        Parameters
        ----------
        sector: int or str
            TESS sector number or "all"
        author: str
            Pipeline author, e.g. "QLP" or "SPOC"
        """
        if kwargs.get("sector") is None:
            sector_orig = None
        else:
            sector_orig = kwargs.pop("sector")
            sector = None if sector_orig in ["all", -1] else int(sector_orig)
            kwargs["sector"] = sector

        if kwargs.get("author") is None:
            kwargs["author"] = "SPOC"

        # Search for TPF files
        search_result_all_tpfs = lk.search_targetpixelfile(self.query_name)
        if len(search_result_all_tpfs) == 0:
            raise NoDataError("No TPF files found.")

        cols = ["author", "mission", "t_exptime"]
        if self.verbose:
            tpf_table = search_result_all_tpfs.table.to_pandas()[cols]
            logger.info(f"All available TPFs:\n{tpf_table.to_string()}")
        tpf_authors = search_result_all_tpfs.table.to_pandas()["author"].unique()
        if kwargs.get("author").upper() not in tpf_authors:
            if self.verbose:
                logger.error(f"No TPF for {kwargs.get('author').upper()} pipeline.")
            kwargs["author"] = [tpf_authors[i] for i in range(len(tpf_authors)) if i != "QLP"][0]
        if self.verbose:
            logger.info(f"Using {kwargs.get('author').upper()} TPF...")

        # Search using the specified author and sector
        search_result = lk.search_targetpixelfile(self.query_name, **kwargs)
        if len(search_result) == 0:
            raise NoDataError(
                f"Search using '{self.query_name}' {kwargs} did not yield any TPF results."
            )
        msg = "Downloading TPF..."
        if self.verbose:
            logger.info(msg)
        idx = sector_orig if sector_orig == -1 else 0
        tpf = self._download_with_retry(
            lambda: search_result[idx].download(quality_bitmask=self.quality_bitmask)
        )
        # FIXME: What is the correct tpf aperture for other pipeline?
        # author = tpf.meta['PROCVER'].split('-')[0]
        author = search_result.author[idx].upper()
        exptime = search_result.exptime[idx].value
        msg = f"Downloaded {author} (exp={exptime} s) TPF "
        msg += f"in sector {tpf.meta['SECTOR']}."
        if self.verbose:
            logger.info(msg)
        return tpf

    def _load_tpf(self) -> bool:
        """Set ``self.tpf`` / ``self.sap_mask`` for the neighbor-check panels.

        Returns False when no TPF can be had: for FFI pipelines both TESScut
        and the SPOC fallback failed (e.g. a TGLC-only sector with no SPOC
        TPF, or MAST timing out). The TPF panel is a diagnostic, so its
        absence must not abort the whole quicklook.
        """
        if self.pipeline == "gsfc-eleanor-lite":
            err_msg = "Pipeline to be added soon."
            logger.info(err_msg)
            raise NotImplementedError(err_msg)
        if self.pipeline not in FULL_FRAME_TESS_PIPELINES:
            self.tpf = self.get_tpf(sector=self._tpf_sector, author=self.pipeline)
            self.sap_mask = "pipeline"
            return True
        if self.verbose:
            logger.info("Getting TPF with tesscut...")
        try:
            self.tpf = self.get_tpf_tesscut()
            self.sap_mask = "square"
            return True
        except Exception as exc:
            logger.warning(
                f"TESScut TPF search/download failed ({exc}); falling back to SPOC TPF search..."
            )
        try:
            self.tpf = self.get_tpf(sector=self.sector, author="SPOC")
            self.sap_mask = "pipeline"
            return True
        except Exception as exc:
            logger.warning(f"No SPOC TPF either ({exc}); skipping the neighbor-check panels.")
            return False

    def get_tpf_tesscut(self):
        """Download a 15x15 TESS postage stamp.

        Returns
        -------
        tpf : lightkurve.targetpixelfile.TargetPixelFile
            The downloaded TESS postage stamp.
        """
        if self.sector is None:
            raise InvalidInputError("Provide sector for TESScut download.")

        def _download():
            return lk.search_tesscut(self.query_name, sector=self._tpf_sector).download(
                cutout_size=(15, 15), quality_bitmask=self.quality_bitmask
            )

        tpf = self._download_with_retry(_download)
        if tpf is None:
            raise NoDataError("No results from Tesscut search.")
        # remove zeros
        zero_mask = (tpf.flux_err == 0).all(axis=(1, 2))
        if zero_mask.sum() > 0:
            tpf = tpf[~zero_mask]
        return tpf

    def get_planet_params(self):
        """
        e.g.
        https://exofop.ipac.caltech.edu/tess/target.php?id=TOI-6715&json
        """
        if self.verbose:
            logger.info(f"Querying ephemeris for {self.target_name}:")
        try:
            # Use TIC latest uploaded ephem as default
            planet_params = get_params_from_exofop(self.exofop_data, "planet_parameters")
        except (KeyError, TypeError, ValueError, IndexError) as e:
            logger.warning(f"Latest planet params unavailable ({e}); falling back to idx=1")
            planet_params = None
        if planet_params is None:
            planet_params = get_params_from_exofop(self.exofop_data, "planet_parameters", idx=1)
        return planet_params

    def get_toi_radius(self) -> tuple:
        planet_params = self.get_planet_params()
        if planet_params is not None:
            r = planet_params.get("rad")
            r = float(r) if r and (r != "") else np.nan
            re = planet_params.get("rad_e")
            re = float(re) if re and (re != "") else np.nan
            if not math.isnan(r) or not math.isnan(re):
                return (r, re)

    def get_toi_ephem(self, params=["epoch", "per", "dur"]) -> list:
        """
        Query TOI ephemeris from TFOP.

        Parameters
        ----------
        params : list
            List of parameter names to query. Default is ["epoch", "per", "dur"].

        Returns
        -------
        list : list
            A list of tuples, each containing the value and error for the
            corresponding parameter.
        """
        planet_params = self.get_planet_params()
        if planet_params is None:
            return (None, None, None, None)
        if self.verbose:
            logger.info(f"Parameters for {planet_params.get('name', 'unknown')}:")

        # Initialize variables
        toi_epoch = None
        toi_period = None
        toi_dur = None
        toi_depth = None

        # Query values and errors
        for p in params:
            unit = "hr" if p == "dur" else "d"
            val = planet_params.get(p)
            val = float(val) if val else 0.1
            err = planet_params.get(p + "_e")
            err = float(err) if err else 0.1
            if self.verbose:
                logger.info(f"{p}: {val}, {err} {unit}")
            if p == "epoch":
                toi_epoch = np.array((val, err))
                toi_epoch[0] -= TESS_TIME_OFFSET
            elif p == "per":
                toi_period = np.array((val, err))
            elif p == "dur":
                toi_dur = np.array((val, err)) / 24

        # Query depth
        d = planet_params.get("dep_p")
        de = planet_params.get("dep_p_e")
        d = float(d) if d and (d != "") else np.nan
        de = float(de) if de and (de != "") else np.nan
        if not math.isnan(d) or not math.isnan(de):
            toi_depth = np.array((d, de)) / 1e3

        return (toi_epoch, toi_period, toi_dur, toi_depth)

    def run_tls(self):
        """
        Run Transit Least Squares (TLS) on the flattened light curve.

        TLS is a method for searching for transiting exoplanets. It fits a
        transit model to the light curve and computes the power of the
        transit signal at each period. The periodogram is then calculated
        by taking the power values at each period and normalizing them by
        the maximum power.

        Returns
        -------
        tls_results : dict
            A dictionary containing the results of the TLS calculation.
            The keys are the periods and the values are the corresponding
            powers.
        """
        # CDIPS light curve has no flux err
        flux_err = self._tls_flux_err(self.flat_lc)
        power_kwargs = self._tls_power_kwargs()
        gpu_tls = _get_gpu_tls()
        engine_args = (
            self.flat_lc.time.value,
            self.flat_lc.flux.value,
            flux_err,
        )
        ran_on_gpu = False
        if gpu_tls:
            logger.info("Running GTLS on GPU")
            try:
                model = gpu_tls(*engine_args, verbose=self.verbose)
                result = model.power(**power_kwargs)
                self.tls_results = _adapt_gtls_result(result, model)
                ran_on_gpu = True
            except Exception as exc:
                logger.warning(f"GTLS failed ({exc}); retrying with CPU TLS")

        if not ran_on_gpu:
            logger.info("Running TLS on CPU")
            self.tls_results = self._tls_search(self.flat_lc)

        # Advanced vetting flags (also covers the GPU/GTLS path above, which
        # does not compute them itself).
        self.compute_advanced_vetting_metrics()

    def _tls_flux_err(self, lc):
        """Flux error array for TLS.

        CDIPS light curves carry NaN flux errors, so fall back to the LC
        scatter for those cadences instead of failing the search.
        """
        if math.isnan(np.median(lc.flux_err.value)):
            flux_err = np.zeros_like(lc.flux_err)
            flux_err += np.nanstd(lc.flux.value)
        else:
            flux_err = lc.flux_err.value
        return flux_err

    def _tls_power_kwargs(self, period_min=None, period_max=None):
        """Build the shared ``power()`` kwargs for a TLS run."""
        power_kwargs = {
            "period_min": self.Porb_min
            if period_min is None
            else period_min,  # Roche limit default
            "period_max": self.Porb_max if period_max is None else period_max,
        }
        if self.tls_use_threads is not None:
            power_kwargs["use_threads"] = self.tls_use_threads
        if self.use_star_priors:
            logger.info(
                "use_priors=True: fetching ExoFOP stellar parameters (R_star, M_star) for TLS"
            )
            power_kwargs.update(self._stellar_prior_kwargs())
        else:
            logger.info("use_priors=False: TLS will use Sun-like defaults (R_star=1, M_star=1)")
        return power_kwargs

    def _tls_search(self, lc, period_min=None, period_max=None):
        """Run the TLS algorithm on ``lc`` and return its results.

        Used both for the primary detection (on the full flattened light
        curve) and for every iterative re-search (on the transit-masked
        remainder), so threads/priors settings are applied consistently.
        """
        time = lc.time.value
        flux, flux_err = lc.flux.value, self._tls_flux_err(lc)
        span = float(np.ptp(time))
        kwargs = self._tls_power_kwargs(period_min, period_max)
        # Beyond half the baseline TLS only searches if told a single transit
        # may define the grid; winners with <2 transits are discarded below.
        extended = kwargs["period_max"] > span / 2
        if extended:
            kwargs["n_transits_min"] = 1
        try:
            results = tls(time, flux, flux_err, verbose=self.verbose).power(**kwargs)
        except ValueError as e:
            # TLS's duration cache can come out empty for sparse data over the
            # extended range (Kepler-1928: 399 cadences over 26 d); the
            # standard range still works.
            if not extended:
                raise
            logger.warning(
                f"Extended-range TLS failed ({e}); falling back to the standard period range."
            )
            kwargs["period_max"] = span / 2
            kwargs.pop("n_transits_min")
            return tls(time, flux, flux_err, verbose=self.verbose).power(**kwargs)
        if extended and needs_standard_fallback(
            getattr(results, "distinct_transit_count", 2), kwargs["period_max"], span
        ):
            if self.verbose:
                logger.info(
                    f"Extended-range TLS winner (P={float(results.period):.2f} d) has "
                    "<2 transits; falling back to the standard period range."
                )
            kwargs["period_max"] = span / 2
            kwargs.pop("n_transits_min")
            results = tls(time, flux, flux_err, verbose=self.verbose).power(**kwargs)
        return results

    def run_epoch_check(self):
        """Stack the flattened flux at the catalogue ephemeris (if one is known).

        Adds ``epoch_snr``, ``epoch_n_covered``, ``epoch_control_snr_p95`` and
        ``epoch_detected`` to ``tls_results`` (and hence the saved HDF5).
        """
        epoch = getattr(self, "toi_epoch", None)
        period = getattr(self, "toi_period", None)
        dur = getattr(self, "toi_dur", None)
        if epoch is None or period is None or dur is None or self.tls_results is None:
            return None
        res = epoch_check(
            self.flat_lc.time.value,
            self.flat_lc.flux.value,
            float(period[0]),
            float(epoch[0]),
            float(dur[0]),
        )
        for key, val in res.items():
            try:
                self.tls_results[key] = val
            except Exception:
                pass
        if self.verbose:
            if res["epoch_n_covered"] == 0:
                logger.info("Epoch check: no predicted transit is covered by the data.")
            else:
                logger.info(
                    f"Epoch check ({self.ephem_source} ephem): {res['epoch_n_covered']} covered "
                    f"transit(s), stacked dip SNR={res['epoch_snr']:.1f} "
                    f"(random-phase 95th pct {res['epoch_control_snr_p95']:.1f}) -> "
                    f"{'DETECTED' if res['epoch_detected'] else 'not detected'}"
                )
        return res

    def compute_advanced_vetting_metrics(self):
        """Compute advanced vetting metrics for young variable star candidates."""
        if not hasattr(self, "tls_results") or self.tls_results is None:
            return
        for key, val in self._compute_vetting_metrics(
            self.tls_results, getattr(self, "fold_lc", None)
        ).items():
            try:
                self.tls_results[key] = val
            except Exception:
                pass

    def _compute_vetting_metrics(self, results, fold_lc):
        """Shared vetting-metric computation for primary and iterative candidates.

        Returns a dict with ``depth_variance_ratio``, ``duration_consistency_ratio``,
        ``secondary_depth`` and ``secondary_sde``. Metrics that cannot be computed
        from the available data are set to NaN rather than raising.
        """
        metrics = {}

        # 1. Depth variance ratio
        try:
            int_depths = getattr(results, "transit_depths", None)
            if int_depths is not None and len(int_depths) > 1 and np.nanmean(int_depths) > 0:
                metrics["depth_variance_ratio"] = float(
                    np.nanstd(int_depths) / np.nanmean(int_depths)
                )
            else:
                metrics["depth_variance_ratio"] = np.nan
        except Exception:
            metrics["depth_variance_ratio"] = np.nan

        # 2. Duration consistency ratio
        try:
            obs_dur = getattr(results, "duration", None)
            exp_dur = getattr(results, "duration_expected", None)
            if obs_dur and exp_dur and exp_dur > 0:
                metrics["duration_consistency_ratio"] = float(obs_dur / exp_dur)
            else:
                metrics["duration_consistency_ratio"] = np.nan
        except Exception:
            metrics["duration_consistency_ratio"] = np.nan

        # 3. Secondary eclipse search in phase range [0.1, 0.9]
        try:
            if fold_lc is not None:
                phase = fold_lc.phase.value
                flux = fold_lc.flux.value
                sec_mask = (phase >= 0.1) & (phase <= 0.9)
                if np.sum(sec_mask) > 10:
                    sec_std = np.nanstd(flux[sec_mask])
                    sec_min = np.nanmin(flux[sec_mask])
                    sec_depth = 1.0 - sec_min
                    metrics["secondary_depth"] = float(sec_depth)
                    metrics["secondary_sde"] = float(sec_depth / sec_std if sec_std > 0 else np.nan)
                else:
                    metrics["secondary_depth"] = np.nan
                    metrics["secondary_sde"] = np.nan
            else:
                metrics["secondary_depth"] = np.nan
                metrics["secondary_sde"] = np.nan
        except Exception:
            metrics["secondary_depth"] = np.nan
            metrics["secondary_sde"] = np.nan

        return metrics

    def run_iterative_tls(self, min_sde: float = None, max_planets: int = None):
        """Iteratively search for additional planets by masking detected transits.

        The primary TLS detection (``self.tls_results``) is candidate 1. If its
        SDE is below ``min_sde`` no additional search is performed. Otherwise the
        primary transits are masked and TLS is re-run on the remaining flattened
        light curve, repeating until (a) no candidate reaches ``min_sde``, (b) too
        few points remain, (c) the detected period duplicates an earlier
        candidate, or (d) ``max_planets`` candidates have been found.

        Each candidate (with its vetting metrics) is stored in
        ``self.iterative_tls_results``; the primary is always element 0.

        Parameters
        ----------
        min_sde : float, optional
            Minimum SDE for the primary and for each additional candidate.
            Defaults to ``self.min_sde_iterative`` (5.0).
        max_planets : int, optional
            Maximum number of planets to search for. Defaults to
            ``self.max_planets`` (None = search until SDE < min_sde).

        Returns
        -------
        list
            All candidates found (primary first).
        """
        if min_sde is None:
            min_sde = self.min_sde_iterative
        if max_planets is None:
            max_planets = self.max_planets
        if not hasattr(self, "tls_results") or self.tls_results is None:
            self.iterative_tls_results = []
            return self.iterative_tls_results

        candidates = [self.tls_results]
        # Ensure the primary carries the advanced vetting metrics (the GPU/GTLS
        # path returns early from run_tls without computing them).
        if getattr(self.tls_results, "depth_variance_ratio", None) is None:
            self.compute_advanced_vetting_metrics()

        primary_sde = getattr(self.tls_results, "SDE", 0)
        if primary_sde < min_sde:
            if self.verbose:
                logger.info(
                    f"Primary SDE={primary_sde:.2f} < {min_sde}; skipping additional searches."
                )
            self.iterative_tls_results = candidates
            return candidates

        remaining = self.flat_lc[~self.primary_transit_mask(self.flat_lc)]

        while True:
            if len(remaining) < 50:
                if self.verbose:
                    logger.info("Too little data remains after masking; stopping iterative search.")
                break
            if max_planets is not None and len(candidates) >= max_planets:
                if self.verbose:
                    logger.info(f"Reached max_planets={max_planets}; stopping iterative search.")
                break

            results = self._tls_search(remaining)
            sde = getattr(results, "SDE", 0)
            if sde < min_sde:
                if self.verbose:
                    logger.info(f"SDE={sde:.2f} < {min_sde}; stopping iterative search.")
                break
            if self._period_is_duplicate(getattr(results, "period", np.nan), candidates):
                if self.verbose:
                    logger.warning(
                        f"Period P={results.period:.4f} d duplicates an earlier candidate; "
                        "stopping iterative search."
                    )
                break

            fold_lc = remaining.fold(
                period=results.period,
                epoch_time=results.T0,
                normalize_phase=False,
                wrap_phase=results.period / 2,
            )
            for key, val in self._compute_vetting_metrics(results, fold_lc).items():
                try:
                    results[key] = val
                except Exception:
                    pass
            candidates.append(results)
            if self.verbose:
                logger.info(
                    f"Found candidate {len(candidates)}: P={results.period:.4f} d, "
                    f"T0={results.T0:.3f} BTJD, SDE={sde:.2f}, "
                    f"depth={_get_frac_depth(results) * 1e3:.2f} ppt"
                )

            remaining = remaining[~self._transit_mask_for_lc(remaining, results)]

        self.iterative_tls_results = candidates
        return candidates

    @staticmethod
    def _period_is_duplicate(period, candidates, rtol=1e-3, atol=1e-4):
        """True if ``period`` (days) matches any already-found candidate period."""
        for cand in candidates:
            if np.isclose(period, getattr(cand, "period", np.inf), rtol=rtol, atol=atol):
                return True
        return False

    @staticmethod
    def _transit_mask_for_lc(lc, results, margin=1.2):
        """Boolean transit mask covering the signal in ``results`` with extra margin."""
        return lc.create_transit_mask(
            transit_time=getattr(results, "T0", 0),
            period=getattr(results, "period", 1),
            duration=getattr(results, "duration", 0.1) * margin,
        )

    def _stellar_prior_kwargs(self):
        """Pull R_star, M_star (and ±1σ bounds) from ExoFOP for the TLS prior.

        Returns a dict of TLS ``.power()`` kwargs. Any value that is missing
        or non-finite is omitted, letting TLS fall back to its Sun-like
        default for that parameter alone.
        """
        try:
            star_params = get_params_from_exofop(self.exofop_data, name="stellar_parameters", idx=1)
        except (KeyError, TypeError, ValueError, IndexError):
            try:
                star_params = get_params_from_exofop(self.exofop_data, name="stellar_parameters")
            except (KeyError, TypeError, ValueError, IndexError) as e:
                logger.warning(f"No ExoFOP stellar params available ({e}); skipping priors")
                return {}

        def _as_float(key):
            try:
                v = float(star_params.get(key))
            except (TypeError, ValueError):
                return None
            return v if np.isfinite(v) else None

        kwargs = {}
        # R_star and ±1σ bounds (R_star ± srad_e), clipped to TLS's allowed range.
        r = _as_float("srad")
        r_err = _as_float("srad_e")
        if r is not None:
            kwargs["R_star"] = r
            if r_err is not None and r_err > 0:
                kwargs["R_star_min"] = max(r - r_err, 0.13)
                kwargs["R_star_max"] = r + r_err
        # M_star and ±1σ bounds (M_star ± mass_e).
        m = _as_float("mass")
        m_err = _as_float("mass_e")
        if m is not None:
            kwargs["M_star"] = m
            if m_err is not None and m_err > 0:
                kwargs["M_star_min"] = max(m - m_err, 0.1)
                kwargs["M_star_max"] = m + m_err

        if not kwargs:
            logger.warning("ExoFOP stellar params present but unusable; skipping priors")
        else:
            # Logged unconditionally (not gated on verbose) so users can always
            # verify which prior values reached TLS when they enabled --use_priors.
            # Round to 2 significant figures for readability (display only; the
            # full-precision kwargs are still returned to TLS).
            pretty = {k: float(f"{v:.2g}") for k, v in kwargs.items()}
            logger.info(f"Using ExoFOP stellar priors for TLS: {pretty}")
        return kwargs

    def init_gls(self):
        masked_lc = self.raw_lc[~self.tmask]
        if self.pipeline == "pathos":
            data = np.array(
                [
                    np.asarray(masked_lc.time.value, dtype=float),
                    np.asarray(masked_lc.flux.value, dtype=float),
                ]
            )
        else:
            data = np.array(
                [
                    np.asarray(masked_lc.time.value, dtype=float),
                    np.asarray(masked_lc.flux.value, dtype=float),
                    np.asarray(masked_lc.flux_err.value, dtype=float),
                ]
            )
        if self.verbose:
            msg = "Estimating rotation period using Generalized Lomb-Scargle (GLS) periodogram..."
            logger.info(msg)
        return Gls(data, Pbeg=self.Porb_min, Pend=self.Porb_max, verbose=self.verbose)

    def flatten_raw_lc(self):
        """
        Flatten the raw light curve using wotan, or Notch/LOCoR for
        ``flatten_method`` in ("notch", "locor").

        Returns
        -------
        flat_lc : lightkurve.lightcurve.TessLightCurve
            The flattened light curve.
        trend_lc : lightkurve.lightcurve.TessLightCurve
            The trend light curve.

        """
        if self.verbose:
            logger.info(f"Using {self.flatten_method} method to flatten raw lc.")
        if float(self.window_length) == 0.0 and self.flatten_method in NOTCH_LOCOR_METHODS:
            # Notch/LOCoR aren't wotan methods, so the injection-recovery grid
            # search below (which calls wotan.flatten internally) can't tune
            # their window/period. Fall back to a sensible default instead.
            if self.flatten_method == "locor":
                try:
                    gls_obj = self.init_gls()
                    self.window_length = float(gls_obj.P)
                except Exception as exc:
                    self.window_length = 1.0
                    if self.verbose:
                        logger.warning(
                            f"Could not compute rotation period via GLS for LOCoR ({exc}); "
                            f"using window_length={self.window_length:.2f} d."
                        )
            else:
                self.window_length = 0.5
            self.window_length_opt = None
            if self.verbose:
                logger.info(
                    f"Auto window-length search isn't available for {self.flatten_method}; "
                    f"using window_length={self.window_length:.2f} d."
                )
        elif float(self.window_length) == 0.0:
            np.random.default_rng(1)

            # Window grid must safely cover the injected duration range below
            # at Hippke & Heller 2019's ~2.2-3x duration ratio (below ~2.2x a
            # wotan window starts distorting the transit it's supposed to
            # preserve). The grid used to be a fixed 0.1-0.6 d regardless of
            # these durations; with duration_max=0.5 d that excluded every
            # window >=3x the longest duration actually being tested, so the
            # grid search could never select a safely wide window for it.
            duration_min, duration_max = 0.05, 0.5
            window_low = 2.2 * duration_min
            window_high = 3.0 * duration_max
            logger.info(
                f"Optimizing window length between {window_low:.2f} and "
                f"{window_high:.2f} d. This may take several minutes."
            )
            window_lengths = np.linspace(window_low, window_high, 10)

            time = self.raw_lc.time.value
            flux = self.raw_lc.flux.value

            baseline = time[-1] - time[0]
            half_baseline = baseline / 2
            Nmodels = 5
            t0s = np.random.uniform(low=0, high=half_baseline, size=Nmodels)
            periods = np.random.uniform(low=0.1, high=half_baseline, size=Nmodels)
            durations = np.random.uniform(low=duration_min, high=duration_max, size=Nmodels)
            depths = np.random.uniform(low=0.005, high=0.05, size=Nmodels)

            injections = []
            for t0, p, dur, depth in zip(t0s, periods, durations, depths):
                injections.append(
                    InjectionParams(t0=time[0] + t0, period=p, duration=dur, depth=depth)
                )

            # Run grid
            df = run_grid(time, flux, window_lengths, injections, method=self.flatten_method)

            # Select optimal window
            self.window_length, valid = select_best_window(df)
            logger.info(
                f"Best window length: {self.window_length:.1f} d using {self.flatten_method}."
            )
            self.window_length_opt = self.window_length
        else:
            self.window_length_opt = None

        # Rotation-aware window and spline tuning if window_length is default or unset.
        # Skipped for LOCoR: its window_length is a rotation *period*, not a
        # sub-day detrending window, so the 0.5/1.0 d heuristic below doesn't apply.
        if (
            getattr(self, "_use_rotation_aware_window", False)
            and getattr(self, "toi_dur", None) is None
            and self.flatten_method != "locor"
        ):
            try:
                gls_obj = self.init_gls()
                prot = float(gls_obj.P)
                if prot <= 2.0:
                    self.window_length = 0.5
                else:
                    self.window_length = 1.0
                if self.verbose:
                    logger.info(
                        f"Rotation-aware tuning: P_rot={prot:.2f} d -> window_length={self.window_length:.2f} d"
                    )
            except Exception as exc:
                if self.verbose:
                    logger.warning(
                        f"Could not compute rotation period via GLS for dynamic window tuning ({exc})"
                    )

        if self.flatten_method in NOTCH_LOCOR_METHODS:
            # Notch and LOCoR (Rizzuto et al. 2017; arizzuto/Notch_and_LOCoR),
            # reimplemented in quicklook.notch_locor -- not part of wotan.
            if self.flatten_method == "notch":
                wflat_lc, wtrend_lc = notch_flatten(
                    self.raw_lc.time.value,
                    self.raw_lc.flux.value,
                    window_length=self.window_length,
                )
            else:
                wflat_lc, wtrend_lc = locor_flatten(
                    self.raw_lc.time.value,
                    self.raw_lc.flux.value,
                    period=self.window_length,
                )
        else:
            # Additional kwargs for wotan flatten
            flatten_kwargs = {
                "method": self.flatten_method,
                "robust": True,
                "kernel": self.gp_kernel,
                "kernel_size": self.gp_kernel_size,
                "window_length": self.window_length,
                "edge_cutoff": self.edge_cutoff,
                "break_tolerance": 1,
                "return_trend": True,
                "cval": 5.0,
            }

            # Dynamic max_splines calculation for spline methods based on P_rot
            if self.flatten_method in ("rspline", "pspline", "spline"):
                try:
                    gls_obj = self.init_gls()
                    prot = float(gls_obj.P)
                    if prot <= 2.0:
                        max_splines = 200
                    elif prot <= 10.0:
                        max_splines = int(200 / prot)
                    else:
                        max_splines = 25
                    flatten_kwargs["max_splines"] = max_splines
                    if self.verbose:
                        logger.info(
                            f"Rotation-aware spline tuning: P_rot={prot:.2f} d -> max_splines={max_splines}"
                        )
                except Exception as exc:
                    if self.verbose:
                        logger.warning(f"Could not set max_splines dynamically ({exc})")

            # https://github.com/hippke/wotan#available-detrending-algorithms
            wflat_lc, wtrend_lc = flatten(
                # Array of time values
                self.raw_lc.time.value,
                # Array of flux values
                self.raw_lc.flux.value,
                **flatten_kwargs,
            )

        if self.sigma_clip_flat is not None:
            msg = "Applying sigma clip on flattened lc with "
            msg += f"(lower,upper)=({self.sigma_clip_flat})"
            if self.verbose:
                logger.info(msg)
            idx = sigma_clip(
                wflat_lc,
                sigma_lower=self.sigma_clip_flat[0],
                sigma_upper=self.sigma_clip_flat[1],
            ).mask
        else:
            idx = np.zeros_like(wflat_lc, dtype=bool)

        valid_mask = ~idx
        flat_lc, trend_lc = self.raw_lc.flatten(return_trend=True)
        flat_lc = flat_lc[valid_mask]
        trend_lc = trend_lc[valid_mask]
        trend_lc.flux = wtrend_lc[valid_mask]
        flat_lc.flux = wflat_lc[valid_mask]
        return flat_lc, trend_lc

    def get_transit_mask(self):
        """
        Generate a mask for the transit based on the user-provided
        transit ephemeris or the ephemeris from the TOI portal.

        Returns
        -------
        tmask : np.ndarray
            A boolean mask where the transit is True and the out-of-transit
            periods are False.
        """
        if np.all([self.toi_epoch, self.toi_period, self.toi_dur]):
            # Use the user-provided transit ephemeris to create the mask
            tmask = self.raw_lc.create_transit_mask(
                transit_time=self.toi_epoch[0],
                period=self.toi_period[0],
                duration=self.toi_dur[0],
            )
        else:
            # If no transit ephemeris is provided, create an empty mask
            tmask = np.zeros_like(self.raw_lc.time.value, dtype=bool)
        return tmask

    def _warn_no_transits_in_sector(self):
        """Warn that the known ephemeris predicts no transit in the data.

        For long-period planets the transit often falls outside the ~27-day
        TESS sector that was downloaded, so there is simply nothing to mask.
        This is expected, not an error: the quicklook (flattening, TLS search,
        plots) still runs. Report the sector's time coverage and the nearest
        predicted transit so the user can pick a sector that actually contains
        a transit.
        """
        time = self.raw_lc.time.value
        t_start, t_end = float(np.nanmin(time)), float(np.nanmax(time))
        epoch, period = float(self.toi_epoch[0]), float(self.toi_period[0])
        msg = (
            f"No transit predicted by the {self.ephem_source} ephem falls "
            f"within sector {self._sector_display()} (BTJD {t_start:.2f}-{t_end:.2f}); "
            "nothing to mask."
        )
        if period > 0:
            # Nearest transit epoch to the middle of the observed baseline.
            t_mid = 0.5 * (t_start + t_end)
            nearest = epoch + round((t_mid - epoch) / period) * period
            gap = max(t_start - nearest, nearest - t_end, 0.0)
            msg += (
                f" Nearest predicted transit at BTJD {nearest:.2f} "
                f"({gap:.1f} d outside coverage, P={period:.4f} d)."
            )
        if self.all_sectors is not None and len(self.all_sectors) > 1:
            msg += f" Try another sector: {self.all_sectors}."
        logger.warning(msg)

    def _summary_sections(self):
        """Build the summary as structured (label, value) rows per section.

        Single source of truth for both the plain-text summary
        (:meth:`make_summary_info`) and the rendered panel
        (:meth:`_render_summary_panel`). Units use mathtext so they render
        cleanly in matplotlib. TOI/TFOP reference values are emitted as their
        own rows (rather than appended to the TLS value) so each row stays
        short and the label/value columns line up.

        Returns
        -------
        list of (str, list of (str, str))
            ``[(section_title, [(label, value), ...]), ...]``
        """
        try:
            # Use the TIC stellar parameters as default
            star_params = get_params_from_exofop(self.exofop_data, name="stellar_parameters", idx=1)
        except (KeyError, TypeError, ValueError, IndexError) as e:
            logger.warning(f"Stellar params at idx=1 unavailable ({e}); using latest")
            star_params = get_params_from_exofop(self.exofop_data, name="stellar_parameters")
        params = {}
        param_names = ["srad", "mass", "teff", "logg", "dist"]
        for name in param_names:
            try:
                # Convert to float or int as needed
                params[name] = (
                    int(float(star_params.get(name)))
                    if name == "teff"
                    else float(star_params.get(name))
                )
            except (TypeError, ValueError) as e:
                logger.debug(f"Could not parse stellar param '{name}': {e}")
                params[name] = np.nan
            try:
                # Convert to float or int as needed
                params[name + "_e"] = (
                    int(float(star_params.get(name + "_e")))
                    if name == "teff"
                    else float(star_params.get(name + "_e"))
                )
            except (TypeError, ValueError) as e:
                logger.debug(f"Could not parse stellar param '{name}_e': {e}")
                params[name + "_e"] = np.nan
        # Get the meta data
        meta = self.raw_lc.meta
        # Calculate the planet radius
        Rp = self.tls_results["rp_rs"] * params["srad"] * u.Rsun.to(u.Rearth)
        # If the pipeline is Spoc or Tess-Spoc, correct for dilution. Keep the
        # applied contamination factor so it can be shown in the summary panel.
        if self.pipeline in ["spoc", "tess-spoc"]:
            contam_factor = np.sqrt(1 + meta["CROWDSAP"])
            Rp_true = Rp * contam_factor
        else:
            # Otherwise, use the raw radius
            contam_factor = None
            Rp_true = Rp

        pm = r"$\pm$"
        r_earth = r"R$_{\oplus}$"
        r_sun = r"R$_{\odot}$"
        m_sun = r"M$_{\odot}$"
        tls = self.tls_results

        # --- Candidate properties ---
        candidate = [
            ("SDE", f"{tls.SDE:.2f}"),
            ("Period (TLS)", f"{tls.period:.4f}{pm}{tls.period_uncertainty:.4f} d"),
        ]
        if self.toi_period is not None:
            candidate.append(
                (
                    f"Period ({self.ephem_source})",
                    f"{self.toi_period[0]:.4f}{pm}{self.toi_period[1]:.4f} d",
                )
            )
        # BTJD = BJD - TESS_TIME_OFFSET (2457000); compact unit keeps the value
        # short enough for the half-width column.
        candidate.append(("T0 (TLS)", f"{tls.T0:.4f} BTJD"))
        if self.toi_period is not None:
            candidate.append(
                (
                    f"T0 ({self.ephem_source})",
                    f"{self.toi_epoch[0]:.4f}{pm}{self.toi_epoch[1]:.4f}",
                )
            )
        candidate.append(("Duration (TLS)", f"{tls.duration * 24:.2f} hr"))
        if self.toi_dur is not None:
            candidate.append(
                (
                    f"Duration ({self.ephem_source})",
                    f"{self.toi_dur[0] * 24:.2f}{pm}{self.toi_dur[1] * 24:.2f} hr",
                )
            )
        candidate.append(("Depth (TLS)", f"{_get_frac_depth(tls) * 1e3:.2f} ppt"))
        if self.toi_depth is not None:
            candidate.append(
                ("Depth (TFOP)", f"{self.toi_depth[0]:.1f}{pm}{self.toi_depth[1]:.1f} ppt")
            )
        # Contamination factor actually applied to undilute Rp (sqrt(1+CROWDSAP)
        # for SPOC-family light curves); shown so the correction is auditable.
        if contam_factor is not None:
            candidate.append(("Contamination", f"{contam_factor:.2g}"))
        is_undiluted = (meta["FLUX_ORIGIN"].lower() in ("pdcsap", "sap")) or (
            self.pipeline == "tglc"
        )
        if is_undiluted:
            candidate.append(("Rp (TLS)", f"{Rp_true:.2f} {r_earth}"))
        else:
            candidate.append(("Rp (diluted)", f"{Rp:.2f} {r_earth}"))
            candidate.append(("Rp (undiluted)", f"{Rp_true:.2f} {r_earth}"))
        if self.toi_rp is not None:
            candidate.append(
                ("Rp (TFOP)", f"{self.toi_rp[0]:.2f}{pm}{self.toi_rp[1]:.2f} {r_earth}")
            )
        candidate.append(("Odd-Even", f"{tls.odd_even_mismatch:.2f} " + r"$\sigma$"))
        # Number of transits sampled by the data. Present on CPU TLS results;
        # the GPU adapter does not expose it, so guard with getattr.
        n_transits = getattr(tls, "distinct_transit_count", None)
        # NaN when TLS fit no transit (HD 114082).
        if n_transits is not None and np.isfinite(n_transits):
            candidate.append(("Num. Transits", f"{int(n_transits)}"))
        candidate.append(("Available sectors", self._format_available_sectors(self.all_sectors)))
        n_planets = len(getattr(self, "iterative_tls_results", None) or [])
        if n_planets > 1:
            candidate.append(("Num. Planets", f"{n_planets}"))

        if self.event_class in ("mono", "duo"):
            candidate = self._single_event_candidate_rows(params["srad"])
        elif self.has_ttv_layout:
            candidate.extend(self._ttv_candidate_rows())

        # --- Stellar properties ---
        per = 2 * np.pi * params["srad"] * u.Rsun.to(u.km)
        t = self.Prot_ls * u.day.to(u.second)
        coords = self.target_coord.to_string("decimal").split()
        mags = self.exofop_data["magnitudes"][0]
        stellar = [
            ("Rstar", f"{params['srad']:.2f}{pm}{params['srad_e']:.2f} " + r_sun),
            ("Mstar", f"{params['mass']:.2f}{pm}{params['mass_e']:.2f} " + m_sun),
            ("Teff", f"{params.get('teff')}{pm}{params.get('teff_e')} K"),
            ("logg", f"{params['logg']:.2f}{pm}{params['logg_e']:.2f} cgs"),
            ("Distance", f"{params['dist']:.1f}{pm}{params['dist_e']:.1f} pc"),
            ("Prot", f"{self.Prot_ls:.2f} d"),
            ("Vsini", f"{per / t:.2f} km/s"),
            ("RA, Dec", f"{float(coords[0])}, {float(coords[1])}"),
            (f"{mags['band']}mag", f"{float(mags['value']):.1f}"),
            ("Gaia DR2 ID", f"{self.gaiaid}"),
        ]
        # RUWE (Gaia DR3 astrometric fit quality; >~1.4 hints at an unresolved
        # companion/blend). Shown when the Gaia catalog query returned a value.
        ruwe = getattr(self, "gaia_ruwe", np.nan)
        if np.isfinite(ruwe):
            stellar.append(("RUWE", f"{ruwe:.2f}"))
        # Unset when the neighbor-check panels were skipped (no TPF).
        nearby_sep = getattr(self, "nearby_star_sep", None)
        if nearby_sep is not None:
            val = nearby_sep.to(u.arcsec) if nearby_sep < 1 * u.arcmin else nearby_sep
            stellar.append(("Nearby sep", f"{val:.1f}"))
        if self.simbad_obj_type is not None:
            stellar.append(("SIMBAD type", f"{self.simbad_obj_type}"))

        return [
            ("Candidate Properties", candidate),
            ("Stellar Properties", stellar),
        ]

    @staticmethod
    def _format_available_sectors(all_sectors, max_listed=8, sectors_per_line=4):
        """Format available sectors for the summary panel.

        Avoid long inline sector lists because they can run into the
        neighbouring summary column. Beyond max_listed sectors the list is
        replaced by a count.
        """
        if all_sectors is None:
            return ""
        try:
            sectors = list(all_sectors)
        except TypeError:
            sectors = [all_sectors]
        if len(sectors) == 0:
            return ""
        if len(sectors) > max_listed:
            return f"{len(sectors)} sectors available"

        sector_chunks = [
            sectors[i : i + sectors_per_line] for i in range(0, len(sectors), sectors_per_line)
        ]
        lines = [", ".join(str(s) for s in chunk) for chunk in sector_chunks]
        return "\n".join(lines)

    @staticmethod
    def _extract_ruwe(gaia_sources, gaiaid):
        """Return the target's Gaia RUWE from a Gaia catalog table, or NaN.

        The target row is matched on ``source_id`` (the same convention the
        Gaia overlay plots use). Returns NaN when the table is missing, lacks a
        ``ruwe``/``source_id`` column, has no matching row, or holds a
        non-numeric value — so callers can guard on ``np.isfinite``.
        """
        if gaia_sources is None:
            return np.nan
        columns = getattr(gaia_sources, "columns", [])
        if "ruwe" not in columns or "source_id" not in columns:
            return np.nan
        try:
            match = gaia_sources.loc[
                gaia_sources["source_id"].astype("int64") == int(gaiaid), "ruwe"
            ]
        except (TypeError, ValueError):
            return np.nan
        if len(match) == 0:
            return np.nan
        try:
            return float(match.iloc[0])
        except (TypeError, ValueError):
            return np.nan

    def make_summary_info(self):
        """
        Generate a plain-text summary of the TLS results, stellar params,
        and other useful information.

        Returns
        -------
        msg : str
            A summary string with the results.
        """
        lines = []
        for title, rows in self._summary_sections():
            lines.append(title)
            lines.append("-" * 30)
            for label, value in rows:
                lines.append(f"{label}: {value}")
            lines.append("")
        return "\n".join(lines)

    def _render_summary_panel(self, ax):
        """Render the summary as aligned label/value rows grouped by section.

        Each section is placed in its own column (side by side) so the two
        stacks stay short enough to fit the panel without overflowing the
        bottom edge. Everything is positioned in axes-fraction coordinates
        (robust to figure size), with bold section headers, a thin divider,
        and a monospace value column so units and uncertainties line up. Font
        size follows the figure default so the panel matches the neighbouring
        plots; the row pitch is sized to the longest column so the content
        always fits regardless of how many optional rows appear.
        """
        sections = self._summary_sections()
        ax.axis("off")
        n_cols = len(sections)
        col_w = 1.0 / n_cols
        header_fs = pl.rcParams["font.size"]
        row_fs = header_fs - 1

        # Pitch is set by the longest column plus generous overhead (header,
        # divider and a comfortable bottom margin) so the final row of the
        # tallest section stays clear of the panel's bottom edge instead of
        # being clipped when the figure is compressed (e.g. by the suptitle).
        def row_units(rows):
            return sum(
                max(str(label).count("\n"), str(value).count("\n")) + 1 for label, value in rows
            )

        max_units = max(row_units(rows) for _, rows in sections) + 5.0
        top, bottom = 0.99, 0.03
        dy = (top - bottom) / max_units
        # Nudge each column leftward: the first (Candidate) column into the left
        # margin so its labels line up with the "DEC" ylabel of the archival
        # panel directly above, and the second (Stellar) column slightly left so
        # it isn't pushed hard against the right edge.
        col_pads = (-0.195, -0.10)
        for col, (title, rows) in enumerate(sections):
            x0 = col * col_w + (col_pads[col] if col < len(col_pads) else 0.0)
            label_x = x0 + 0.01
            value_x = x0 + 0.25
            y = top
            ax.text(
                x0,
                y,
                title,
                transform=ax.transAxes,
                fontsize=header_fs,
                fontweight="bold",
                va="top",
            )
            y -= dy * 0.8
            ax.plot(
                [x0, x0 + col_w - 0.03],
                [y, y],
                transform=ax.transAxes,
                color="0.7",
                lw=0.8,
            )
            y -= dy * 0.9
            for label, value in rows:
                label_str = str(label)
                value_str = str(value)
                units = max(label_str.count("\n"), value_str.count("\n")) + 1
                ax.text(
                    label_x,
                    y,
                    label_str,
                    transform=ax.transAxes,
                    fontsize=row_fs,
                    color="0.30",
                    va="top",
                )
                # Keep the first value line beside its label, but wrap any
                # continuation lines (e.g. a wrapped available-sector list)
                # back to the label column so they line up under the label
                # instead of trailing beneath the value column. The leading
                # newline reserves the first (label) row so the wrapped lines
                # fall on the same baselines they would as a single block.
                first_value, _, rest_value = value_str.partition("\n")
                ax.text(
                    value_x,
                    y,
                    first_value,
                    transform=ax.transAxes,
                    fontsize=row_fs,
                    family="monospace",
                    va="top",
                )
                if rest_value:
                    ax.text(
                        label_x,
                        y,
                        "\n" + rest_value,
                        transform=ax.transAxes,
                        fontsize=row_fs,
                        family="monospace",
                        va="top",
                    )
                y -= dy * units

    def append_tls_results(self):
        """
        Append TLS results to the TLS results dictionary.

        This will add the raw and flattened light curves, period limits, and
        TFOP parameters to the TLS results dictionary.

        Returns
        -------
        None
        """
        # Append meta
        self.tls_results["pipeline"] = self.pipeline
        self.tls_results["flux_type"] = self.flux_type
        self.tls_results["exptime"] = self.exptime
        self.tls_results["cadence"] = self.cadence

        # Append the raw light curve
        self.tls_results["time_raw"] = self.raw_lc.time.value
        self.tls_results["flux_raw"] = self.raw_lc.flux.value
        self.tls_results["err_raw"] = self.raw_lc.flux_err.value

        # Append the flattened light curve
        self.tls_results["time_flat"] = self.flat_lc.time.value
        self.tls_results["flux_flat"] = self.flat_lc.flux.value
        self.tls_results["err_flat"] = self.flat_lc.flux_err.value

        # Append the period limits
        self.tls_results["Porb_min"] = self.Porb_min
        self.tls_results["Porb_max"] = self.Porb_max
        self.tls_results["depth_ppt"] = _get_frac_depth(self.tls_results) * 1e3

        # Append the TFOP parameters (also available in meta)
        self.tls_results["period_toi"] = self.toi_period
        self.tls_results["T0_toi"] = self.toi_epoch
        self.tls_results["duration_toi"] = self.toi_dur
        self.tls_results["depth_toi"] = self.toi_depth
        self.toi_rp = self.get_toi_radius()
        self.tls_results["Rp_toi"] = self.toi_rp

        # A detection built on fewer than 2 transits with actual in-transit
        # data can't establish periodicity at all -- TLS can still report a
        # high SDE by pairing one real dip with an empty data gap. Flag that
        # here instead of leaving it to be discovered downstream.
        counts = self.tls_results.get("per_transit_count")
        depths = self.tls_results.get("transit_depths")
        if counts is not None and depths is not None:
            n_valid = _count_valid_transits(counts, depths)
            self.tls_results["n_valid_transits"] = n_valid
            if n_valid < 2:
                logger.warning(
                    f"Only {n_valid} transit(s) have both in-transit data and a "
                    f"measurable depth (SDE={self.tls_results.get('SDE')}); "
                    "periodicity can't be established from this alone."
                )
        else:
            self.tls_results["n_valid_transits"] = None

        # Append the Gaia ID, TIC ID, and TOI ID
        self.tls_results["gaiaid"] = int(self.gaiaid)
        self.tls_results["ticid"] = int(self.ticid)
        self.tls_results["toiid"] = self.toiid
        self.tls_results["sector"] = self.sector

        # Append baseline model used to flatten raw flux
        self.tls_results["flatten_method"] = self.flatten_method
        self.tls_results["window_length"] = self.window_length

        # Iterative multi-planet search metadata: whether it was requested and,
        # if any additional candidates were found, their parameters and vetting
        # metrics (as h5io-serializable summary rows).
        self.tls_results["iterative_search"] = bool(getattr(self, "iterative_search", False))
        n_planets = len(getattr(self, "iterative_tls_results", None) or [])
        self.tls_results["n_planets"] = n_planets if n_planets else 1
        if n_planets > 1:
            self.tls_results["iterative_candidates"] = self._summarize_iterative_candidates()

        # Record whether ExoFOP stellar parameters were used as TLS priors,
        # and (if so) which prior values reached the search. Stored in the
        # HDF5 output so downstream analysis can tell prior-informed runs
        # apart from Sun-like-default runs.
        self.tls_results["use_star_priors"] = bool(self.use_star_priors)
        if self.use_star_priors:
            prior_kwargs = self._stellar_prior_kwargs()
            self.tls_results["star_prior_kwargs"] = prior_kwargs
            # Also expose R_star/M_star as plain top-level fields -- nested
            # only inside star_prior_kwargs, a downstream reader has to know
            # that dict and its TLS-kwarg-style key names (R_star, not
            # Rstar) exist to find them. Without this, computing a planet
            # radius in physical units (R_Earth) downstream was only
            # possible from the dimensionless rp_rs ratio, which is how a
            # 3.6 R_Jup eclipsing binary (rp_rs=0.227 around a 1.62 R_sun
            # star) passed a rp_rs-only size check during vetting work on
            # this pipeline's output.
            self.tls_results["Rstar"] = prior_kwargs.get("R_star")
            self.tls_results["Mstar"] = prior_kwargs.get("M_star")

        # Append exofop data (TICv8)
        self.tls_results["exofop_data"] = self.exofop_data

        # Prot_ls is set by whichever periodogram method actually succeeded
        # (gls, or the astropy lombscargle fallback above, or the plain
        # lombscargle/bls branch when pg_method != "gls") -- save it
        # unconditionally so a rotation-period estimate survives even when
        # self.gls itself is None.
        self.tls_results["Prot_ls"] = getattr(self, "Prot_ls", None)

        # Append the Gls results
        if self.gls is not None:
            self.tls_results["power_gls"] = (
                self.gls.power.max(),
                self.gls.power.std(),
            )
            self.tls_results["Prot_gls"] = (
                self.gls.hpstat["P"],
                self.gls.hpstat["e_P"],
            )
            self.tls_results["amp_gls"] = (
                self.gls.hpstat["amp"],
                self.gls.hpstat["e_amp"],
            )

        # Flag when the TLS period aliases the star's own rotation period --
        # manually re-derived by hand several times during vetting work on
        # this pipeline's output (e.g. a candidate at period=0.806d next to
        # Prot_ls=0.804d, an almost exact 1:1 alias masquerading as a
        # transit). Prefer Prot_gls when the GLS periodogram succeeded,
        # otherwise fall back to Prot_ls.
        prot = None
        if self.tls_results.get("Prot_gls") is not None:
            prot = self.tls_results["Prot_gls"][0]
        elif self.tls_results.get("Prot_ls") is not None:
            prot = self.tls_results["Prot_ls"]
        porb = self.tls_results.get("period")
        if prot is not None and porb is not None and prot > 0:
            rel_tol = 0.03
            self.tls_results["rotation_alias_risk"] = any(
                abs(porb - n * prot) <= rel_tol * n * prot for n in (1, 2, 3, 1 / 2, 1 / 3)
            )
        else:
            self.tls_results["rotation_alias_risk"] = None

        # Append the Simbad object type
        if self.simbad_obj_type is not None:
            self.tls_results["simbad_obj"] = self.simbad_obj_type

    def _summarize_iterative_candidates(self):
        """Serialize each found planet into an h5io-safe summary row dict."""
        rows = []
        for i, res in enumerate(getattr(self, "iterative_tls_results", []), start=1):
            rows.append(
                {
                    "planet": i,
                    "period": getattr(res, "period", np.nan),
                    "period_uncertainty": getattr(res, "period_uncertainty", np.nan),
                    "T0": getattr(res, "T0", np.nan),
                    "duration": getattr(res, "duration", np.nan),
                    "depth": getattr(res, "depth", np.nan),
                    "rp_rs": getattr(res, "rp_rs", np.nan),
                    "SDE": getattr(res, "SDE", np.nan),
                    "FAP": getattr(res, "FAP", np.nan),
                    "odd_even_mismatch": getattr(res, "odd_even_mismatch", np.nan),
                    "per_transit_count": getattr(res, "distinct_transit_count", np.nan),
                    "depth_variance_ratio": getattr(res, "depth_variance_ratio", np.nan),
                    "duration_consistency_ratio": getattr(
                        res, "duration_consistency_ratio", np.nan
                    ),
                    "secondary_depth": getattr(res, "secondary_depth", np.nan),
                    "secondary_sde": getattr(res, "secondary_sde", np.nan),
                }
            )
        return rows

    def _mask_all_except(self, planet_index, margin=1.2):
        """Boolean mask over ``self.flat_lc`` covering every candidate but ``planet_index``."""
        mask = np.zeros(len(self.flat_lc), dtype=bool)
        for j, res in enumerate(getattr(self, "iterative_tls_results", []), start=1):
            if j == planet_index:
                continue
            mask |= self._transit_mask_for_lc(self.flat_lc, res, margin=margin)
        return mask

    def _save_iterative_candidates(self, fp):
        """Save a periodogram + odd-even PNG and a TLS HDF5 for every extra planet.

        Only the candidate's own transits are retained in the fold (all other
        found candidates are masked out) so the panel shows the signal as it
        appears in the raw flattened data.
        """
        for idx, results in enumerate(getattr(self, "iterative_tls_results", [])[1:], start=2):
            cand_png = fp.with_stem(fp.stem + f"_p{idx}").with_suffix(".png")
            cand_h5 = Path(self.outdir, fp.name + f"_tls_p{idx}").with_suffix(".h5")

            fold_lc = self.flat_lc[~self._mask_all_except(idx)].fold(
                period=results.period,
                epoch_time=results.T0,
                normalize_phase=False,
                wrap_phase=results.period / 2,
            )
            fig2, axes2 = pl.subplots(1, 2, figsize=(14, 5), tight_layout=True)
            plot_tls(
                results,
                toi_period=None,
                period_min=self.Porb_min,
                period_max=self.Porb_max,
                ax=axes2[0],
            )
            plot_odd_even_transit(
                fold_lc,
                results,
                bin_mins=10,
                markersize=6,
                ax=axes2[1],
                phase_xlim=self.phase_xlim,
            )
            fig2.suptitle(
                f"Planet {idx}: P={results.period:.4f} d | SDE={results.SDE:.2f} | "
                f"Rp/Rs={results.rp_rs:.3f} | depth={_get_frac_depth(results) * 1e3:.2f} ppt",
                y=1.0,
                fontsize=14,
            )
            if self.savefig:
                fig2.savefig(cand_png, dpi=150, bbox_inches="tight")
                logger.info(f"Saved: {cand_png}")
            if self.savetls:
                h5io.save(cand_h5, results)
                logger.info(f"Saved: {cand_h5}")
            pl.close(fig2)

    def plot_tql(self, return_fig_and_paths=False, **kwargs: dict) -> pl.Figure:
        """
        Plot a TQL report.

        Parameters
        ----------
        **kwargs: dict
            Keyword arguments passed to `TessQuickLook`.

        Returns
        -------
        fig: pl.Figure
            The plotted figure.
        """

        if self.verbose:
            logger.info("Creating panels...")
        fig, axes = pl.subplots(3, 3, figsize=(16, 12), tight_layout=True)

        if self.verbose:
            logger.info("Plotting raw light curve...")
        # Mono/duo layouts break the time axis between separate sectors.
        self._panel0_axes = self._time_panel_axes(axes.flatten()[0])
        label = f"baseline trend\nmethod={self.flatten_method}"
        if self.window_length_opt is None:
            label += "(window_size="
        else:
            label += "(optim. window_size="
        label += f"{self.window_length:.2f})"
        for ax in self._panel0_axes:
            self.raw_lc.scatter(ax=ax, label=f"raw (exp={self.exptime} s)")
            self.trend_lc.plot(ax=ax, color="r", lw=2, label=label)
        finish_broken_axes(self._panel0_axes)

        if self.verbose:
            logger.info("Running Lomb-Scargle periodogram...")

        ref_period = self.toi_period[0] if self.toi_period is not None else None
        ax = axes.flatten()[1]
        if self.pg_method == "gls":
            self.gls = self.init_gls()
            ax = plot_gls_periodogram(
                self.gls,
                offset=0.1,
                toi_period=ref_period,
                N_peaks=3,
                relative_height=10,
                FAP_levels=[0.1, 0.01, 0.001],
                verbose=self.verbose,
                ax=ax,
            )
            self.Prot_ls = self.gls.best["P"]
            if math.isnan(self.gls.power.max()):
                logger.error("GLS power is NaN, switching to astropy's Lomb-Scargle...")
                ax.clear()
                self.pg_method = "lombscargle"
                pg = plot_periodogram(
                    self.raw_lc[~self.tmask],
                    toi_period=ref_period,
                    method="lombscargle",
                    verbose=self.verbose,
                    ax=ax,
                )
                self.Prot_ls = pg.period_at_max_power.value
                # self.gls's periodogram is all-NaN (failed) -- append_tls_results()
                # reads self.gls.hpstat directly, so leaving it set would save a
                # degenerate Prot_gls/amp_gls/power_gls (observed landing at
                # Porb_max, a boundary artifact) under the GLS field names while
                # the plot and self.Prot_ls correctly show the astropy fallback
                # period. Clear it so those fields are honestly left unset
                # instead of silently wrong.
                self.gls = None
            else:
                pg = self.raw_lc[~self.tmask].to_periodogram(method="lombscargle")
        else:
            self.gls = None
            pg = plot_periodogram(
                self.raw_lc[~self.tmask],
                toi_period=ref_period,
                method=self.pg_method,
                verbose=self.verbose,
                ax=ax,
            )
            self.Prot_ls = pg.period_at_max_power.value
        if self.verbose:
            logger.info("Appending TLS results...")
        self.append_tls_results()

        if self.verbose:
            logger.info("Plotting phase-folded light curve...")
        ax = axes.flatten()[2]
        label = f"data folded at Prot={self.Prot_ls:.2f} d\n"
        if self.ephem_source:
            label += f"(masked transits using {self.ephem_source} ephem)"
        _ = (
            self.raw_lc[~self.tmask]
            .fold(
                self.Prot_ls,
                normalize_phase=False,
                wrap_phase=self.Prot_ls / 2,
            )
            .scatter(
                ax=ax,
                c=self.raw_lc[~self.tmask].time.value,
                cmap=pl.get_cmap("Blues_r"),
                label=label,
                zorder=0,
                show_colorbar=False,
            )
        )
        # Evaluate the model sinusoid at the *fold* period. In the GLS path the
        # fold period (self.Prot_ls) comes from the Gls object, which can differ
        # from this lightkurve periodogram's own peak frequency; without pinning
        # the frequency here pg.model() builds a sine at pg.period_at_max_power,
        # which then drifts across phase when folded at Prot_ls instead of
        # closing into a single clean sine wave.
        model_folded = pg.model(
            self.raw_lc[~self.tmask].time,
            frequency=(1.0 / self.Prot_ls) / u.day,
        ).fold(
            self.Prot_ls,
            normalize_phase=False,
            wrap_phase=self.Prot_ls / 2,
        )
        # Sort by phase so the model renders as a smooth curve, not zigzag lines
        sort_idx = np.argsort(model_folded.phase.value)
        ax.plot(
            model_folded.phase.value[sort_idx],
            model_folded.flux.value[sort_idx],
            label=f"{self.pg_method.upper()} model",
            color="r",
            ls="--",
            lw=2,
            zorder=1,
        )
        ax.set_xlabel("Rotation Phase [days]")

        single_layout = self.event_class in ("mono", "duo")
        if single_layout:
            if self.verbose:
                logger.info("Plotting single-event panels...")
            self._plot_single_event_panels(axes.flatten())
        else:
            if self.verbose:
                logger.info("Plotting TLS periodogram...")
            ax = axes.flatten()[4]
            ax = plot_tls(
                self.tls_results,
                toi_period=ref_period,
                period_min=self.Porb_min,
                period_max=self.Porb_max,
                ax=ax,
            )

            if self.verbose:
                logger.info("Plotting flattened light curve...")
            ax = axes.flatten()[3]
            self.flat_lc.scatter(ax=ax, label="flat")
            tmask2 = self.flat_lc.create_transit_mask(
                transit_time=self.tls_results.T0,
                period=self.tls_results.period,
                duration=self.tls_results.duration,
            )
            if self.has_ttv_layout:
                # Transits at their individually fitted times, not the linear ephemeris.
                ttv_lc_mask = self.ttv_transit_mask(self.flat_lc)
                self.flat_lc[ttv_lc_mask].scatter(ax=ax, color="r", label="transit (fitted times)")
            else:
                self.flat_lc[tmask2].scatter(ax=ax, color="r", label="transit")

        if self.verbose:
            logger.info("Plotting TPF...")
        ax = axes.flatten()[5]
        if not self.neighbor_check or not self._load_tpf():
            ax.text(
                0.5,
                0.5,
                "Neighbor check skipped\n(no TESS TPF for this target)",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.set_axis_off()
        else:
            ny, nx = self.tpf.flux.shape[1:]
            diag = np.sqrt(nx**2 + ny**2)
            fov_rad = (0.4 * diag * TESS_pix_scale).to(u.arcmin).round(2)
            tab = Catalogs.query_region(self.target_coord, radius=fov_rad, catalog="gaiadr3")
            self.gaia_sources = tab.to_pandas()
            # RUWE (renormalised unit weight error) of the target from Gaia DR3.
            # Values >~1.4 flag a poor single-star astrometric fit — often an
            # unresolved companion or blend. Matched to the target by source_id.
            self.gaia_ruwe = self._extract_ruwe(self.gaia_sources, self.gaiaid)
            self.tls_results["gaia_ruwe"] = float(self.gaia_ruwe)
            if len(self.gaia_sources) > 1:
                sep = self.gaia_sources.sort_values(by="distance", ascending=True)["distance"]
                self.nearby_star_sep = sep.iloc[1] * u.arcmin
            else:
                self.nearby_star_sep = None

            try:
                if self.verbose:
                    logger.info("Querying DSS data...")
                hdu = get_dss_data(
                    ra=self.target_coord.ra.deg,
                    dec=self.target_coord.dec.deg,
                    survey=self.archival_survey,
                    width=fov_rad.value,
                    height=fov_rad.value,
                )
                if hdu is None:
                    raise ValueError("DSS archival image not found")
                ax.remove()
                ax = fig.add_subplot(3, 3, 6, projection=WCS(hdu.header))
                _ = plot_gaia_sources_on_survey(
                    tpf=self.tpf,
                    target_gaiaid=self.gaiaid,
                    hdu=hdu,
                    fov_rad=fov_rad,
                    gaia_sources=self.gaia_sources,
                    kmax=1,
                    depth=_get_frac_depth(self.tls_results),
                    sap_mask=self.sap_mask,
                    aper_radius=2,
                    survey=self.archival_survey,
                    show_simbad=self.show_simbad,
                    verbose=self.verbose,
                    ax=ax,
                )
            except (OSError, FileNotFoundError, ValueError) as e:
                logger.warning(f"Archival image failed ({e}); plotting TPF instead.")
                _ = plot_gaia_sources_on_tpf(
                    tpf=self.tpf,
                    target_gaiaid=self.gaiaid,
                    fov_rad=fov_rad,
                    gaia_sources=self.gaia_sources,
                    kmax=1,
                    depth=_get_frac_depth(self.tls_results),
                    sap_mask=self.sap_mask,
                    aper_radius=2,
                    cmap="viridis",
                    dmag_limit=8,
                    verbose=self.verbose,
                    ax=ax,
                )

        if not single_layout and self.has_ttv_layout:
            if self.verbose:
                logger.info("Plotting TTV panels...")
            self._plot_ttv_panels(axes.flatten())
        elif not single_layout:
            if self.verbose:
                logger.info("Plotting odd-even transit...")
            ax = axes.flatten()[6]
            _ = plot_odd_even_transit(
                self.fold_lc,
                self.tls_results,
                bin_mins=10,
                markersize=6,
                ax=ax,
                phase_xlim=self.phase_xlim,
            )

            if self.verbose:
                logger.info("Plotting secondary eclipse...")
            ax = axes.flatten()[7]
            _ = plot_secondary_eclipse(
                self.flat_lc,
                self.tls_results,
                tmask2,
                bin_mins=10,
                markersize=6,
                ax=ax,
                phase_xlim=self.phase_xlim,
            )

        if self.verbose:
            logger.info("Creating summary panel...")
        ax = axes.flatten()[8]
        self._render_summary_panel(ax)
        title = ""
        if self.toiid is not None:
            title = f"TOI {self.toiid} | "
        title += f"TIC {self.ticid} | sector {self._sector_display()} | "
        lctype = (
            f"{self.pipeline.upper()}/{self.flux_type}"
            if self.pipeline == "spoc"
            else self.pipeline.upper()
        )
        title += f"{lctype.upper()} lightcurve"
        fig.suptitle(title, y=1.0, fontsize=20)

        if (self.outdir is not None) & (not Path(self.outdir).exists()):
            Path(self.outdir).mkdir()
            logger.info(f"Created output directory: {self.outdir}.")

        fp = self.check_output_file_exists()
        png_file = fp.with_suffix(".png")
        if self.savefig:
            fig.savefig(png_file, dpi=150, bbox_inches="tight")
            logger.info(f"Saved: {png_file}")
            self._save_single_event_figure(fp)
            self._save_ttv_figures(fp)

        h5_file = Path(self.outdir, fp.name + "_tls").with_suffix(".h5")
        if self.savetls:
            h5io.save(h5_file, self.tls_results)
            logger.info(f"Saved: {h5_file}")

        if (self.savefig or self.savetls) and len(
            getattr(self, "iterative_tls_results", None) or []
        ) > 1:
            self._save_iterative_candidates(fp)

        self.timer_end = timer()
        elapsed_time = self.timer_end - self.timer_start
        logger.info(f"#----------Runtime: {elapsed_time:.2f} s----------#\n")
        if not self.show_plot:
            pl.clf()
        if return_fig_and_paths:
            return fig, str(png_file), str(h5_file)
        else:
            return fig


if __name__ == "__main__":
    try:
        ql = TessQuickLook(
            "TOI-6043",
            sector=-1,
            # pipeline="qlp",
            # exptime=1800,
            pg_method="gls",
            savefig=True,
            savetls=True,
            overwrite=True,
            plot=True,
            flatten_method="biweight",
            gp_kernel="matern",  # squared_exp, matern, periodic, periodic_auto
            # gp_kernel_size=0.1,
            # window_length=0.1,
            # Porb_min=2,
            # Porb_max=4,
            outdir="../tests",
            # author='cdips'
            verbose=True,
        )
        if False:
            import cProfile
            import pstats

            profiler = cProfile.Profile()
            profiler.enable()
            fig = ql.plot_tql()
            profiler.disable()

            # Print stats to the terminal
            stats = pstats.Stats(profiler)
            stats.sort_stats("time").print_stats(10)

        if False:
            import tracemalloc

            tracemalloc.start()

            # Run the function or method
            fig = ql.plot_tql()

            # Take a snapshot and find any unclosed resources
            snapshot = tracemalloc.take_snapshot()
            for stat in snapshot.statistics("lineno"):
                logger.debug(stat)

        else:
            fig = ql.plot_tql()

        # warnings.resetwarnings()

    except Exception:
        # Get current system exception
        ex_type, ex_value, ex_traceback = sys.exc_info()
        # Extract unformatter stack traces as tuples
        trace_back = traceback.extract_tb(ex_traceback)

        logger.error(f"\nException type: {ex_type.__name__}")
        logger.error(f"Exception message: {ex_value}")
        # pdb.post_mortem(ex_traceback)
        # Format stacktrace
        for trace in trace_back:
            logger.error(f"Line : {trace[1]}")
            logger.error(f"Func : {trace[2]}")
            # logger.info(f"Message : {trace[3]}")
            logger.error(f"File : {trace[0]}")
        # traceback.print_exc()
