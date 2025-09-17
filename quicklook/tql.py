"""
1. See:
* https://ui.adsabs.harvard.edu/abs/2021ascl.soft01011R/abstract
* https://deep-lightcurve.readthedocs.io/en/latest/notebooks/Quickstart.html
* https://ui.adsabs.harvard.edu/abs/2022MNRAS.516.4432M/abstract
2. Add momentum dumps as in TESSLatte:
https://github.com/noraeisner/LATTE/blob/7ac35c8a51949345bc076fd30a456e74fce70c51/LATTE/LATTEutils.py#L3501C13-L3501C63
3. Add RUWE in plots
4. ingest functions from target.py:
http://localhost:9995/lab/workspaces/auto-q/tree/chronos/chronos/target.py

"""

import sys
import math
import traceback
import textwrap
import warnings
from pathlib import Path
from time import time as timer
from loguru import logger
from quicklook.compat import get_data_path
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
import lightkurve as lk
import flammkuchen as fk
from quicklook.utils import (
    get_tfop_info,
    get_params_from_tfop,
    TESS_TIME_OFFSET,
    TESS_pix_scale,
)
from quicklook.gls import Gls
from quicklook.plot import (
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

# FITSFixedWarning: 'datfix' made the change 'Invalid time in DATE-OBS
warnings.filterwarnings("ignore", category=Warning, message=".*datfix.*")
warnings.filterwarnings("ignore", category=Warning, message=".*obsfix.*")


__all__ = ["TessQuickLook"]

DATA_PATH = get_data_path("quicklook").joinpath("../data")
simbad_obj_list_file = Path(DATA_PATH, "simbad_obj_types.csv")
use_style("science")


class TessQuickLook:
    def __init__(
        self,
        target_name: str,
        sector=-1,
        pipeline: str = "SPOC",
        flux_type="pdcsap",
        exptime: float = None,
        pg_method: str = "gls",
        flatten_method: str = "biweight",
        gp_kernel: str = "matern",
        gp_kernel_size: float = 1,
        window_length: float = None,
        edge_cutoff: float = 0.1,
        sigma_clip_raw: tuple = None,
        sigma_clip_flat: tuple = None,
        custom_ephem: list = None,
        mask_ephem: bool = False,
        Porb_limits: tuple = None,
        archival_survey="dss1",
        show_plot: bool = True,
        verbose: bool = True,
        savefig: bool = False,
        savetls: bool = False,
        overwrite: bool = False,
        quality_bitmask: str = "default",
        outdir: str = ".",
    ):
        # start timer
        self.timer_start = timer()
        self.target_name = target_name
        logger.info(f"Generating quicklook for {self.target_name}...")
        self.verbose = verbose
        self.show_plot = show_plot
        self.tfop_info = get_tfop_info(target_name)
        self.parse_tfop_info()
        self.simbad_obj_type = self.get_simbad_obj_type()
        self.flux_type = flux_type
        self.exptime = exptime
        self.pg_method = pg_method
        self.sigma_clip_raw = sigma_clip_raw
        self.sigma_clip_flat = sigma_clip_flat
        self.quality_bitmask = quality_bitmask
        self.raw_lc = self.get_lc(
            author=pipeline,
            sector=sector,
            exptime=self.exptime,  # cadence=cadence
        )
        self.overwrite = overwrite
        self.outdir = outdir
        self.mask_ephem = mask_ephem
        _ = self.check_output_file_exists()
        self.flatten_method = flatten_method
        self.gp_kernel = (
            gp_kernel  # squared_exp, matern, periodic, periodic_auto
        )
        self.gp_kernel_size = gp_kernel_size
        self.edge_cutoff = edge_cutoff
        self.custom_ephem = custom_ephem
        self.parse_custom_ephem()

        if window_length is None:
            self.window_length = (
                self.tfop_dur[0] * 3
                if (self.tfop_dur is not None)
                and (self.tfop_dur[0] * 3 >= 0.1)
                else 0.5
            )
        else:
            self.window_length = window_length

        self.tmask = self.get_transit_mask()
        err_msg = "No masked transits"
        if self.tfop_epoch is not None and self.tmask.sum() == 0:
            logger.error(f"Error: {err_msg}")
            sys.exit()
        if self.mask_ephem:
            if self.verbose:
                logger.info(
                    f"Masking transits in raw lightcurve using {self.ephem_source} ephem..."
                )
            self.raw_lc = self.raw_lc[~self.tmask]
            # update tmask
            self.tmask = self.get_transit_mask()
        # self.flat_lc, self.trend_lc = self.raw_lc.flatten(return_trend=True)
        self.flat_lc, self.trend_lc = self.flatten_raw_lc()
        self.Porb_min = 0.1 if Porb_limits is None else Porb_limits[0]
        self.Porb_max = (
            (max(self.flat_lc.time.value) - min(self.flat_lc.time.value)) / 2
            if Porb_limits is None
            else Porb_limits[1]
        )
        self.run_tls()
        self.fold_lc = self.flat_lc.fold(
            period=self.tls_results.period,
            epoch_time=self.tls_results.T0,
            normalize_phase=False,
            wrap_phase=self.tls_results.period / 2,
        )
        self.savefig = savefig
        self.savetls = savetls
        self.archival_survey = archival_survey

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
                    args.append(f"{key}=({coord.replace(' ',',')})")
                elif val is not None:
                    args.append(f"{key}={val}")
        # Join the args with commas
        args = ", ".join(args)
        return f"{type(self).__name__}({args})"

    def parse_tfop_info(self):
        """
        Parse the TFOP info to get the star names, Gaia name, Gaia ID, and target coordinates.
        """
        self.star_names = np.array(
            self.tfop_info.get("basic_info")["star_names"].split(", ")
        )
        if self.verbose:
            print("Catalog names:")
            for n in self.star_names:
                print(f"\t{n}")
        self.gaia_name = self.star_names[
            np.array([i[:4].lower() == "gaia" for i in self.star_names])
        ][0]
        self.gaiaid = int(self.gaia_name.split()[-1])
        ra, dec = (
            self.tfop_info.get("coordinates")["ra"],
            self.tfop_info.get("coordinates")["dec"],
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
        self.ticid = int(self.tfop_info.get("basic_info")["tic_id"])
        if self.ticid is not None:
            self.query_name = f"TIC{self.ticid}"
        else:
            self.query_name = self.target_name.replace("-", " ")

    def parse_custom_ephem(self):
        """
        Parse the custom ephem to get the ephemeris source and parameters.
        """
        if self.custom_ephem:
            self.ephem_source = "custom"
            err_msg = (
                "Custom ephem must be a tuple: (t0,t0err,P,Perr,t14,t14err)"
            )
            if len(self.custom_ephem) != 6:
                logger.error(f"Error: {err_msg}")
                sys.exit()
            if self.verbose:
                logger.info(
                    f"Using ephemeris mask:\nP={self.custom_ephem[0]}d\nt0={self.custom_ephem[2]}BJD\nt14={self.custom_ephem[4]}d"
                )
            # TODO: using tfop in variable name is misleading
            if self.custom_ephem[0] > TESS_TIME_OFFSET:
                if self.verbose:
                    logger.info(
                        f"Custom transit epoch given in JD. Converting to BTJD = JD-{TESS_TIME_OFFSET:,}."
                    )
                self.custom_ephem[0] -= TESS_TIME_OFFSET
            self.tfop_epoch = (self.custom_ephem[0], self.custom_ephem[1])
            self.tfop_period = (self.custom_ephem[2], self.custom_ephem[3])
            if self.custom_ephem[4] > 1:
                if self.verbose:
                    logger.info(
                        "Custom transit duration given in hours. Converting to days."
                    )
                self.custom_ephem[4] /= 24
                self.custom_ephem[5] /= 24
            self.tfop_dur = (self.custom_ephem[4], self.custom_ephem[5])
            self.tfop_depth = None
        else:
            # use tfop ephem if available
            (
                self.tfop_epoch,
                self.tfop_period,
                self.tfop_dur,
                self.tfop_depth,
            ) = (
                self.get_toi_ephem()
                if len(self.tfop_info.get("planet_parameters")) != 0
                else (None, None, None, None)
            )
            self.ephem_source = "tfop" if self.tfop_epoch is not None else None

    def check_output_file_exists(self):
        name = self.target_name.replace(" ", "")
        if name.lower()[:3] == "toi":
            name = f"TOI{str(self.toiid).zfill(4)}"
        lctype = self.flux_type if self.pipeline == "spoc" else self.pipeline
        fp = Path(
            self.outdir,
            f"{name}_s{str(self.sector).zfill(2)}_{lctype}_{self.cadence[0]}c",
        )
        if self.mask_ephem:
            fp = fp.with_stem(fp.stem + "_mask_ephem")
        png_file = fp.with_suffix(".png")
        if png_file.exists() and not self.overwrite:
            raise FileExistsError(
                f"{png_file} already exists! Set overwrite=True."
            )
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
                    logger.info(
                        f"Simbad classifies {self.target_name} as {oid}={desc}!"
                    )
                    logger.info("***" * 15)
                else:
                    logger.info(
                        f"Simbad classifies {self.target_name} as {oid}={desc}!"
                    )

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
        if kwargs.get("sector") is None:
            sector_orig = None
        else:
            sector_orig = kwargs.pop("sector")
            sector = None if sector_orig in ["all", -1] else int(sector_orig)
            kwargs["sector"] = sector

        if sector_orig == "all":
            # assure data comes from single pipeline
            if kwargs["exptime"] is None:
                logger.error("Error: Supply exptime.")
                sys.exit()

        # Search for light curves related to the target
        search_result_all_lcs = lk.search_lightcurve(self.query_name)
        err_msg = f"Search using '{self.query_name}' did not yield any lightcurve results."
        if len(search_result_all_lcs) == 0:
            logger.error(f"Error: {err_msg}")
            sys.exit()

        # Extract all available sectors
        self.all_sectors = self.get_unique_sectors(search_result_all_lcs)
        # Display available light curves
        cols = ["author", "mission", "t_exptime"]
        if self.verbose:
            print("All available lightcurves:")
            print(search_result_all_lcs.table.to_pandas()[cols])

        # Validate the requested sector
        if kwargs.get("sector") is None:
            if self.verbose:
                print(f"Available sectors: {self.all_sectors}")
        else:
            err_msg = (
                f"sector={kwargs.get('sector')} not in {self.all_sectors}"
            )
            if kwargs.get("sector") not in self.all_sectors:
                logger.error(f"Error: {err_msg}")
                sys.exit()

        # Validate the requested author
        if kwargs.get("author") is None:
            kwargs["author"] = "SPOC"
        else:
            all_authors = set(
                search_result_all_lcs.table["provenance_name"].tolist()
            )
            err_msg = f"author={kwargs.get('author')} not in {all_authors}"
            if kwargs.get("author").upper() not in all_authors:
                logger.error(f"Error: {err_msg}")
                sys.exit()
        self.pipeline = kwargs["author"].lower()
        self.all_pipelines = all_authors

        # Search for specific light curve with given parameters
        search_result = lk.search_lightcurve(self.query_name, **kwargs)
        err_msg = f"Search using '{self.query_name}' {kwargs} did not yield any lightcurve results."
        if len(search_result) == 0:
            logger.error(f"Error: {err_msg}")
            sys.exit()

        # Download and return light curve
        if sector_orig == "all":
            if self.verbose:
                print(f"Filtered lightcurves based on query ({kwargs}):")
                print(search_result.table.to_pandas()[cols])
            msg = f"Downloading all {kwargs.get('author')} lcs..."
            if self.verbose:
                logger.info(msg)
            filtered_sectors = self.get_unique_sectors(search_result)
            if len(filtered_sectors) <= 1:
                logger.error(
                    f"Only {len(filtered_sectors)} sector is available."
                )
                sys.exit()
            lc = search_result.download_all(
                quality_bitmask=self.quality_bitmask
            ).stitch()
            self.sector = self.all_sectors
            # import pdb; pdb.set_trace()
            if self.pipeline in ["spoc"]:
                exptime = int(lc.meta["EXPOSURE"] / 10) * 10
            else:
                # estimate exp time
                exptime = round(np.diff(lc.time.jd).mean() * 24 * 60 * 60, -2)
            msg = f"Downloaded all {kwargs.get('author')} (exp={exptime} s) lcs in sectors {', '.join([str(s) for s in self.all_sectors])}."
            if self.verbose:
                logger.info(msg)
        else:
            # Download the light curve for the specified sector
            msg = f"Downloading {kwargs.get('author').upper()} lc..."
            if self.verbose:
                logger.info(msg)
            idx = sector_orig if sector_orig == -1 else 0
            lc = search_result[idx].download(
                quality_bitmask=self.quality_bitmask
            )
            # import pdb; pdb.set_trace()
            if self.pipeline in ["spoc"]:
                exptime = int(lc.meta["EXPOSURE"] / 10) * 10
            else:
                # estimate exp time
                exptime = round(np.diff(lc.time.jd).mean() * 24 * 60 * 60, -2)
            msg = f"Downloaded {lc.meta['AUTHOR'].upper()} (exp={exptime} s) lc in sector {lc.sector}."
            if self.verbose:
                logger.info(msg)
            self.sector = lc.sector

        # Select flux type for SPOC data
        if lc.meta["AUTHOR"].lower() == "spoc":
            lc = lc.select_flux(self.flux_type + "_flux")

        # Set exposure time and cadence
        if self.exptime is None:
            self.exptime = exptime
        # assert self.exptime == search_result.exptime[idx].value
        self.cadence = "short" if self.exptime < 1800 else "long"

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
        all_sectors = []
        for i in search_result.table["mission"].tolist():
            x = i.split()
            if len(x) == 3:
                s = int(x[-1])
                all_sectors.append(s)
        unique_sectors = sorted(set(all_sectors))
        return unique_sectors

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
        err_msg = "No tpf files found."
        if len(search_result_all_tpfs) == 0:
            logger.error(f"Error: {err_msg}")
            sys.exit()

        cols = ["author", "mission", "t_exptime"]
        if self.verbose:
            print("All available TPFs:")
            print(search_result_all_tpfs.table.to_pandas()[cols])
        tpf_authors = search_result_all_tpfs.table.to_pandas()[
            "author"
        ].unique()
        if kwargs.get("author").upper() not in tpf_authors:
            if self.verbose:
                logger.error(
                    f"No TPF for {kwargs.get('author').upper()} pipeline."
                )
            kwargs["author"] = [
                tpf_authors[i] for i in range(len(tpf_authors)) if i != "QLP"
            ][0]
        if self.verbose:
            logger.info(f"Using {kwargs.get('author').upper()} TPF...")

        # Search using the specified author and sector
        search_result = lk.search_targetpixelfile(self.query_name, **kwargs)
        err_msg = f"Search using '{self.query_name}' {kwargs} "
        err_msg += "did not yield any TPF results."
        if len(search_result) == 0:
            logger.error(f"Error: {err_msg}")
            sys.exit()
        msg = "Downloading TPF..."
        if self.verbose:
            logger.info(msg)
        idx = sector_orig if sector_orig == -1 else 0
        tpf = search_result[idx].download(quality_bitmask=self.quality_bitmask)
        # FIXME: What is the correct tpf aperture for other pipeline?
        # author = tpf.meta['PROCVER'].split('-')[0]
        author = search_result.author[idx].upper()
        exptime = search_result.exptime[idx].value
        msg = f"Downloaded {author} (exp={exptime} s) TPF "
        msg += f"in sector {tpf.meta['SECTOR']}."
        if self.verbose:
            logger.info(msg)
        return tpf

    def get_tpf_tesscut(self):
        """Download a 15x15 TESS postage stamp.

        Returns
        -------
        tpf : lightkurve.targetpixelfile.TargetPixelFile
            The downloaded TESS postage stamp.
        """
        if self.sector is None:
            logger.error("Provide sector.")
            sys.exit()
        tpf = lk.search_tesscut(self.query_name, sector=self.sector).download(
            cutout_size=(15, 15), quality_bitmask=self.quality_bitmask
        )
        if tpf is None:
            logger.error("No results from Tesscut search.")
            sys.exit()
        # remove zeros
        zero_mask = (tpf.flux_err == 0).all(axis=(1, 2))
        if zero_mask.sum() > 0:
            tpf = tpf[~zero_mask]
        return tpf

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
        if self.verbose:
            logger.info(f"Querying ephemeris for {self.target_name}:")
        try:
            # Use TIC latest uploaded ephem as default
            planet_params = get_params_from_tfop(
                self.tfop_info, "planet_parameters"
            )
        except Exception as e:
            logger.error(e)
            # If latest uploaded ephem is not available, use the first one
            planet_params = get_params_from_tfop(
                self.tfop_info, "planet_parameters", idx=1
            )
        if self.verbose:
            print(f"Parameters for {planet_params['name']}:")

        # Initialize variables
        tfop_epoch = None
        tfop_period = None
        tfop_dur = None
        tfop_depth = None

        # Query values and errors
        for p in params:
            unit = "hr" if p == "dur" else "d"
            val = planet_params.get(p)
            val = float(val) if val else 0.1
            err = planet_params.get(p + "_e")
            err = float(err) if err else 0.1
            if self.verbose:
                print(f"{p}: {val}, {err} {unit}")
            if p == "epoch":
                tfop_epoch = np.array((val, err))
                tfop_epoch[0] -= TESS_TIME_OFFSET
            elif p == "per":
                tfop_period = np.array((val, err))
            elif p == "dur":
                tfop_dur = np.array((val, err)) / 24

        # Query depth
        d = planet_params.get("dep_p")
        de = planet_params.get("dep_p_e")
        d = float(d) if d and (d != "") else np.nan
        de = float(de) if de and (de != "") else np.nan
        if not math.isnan(d) or not math.isnan(de):
            tfop_depth = np.array((d, de)) / 1e3

        return (tfop_epoch, tfop_period, tfop_dur, tfop_depth)

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
        if math.isnan(np.median(self.flat_lc.flux_err.value)):
            flux_err = np.zeros_like(self.flat_lc.flux_err)
            flux_err += np.nanstd(self.flat_lc.flux)
        else:
            flux_err = self.flat_lc.flux_err.value
        # Run TLS
        self.tls_results = tls(
            self.flat_lc.time.value,
            self.flat_lc.flux.value,
            flux_err,
            verbose=self.verbose,
        ).power(
            period_min=self.Porb_min,  # Roche limit default
            period_max=self.Porb_max,
        )

    def init_gls(self):
        if self.pipeline == "pathos":
            # pathos do not have flux_err
            cols = ["time", "flux"]
        else:
            cols = ["time", "flux", "flux_err"]
        data = (
            self.raw_lc.to_pandas().reset_index()[cols][~self.tmask].values.T
        )
        if self.verbose:
            logger.info(
                "Estimating rotation period using Generalized Lomb-Scargle (GLS) periodogram..."
            )
        return Gls(
            data, Pbeg=self.Porb_min, Pend=self.Porb_max, verbose=self.verbose
        )

    def flatten_raw_lc(self):
        """
        Flatten the raw light curve using WOTAN.

        Returns
        -------
        flat_lc : lightkurve.lightcurve.TessLightCurve
            The flattened light curve.
        trend_lc : lightkurve.lightcurve.TessLightCurve
            The trend light curve.

        """
        if self.verbose:
            logger.info(
                f"Using wotan's {self.flatten_method} method to flatten raw lc."
            )
        wflat_lc, wtrend_lc = flatten(
            # Array of time values
            self.raw_lc.time.value,
            # Array of flux values
            self.raw_lc.flux.value,
            # The method to use for detrending
            method=self.flatten_method,
            # The kernel to use for the Gaussian process
            kernel=self.gp_kernel,
            # The size of the kernel (if applicable)
            kernel_size=self.gp_kernel_size,
            # The length of the filter window in units of ``time``
            window_length=self.window_length,
            # The fraction of the window to cut off at the edges
            edge_cutoff=self.edge_cutoff,
            # The tolerance for breaks in the data
            break_tolerance=1,
            # Return the trend as well
            return_trend=True,
            # Tuning parameter for the robust estimators
            cval=5.0,
        )
        if self.sigma_clip_flat is not None:
            # Apply sigma clipping to the flattened light curve
            msg = "Applying sigma clip on flattened lc with "
            msg += f"(lower,upper)=({self.sigma_clip_flat})"
            if self.verbose:
                logger.info(msg)
            idx = sigma_clip(
                wflat_lc,
                # The lower and upper sigma limits
                sigma_lower=self.sigma_clip_flat[0],
                sigma_upper=self.sigma_clip_flat[1],
            ).mask
        else:
            # No sigma clipping
            idx = np.zeros_like(wflat_lc, dtype=bool)
        # Get the flattened and trend light curves
        flat_lc, trend_lc = self.raw_lc.flatten(return_trend=True)
        # Replace flux values with that from wotan
        flat_lc = flat_lc[~idx]
        trend_lc = trend_lc[~idx]
        trend_lc.flux = wtrend_lc[~idx]
        flat_lc.flux = wflat_lc[~idx]
        return flat_lc, trend_lc

    def get_transit_mask(self):
        """
        Generate a mask for the transit based on the user-provided
        transit ephemeris or the ephemeris from the TFOP portal.

        Returns
        -------
        tmask : np.ndarray
            A boolean mask where the transit is True and the out-of-transit
            periods are False.
        """
        if np.all([self.tfop_epoch, self.tfop_period, self.tfop_dur]):
            # Use the user-provided transit ephemeris to create the mask
            tmask = self.raw_lc.create_transit_mask(
                transit_time=self.tfop_epoch[0],
                period=self.tfop_period[0],
                duration=self.tfop_dur[0],
            )
        else:
            # If no transit ephemeris is provided, create an empty mask
            tmask = np.zeros_like(self.raw_lc.time.value, dtype=bool)
        return tmask

    def make_summary_info(self):
        """
        Generate a summary string with the TLS results, stellar params,
        and other useful information.

        Returns
        -------
        msg : str
            A summary string with the results.
        """
        try:
            # Use the TIC stellar parameters as default
            star_params = get_params_from_tfop(
                self.tfop_info, name="stellar_parameters", idx=1
            )
        except Exception as e:
            logger.error(e)
            star_params = get_params_from_tfop(
                self.tfop_info, name="stellar_parameters"
            )
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
            except Exception as e:
                logger.error(e)
                # Set to NaN if there's an error
                params[name] = np.nan
            try:
                # Convert to float or int as needed
                params[name + "_e"] = (
                    int(float(star_params.get(name + "_e")))
                    if name == "teff"
                    else float(star_params.get(name + "_e"))
                )
            except Exception as e:
                logger.error(e)
                # Set to NaN if there's an error
                params[name + "_e"] = np.nan
        # Get the meta data
        meta = self.raw_lc.meta
        # Calculate the planet radius
        Rp = self.tls_results["rp_rs"] * params["srad"] * u.Rsun.to(u.Rearth)
        # If the pipeline is Spoc or Tess-Spoc, correct for dilution
        if self.pipeline in ["spoc", "tess-spoc"]:
            Rp_true = Rp * np.sqrt(1 + meta["CROWDSAP"])
        else:
            # Otherwise, use the raw radius
            Rp_true = Rp
        # Create the summary string
        msg = "\nCandidate Properties\n"
        msg += "-" * 30 + "\n"
        text = f"SDE={self.tls_results.SDE:.4f} (sector={self.sector}"
        text += f" in {self.all_sectors})"
        msg += "\n".join(textwrap.wrap(text, 60))
        msg += f"\nPeriod={self.tls_results.period:.4f}" + r"$\pm$"
        msg += f"{self.tls_results.period_uncertainty:.4f} d (TLS)"
        if self.tfop_period is not None:
            msg += f", {self.tfop_period[0]:.4f}" + r"$\pm$"
            msg += f"{self.tfop_period[1]:.4f} d ({self.ephem_source})\n"
        else:
            msg += "\n"
        msg += f"T0={self.tls_results.T0:.4f} "
        if self.tfop_period is not None:
            msg += f"(TLS), {self.tfop_epoch[0]:.4f}" + r"$\pm$"
            msg += f"{self.tfop_epoch[1]:.4f} "
            msg += f"BJD-{TESS_TIME_OFFSET} ({self.ephem_source})\n"
        else:
            msg += f"BJD-{TESS_TIME_OFFSET} (TLS)\n"
        msg += f"Duration={self.tls_results.duration*24:.2f} hr (TLS)"
        if self.tfop_dur is not None:
            msg += f", {self.tfop_dur[0]*24:.2f}" + r"$\pm$"
            msg += f"{self.tfop_dur[1]*24:.2f} hr ({self.ephem_source})\n"
        else:
            msg += "\n"
        msg += f"Depth={(1-self.tls_results.depth)*1e3:.2f} ppt (TLS)"
        if self.tfop_depth is not None:
            msg += f", {self.tfop_depth[0]:.1f}" + r"$\pm$"
            msg += f"{self.tfop_depth[1]:.1f} ppt (tfop)\n"
        else:
            msg += "\n"
        if (meta["FLUX_ORIGIN"].lower() == "pdcsap") or (
            meta["FLUX_ORIGIN"].lower() == "sap"
        ):
            # msg += f"Rp={Rp:.2f} " + r"R$_{\oplus}$" + "(diluted)" + " " * 5
            msg += f"Rp={Rp_true:.2f} " + r"R$_{\oplus}$ "
            msg += (
                f"= {Rp_true*u.Rearth.to(u.Rjup):.2f}"
                + r"R$_{\rm{Jup}}$"
                + "\n"
            )
        else:
            msg += f"Rp={Rp:.2f} " + r"R$_{\oplus}$" + "(diluted), "
            msg += f"Rp={Rp_true:.2f} " + r"R$_{\oplus}$" + "(undiluted)\n"
        msg += (
            f"Odd-Even mismatch={self.tls_results.odd_even_mismatch:.2f}"
            + r"$\sigma$"
        )
        msg += "\n" * 2
        msg += "Stellar Properties\n"
        msg += "-" * 30 + "\n"
        msg += (
            f"Rstar={params['srad']:.2f}"
            + r"$\pm$"
            + f"{params['srad_e']:.2f} "
            + r"R$_{\odot}$"
            + " " * 5
        )
        msg += (
            f"Teff={params.get('teff')}"
            + r"$\pm$"
            + f"{params.get('teff_e')} K"
            + "\n"
        )
        msg += (
            f"Mstar={params['mass']:.2f}"
            + r"$\pm$"
            + f"{params['mass_e']:.2f} "
            + r"M$_{\odot}$"
            + " " * 5
        )
        msg += (
            f"logg={params['logg']:.2f}"
            + r"$\pm$"
            + f"{params['logg_e']:.2f} cgs\n"
        )
        msg += f"Rotation period={self.Prot_ls:.2f} d" + " " * 5
        per = 2 * np.pi * params["srad"] * u.Rsun.to(u.km)
        t = self.Prot_ls * u.day.to(u.second)
        msg += f"Vsini={per/t:.2f} km/s\n"
        msg += f"Gaia DR2 ID={self.gaiaid}\n"
        # msg += f"TIC ID={self.ticid}" + " " * 5
        coords = self.target_coord.to_string("decimal").split()
        msg += f"RA,Dec={float(coords[0]), float(coords[1])}\n"
        msg += (
            f"Distance={params['dist']:.1f}"
            + r"$\pm$"
            + f"{params['dist_e']:.1f} pc"
        )
        mags = self.tfop_info["magnitudes"][0]
        msg += f", {mags['band']}mag={float(mags['value']):.1f}\n"
        # msg += f"GOF_AL={astrometric_gof_al:.2f} (hints binarity if >20)\n"
        # D = gp.astrometric_excess_noise_sig
        # msg += f"astro. excess noise sig={D:.2f} (hints binarity if >5)\n"
        if self.simbad_obj_type is not None:
            msg += f"Simbad Object: {self.simbad_obj_type}"
        # msg += f"met={feh:.2f}"+r"$\pm$"+f"{feh_err:.2f} dex " + " " * 6
        return msg

    def append_tls_results(self):
        """
        Append TLS results to the TLS results dictionary.

        This will add the raw and flattened light curves, period limits, and
        TFOP parameters to the TLS results dictionary.

        Returns
        -------
        None
        """
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

        # Append the TFOP parameters
        self.tls_results["period_tfop"] = self.tfop_period
        self.tls_results["T0_tfop"] = self.tfop_epoch
        self.tls_results["duration_tfop"] = self.tfop_dur
        self.tls_results["depth_tfop"] = self.tfop_depth

        # Append the Gaia ID, TIC ID, and TOI ID
        self.tls_results["gaiaid"] = int(self.gaiaid)
        self.tls_results["ticid"] = int(self.ticid)
        self.tls_results["toiid"] = self.toiid
        self.tls_results["sector"] = self.sector

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

        # Append the Simbad object type
        if self.simbad_obj_type is not None:
            self.tls_results["simbad_obj"] = self.simbad_obj_type

    def plot_tql(self, **kwargs: dict) -> pl.Figure:
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
        ax = axes.flatten()[0]
        self.raw_lc.scatter(ax=ax, label=f"raw (exp={self.exptime} s)")

        if self.verbose:
            logger.info("Plotting trend...")
        label = f"baseline trend\nmethod={self.flatten_method}(window_size={self.window_length:.2f})"
        self.trend_lc.plot(ax=ax, color="r", lw=2, label=label)

        if self.verbose:
            logger.info("Running Lomb-Scargle periodogram...")
        ax = axes.flatten()[1]
        if self.pg_method == "gls":
            self.gls = self.init_gls()
            ax = plot_gls_periodogram(
                self.gls,
                offset=0.1,
                N_peaks=3,
                relative_height=10,
                FAP_levels=[0.1, 0.01, 0.001],
                verbose=self.verbose,
                ax=ax,
            )
            pg = self.raw_lc[~self.tmask].to_periodogram(method="lombscargle")
            self.Prot_ls = self.gls.best["P"]
            if math.isnan(self.gls.power.max()):
                logger.error(
                    "GLS power is NaN, switching to astropy's Lomb-Scargle..."
                )
                ax.clear()
                self.pg_method = "lombscargle"
                pg = plot_periodogram(
                    self.raw_lc[~self.tmask],
                    method="lombscargle",
                    verbose=self.verbose,
                    ax=ax,
                )
                self.Prot_ls = pg.period_at_max_power.value
        else:
            self.gls = None
            pg = plot_periodogram(
                self.raw_lc[~self.tmask],
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
        _ = (
            pg.model(self.raw_lc[~self.tmask].time)
            .fold(
                self.Prot_ls,
                normalize_phase=False,
                wrap_phase=self.Prot_ls / 2,
            )
            .plot(
                label=f"{self.pg_method.upper()} model",
                color="r",
                ls="--",
                lw=2,
                zorder=1,
                ax=ax,
            )
        )
        ax.set_xlabel("Rotation Phase [days]")

        if self.verbose:
            logger.info("Plotting TLS periodogram...")
        ax = axes.flatten()[4]
        ax = plot_tls(
            self.tls_results,
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
        self.flat_lc[tmask2].scatter(ax=ax, color="r", label="transit")

        if self.verbose:
            logger.info("Plotting TPF...")
        ax = axes.flatten()[5]
        if self.pipeline in [
            "cdips",
            "gsfc-eleanor-lite",
        ]:
            err_msg = "Pipeline to be added soon."
            logger.info(err_msg)
            raise NotImplementedError(err_msg)
        elif self.pipeline in ["qlp", "cdips", "tasoc", "pathos", "tglc"]:
            if self.verbose:
                logger.info("Getting TPF with tesscut...")
            self.tpf = self.get_tpf_tesscut()
            self.sap_mask = "square"
        else:
            self.tpf = self.get_tpf(
                sector=self.sector,
                author=self.pipeline,
            )
            self.sap_mask = "pipeline"
        try:
            if self.verbose:
                logger.info("Querying DSS data...")
            ny, nx = self.tpf.flux.shape[1:]
            diag = np.sqrt(nx**2 + ny**2)
            fov_rad = (0.4 * diag * TESS_pix_scale).to(u.arcmin).round(2)
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
                gaia_sources=None,
                kmax=1,
                depth=1 - self.tls_results.depth,
                sap_mask=self.sap_mask,
                aper_radius=2,
                survey=self.archival_survey,
                fov_rad=fov_rad,
                verbose=self.verbose,
                ax=ax,
            )
        except Exception as e:
            logger.error(f"Error: {e}")
            logger.info(
                "Querying archival image failed. Plotting TPF instead."
            )
            _ = plot_gaia_sources_on_tpf(
                tpf=self.tpf,
                target_gaiaid=self.gaiaid,
                gaia_sources=None,
                kmax=1,
                depth=1 - self.tls_results.depth,
                sap_mask=self.sap_mask,
                aper_radius=2,
                cmap="viridis",
                dmag_limit=8,
                verbose=self.verbose,
                ax=ax,
            )

        if self.verbose:
            logger.info("Plotting odd-even transit...")
        ax = axes.flatten()[6]
        _ = plot_odd_even_transit(
            self.fold_lc, self.tls_results, bin_mins=10, ax=ax
        )

        if self.verbose:
            logger.info("Plotting secondary eclipse...")
        ax = axes.flatten()[7]
        _ = plot_secondary_eclipse(
            self.flat_lc, self.tls_results, tmask2, bin_mins=10, ax=ax
        )

        if self.verbose:
            logger.info("Plotting summary panel...")
        ax = axes.flatten()[8]
        ax.axis([0, 10, 0, 10])
        msg = self.make_summary_info()
        ax.text(-1, 11, msg, ha="left", va="top", fontsize=10, wrap=True)
        ax.axis("off")
        title = ""
        if self.toiid is not None:
            title = f"TOI {self.toiid} | "
        title += f"TIC {self.ticid} | sector {self.sector} | "
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
        if self.target_name.lower()[:4] == "gaia":
            name = self.target_name.replace(" ", "_")
        else:
            name = self.target_name.replace(" ", "")
        fp = Path(
            self.outdir,
            f"{name}_s{str(self.sector).zfill(2)}_{self.pipeline}_{self.cadence[0]}c",
        )
        fp = self.check_output_file_exists()
        png_file = fp.with_suffix(".png")
        if self.savefig:
            fig.savefig(png_file, dpi=150, bbox_inches="tight")
            logger.info("Saved: ", png_file)

        if self.savetls:
            h5_file = Path(self.outdir, fp.name + "_tls").with_suffix(".h5")
            fk.save(h5_file, self.tls_results)
            logger.info("Saved: ", h5_file)

        self.timer_end = timer()
        elapsed_time = self.timer_end - self.timer_start
        logger.info(f"#----------Runtime: {elapsed_time:.2f} s----------#\n")
        if not self.show_plot:
            pl.clf()
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
                print(stat)

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
