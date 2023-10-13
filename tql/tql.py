import sys
import math
import traceback
import textwrap
from pkg_resources import resource_filename
from pathlib import Path
from time import time as timer
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
from aesthetic.plot import set_style
from aesthetic.plot import savefig as save_figure
import flammkuchen as fk
from tql.utils import get_tfop_info, TESS_TIME_OFFSET, TESS_pix_scale
from tql.gls import Gls
from tql.plot import (
    get_dss_data,
    plot_gaia_sources_on_survey,
    plot_gaia_sources_on_tpf,
    plot_odd_even_transit,
    plot_secondary_eclipse,
    plot_tls,
    plot_periodogram,
    plot_gls_periodogram,
)

set_style("science")

__all__ = ["TessQuickLook"]

simbad_obj_list_file = Path(
    resource_filename("tql", "../data/simbad_obj_types.csv")
).resolve()


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
        ephem_mask: list = None,
        Porb_limits: tuple = None,
        archival_survey="dss1",
        plot: bool = True,
        savefig: bool = False,
        savetls: bool = False,
        overwrite: bool = False,
        outdir: str = ".",
    ):
        self.timer_start = timer()

        self.target_name = target_name
        print(f"Generating TQL for {self.target_name}...")
        self.tfop_info = get_tfop_info(target_name)
        self.star_names = np.array(
            self.tfop_info.get("basic_info")["star_names"].split(", ")
        )
        self.gaia_name = self.star_names[
            np.array([i[:4].lower() == "gaia" for i in self.star_names])
        ][0]
        self.gaiaid = int(self.gaia_name.split()[-1])
        ra, dec = (
            self.tfop_info.get("coordinates")["ra"],
            self.tfop_info.get("coordinates")["dec"],
        )
        self.target_coord = SkyCoord(ra=ra, dec=dec, unit="degree")

        if target_name.lower()[:3] == "toi":
            toiid = int(target_name.split("-")[-1])
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
            self.query_name = self.target_name.replace(" ", "")
        self.simbad_obj_type = self.get_simbad_obj_type()
        self.flux_type = flux_type
        self.exptime = exptime
        self.pg_method = pg_method
        self.sigma_clip_raw = sigma_clip_raw
        self.sigma_clip_flat = sigma_clip_flat
        self.raw_lc = self.get_lc(
            author=pipeline,
            sector=sector,
            exptime=self.exptime,  # cadence=cadence
        )
        self.overwrite = overwrite
        self.outdir = outdir
        _ = self.check_file_exists()
        self.flatten_method = flatten_method
        self.gp_kernel = (
            gp_kernel  # squared_exp, matern, periodic, periodic_auto
        )
        self.gp_kernel_size = gp_kernel_size
        self.edge_cutoff = edge_cutoff
        self.ephem_mask = ephem_mask
        _ = self.get_toi_ephem()
        if window_length is None:
            self.window_length = (
                self.tfop_dur[0] * 3
                if (self.tfop_dur is not None)
                and (self.tfop_dur[0] * 3 >= 0.1)
                else 0.5
            )
        else:
            self.window_length = window_length
        # self.flat_lc, self.trend_lc = self.raw_lc.flatten(return_trend=True)
        self.flat_lc, self.trend_lc = self.flatten_raw_lc()
        self.tmask = self.get_transit_mask()
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
        self.plot = plot
        self.savefig = savefig
        self.savetls = savetls
        self.archival_survey = archival_survey

    def __repr__(self):
        """Override to print a readable string representation of class"""

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
        args = []
        for key in self.__dict__.keys():
            val = self.__dict__.get(key)
            if key in included_args:
                if key == "target_coord":
                    # format coord
                    coord = self.target_coord.to_string("decimal")
                    args.append(f"{key}=({coord.replace(' ',',')})")
                elif val is not None:
                    args.append(f"{key}={val}")
        args = ", ".join(args)
        return f"{type(self).__name__}({args})"

    def check_file_exists(self):
        name = self.target_name.replace(" ", "")
        if name.lower()[:3] == "toi":
            name = f"TOI{str(self.toiid).zfill(4)}"
        lctype = self.flux_type if self.pipeline == "spoc" else self.pipeline
        fp = Path(
            self.outdir,
            f"{name}_s{str(self.sector).zfill(2)}_{lctype}_{self.cadence[0]}c",
        )
        png_file = fp.with_suffix(".png")
        if png_file.exists() and not self.overwrite:
            raise FileExistsError(
                f"{png_file} already exists! Set overwrite=True."
            )
        return fp

    def get_simbad_obj_type(self):
        """See also: https://simbad.cds.unistra.fr/guide/otypes.htx"""
        Simbad.add_votable_fields("otype")
        try:
            r = Simbad.query_object(self.target_name.replace("-", ""))
            category = r.to_pandas().squeeze()["OTYPE"]
            df = pd.read_csv(simbad_obj_list_file)
            dd = df.query("Id==@category")
            desc = dd["Description"].squeeze()
            oid = dd["Id"].squeeze()
            if dd["Description"].str.contains("(?i)binary").any():
                print("***" * 5)
                print(f"Simbad classifies {self.target_name} as {oid}={desc}!")
                print("***" * 5)
            return desc
        except Exception as e:
            print(f"Simbad cannot resolve {self.target_name}.\n{e}")
            return None

    def get_lc(self, **kwargs: dict) -> lk.TessLightCurve:
        """
        Parameters
        ----------
        target_name : str
        kwargs : dict
            radius, sector, author, cadence, exptime

        Returns
        -------
        lc : lk.Lightcurve
        """
        if kwargs.get("sector") is None:
            sector_orig = None
        else:
            sector_orig = kwargs.pop("sector")
            sector = None if sector_orig in ["all", -1] else sector_orig
            kwargs["sector"] = sector

        # get all info
        search_result_all_lcs = lk.search_lightcurve(self.query_name)
        errmsg = f"Search using '{self.query_name}' "
        errmsg += f"did not yield any lightcurve results."
        assert len(search_result_all_lcs) > 0, errmsg
        cols = ["author", "mission", "t_exptime"]
        print("All available lightcurves:")
        print(search_result_all_lcs.table.to_pandas()[cols])

        all_sectors = []
        for i in search_result_all_lcs.table["mission"].tolist():
            x = i.split()
            if len(x) == 3:
                s = int(x[-1])
                all_sectors.append(s)
        all_sectors = sorted(set(all_sectors))

        if kwargs.get("sector") is None:
            print(f"Available sectors: {all_sectors}")
        else:
            err_msg = f"sector={kwargs.get('sector')} not in {all_sectors}"
            assert kwargs.get("sector") in all_sectors, err_msg

        if kwargs.get("author") is None:
            kwargs["author"] = "SPOC"
        else:
            all_authors = set(
                search_result_all_lcs.table["provenance_name"].tolist()
            )
            err_msg = f"author={kwargs.get('author')} not in {all_authors}"
            assert kwargs.get("author").upper() in all_authors, err_msg
        self.pipeline = kwargs["author"].lower()
        self.all_pipelines = all_authors

        # get specific lc in cache
        search_result = lk.search_lightcurve(self.query_name, **kwargs)
        errmsg = f"Search using '{self.query_name}' {kwargs} "
        errmsg += "did not yield any lightcurve results."
        assert len(search_result) > 0, errmsg
        # print(search_result)

        self.all_sectors = all_sectors
        # FIXME: no use case for running tls on all sectors
        if sector_orig == "all":
            lcs = search_result.download_all().stitch()
            print(
                f"Downloaded all {kwargs.get('author')} lc \
                  in sectors {', '.join([str(s) for s in all_sectors])}."
            )
            self.sector = all_sectors
            return lcs
        else:
            idx = sector_orig if sector_orig == -1 else 0
            lc = search_result[idx].download()
            exptime = search_result.exptime[idx].value
            msg = f"\nDownloaded {kwargs.get('author').upper()} "
            msg += f"(exp={exptime} s) lc in sector {lc.sector}.\n"
            print(msg)
            self.sector = lc.sector
        if lc.meta["AUTHOR"].lower() == "spoc":
            lc = lc.select_flux(self.flux_type + "_flux")
        if self.exptime is None:
            self.exptime = exptime
        self.cadence = "short" if self.exptime < 1800 else "long"
        if self.sigma_clip_raw is not None:
            print("Applying sigma clip on raw lc with ")
            print(f"(lower,upper)={self.sigma_clip_raw}")
            return lc.normalize().remove_outliers(
                sigma_lower=self.sigma_clip_raw[0],
                sigma_upper=self.sigma_clip_raw[1],
            )
        else:
            return lc.normalize()

    def get_tpf(self, **kwargs: dict) -> lk.targetpixelfile.TargetPixelFile:
        if kwargs.get("sector") is None:
            sector_orig = None
        else:
            sector_orig = kwargs.pop("sector")
            sector = None if sector_orig in ["all", -1] else sector_orig
            kwargs["sector"] = sector

        if kwargs.get("author") is None:
            kwargs["author"] = "SPOC"

        search_result_all_tpfs = lk.search_targetpixelfile(self.query_name)
        errmsg = "No tpf files found."
        assert len(search_result_all_tpfs) > 0, errmsg

        cols = ["author", "mission", "t_exptime"]
        print("All available TPFs:")
        print(search_result_all_tpfs.table.to_pandas()[cols])
        tpf_authors = search_result_all_tpfs.table.to_pandas()[
            "author"
        ].unique()
        if kwargs.get("author").upper() not in tpf_authors:
            print(f"No TPF for {kwargs.get('author').upper()} pipeline.")
            kwargs["author"] = [
                tpf_authors[i] for i in range(len(tpf_authors)) if i != "QLP"
            ][0]
        print(f"\nUsing {kwargs.get('author').upper()} TPF.\n")

        search_result = lk.search_targetpixelfile(self.query_name, **kwargs)
        errmsg = f"Search using '{self.query_name}' {kwargs} "
        errmsg += f"did not yield any TPF results."
        assert len(search_result) > 0, errmsg
        idx = sector_orig if sector_orig == -1 else 0
        tpf = search_result[idx].download()
        # FIXME: What is the correct tpf aperture for other pipeline?
        # author = tpf.meta['PROCVER'].split('-')[0]
        author = search_result.author[idx].upper()
        exptime = search_result.exptime[idx].value
        msg = f"Downloaded {author.upper()} (exp={exptime} s) TPF "
        msg += f"in sector {tpf.meta['SECTOR']}."
        print(msg)
        return tpf

    def get_tpf_tesscut(self):
        """ """
        if self.sector is None:
            errmsg = "Provide sector."
            raise ValueError(errmsg)
        tpf = lk.search_tesscut(self.query_name, sector=self.sector).download(
            cutout_size=(15, 15)
        )
        assert tpf is not None, "No results from Tesscut search."
        # remove zeros
        zero_mask = (tpf.flux_err == 0).all(axis=(1, 2))
        if zero_mask.sum() > 0:
            tpf = tpf[~zero_mask]
        return tpf

    def get_toi_ephem(self, idx=None, params=["epoch", "per", "dur"]) -> list:
        print(f"Querying ephemeris for {self.target_name}:")
        params_dict = self.tfop_info.get("planet_parameters")
        if idx is None:
            idx = np.argmax(
                [len(params_dict[x]) for x in range(len(params_dict))]
            )
        planet_params = params_dict[idx]
        vals = []
        for p in params:
            val = planet_params.get(p)
            val = float(val) if val else 0.1
            err = planet_params.get(p + "_e")
            err = float(err) if err else 0.1
            print(f"{p}: {val}, {err}")
            vals.append((val, err))
        print("")
        if len(vals) > 0:
            self.tfop_epoch = np.array(vals[0]) - TESS_TIME_OFFSET
            self.tfop_period = np.array(vals[1])
            self.tfop_dur = np.array(vals[2]) / 24
        else:
            self.tfop_epoch = None
            self.tfop_period = None
            self.tfop_dur = None

        d = planet_params.get("dep_p")
        de = planet_params.get("dep_p_e")
        d = float(d) if d and (d != "") else np.nan
        de = float(de) if de and (de != "") else np.nan
        if not math.isnan(d) or not math.isnan(de):
            self.tfop_depth = np.array((d, de)) / 1e3
        else:
            self.tfop_depth = None
        return vals

    def run_tls(self):
        # CDIPS light curve has no flux err
        if math.isnan(np.median(self.flat_lc.flux_err.value)):
            flux_err = np.zeros_like(self.flat_lc.flux_err)
            flux_err += np.nanstd(self.flat_lc.flux)
        else:
            flux_err = self.flat_lc.flux_err.value
        # Run TLS
        self.tls_results = tls(
            self.flat_lc.time.value, self.flat_lc.flux.value, flux_err
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
        print(
            "Estimating rotation period using Generalized Lomb-Scargle (GLS) periodogram."
        )
        return Gls(data, Pbeg=self.Porb_min, Pend=self.Porb_max, verbose=True)

    def flatten_raw_lc(self):
        print(f"Using wotan's {self.flatten_method} method to flatten raw lc.")
        wflat_lc, wtrend_lc = flatten(
            self.raw_lc.time.value,  # Array of time values
            self.raw_lc.flux.value,  # Array of flux values
            method=self.flatten_method,
            kernel=self.gp_kernel,
            # FIXME: might be useful for method=gp
            kernel_size=self.gp_kernel_size,
            # The length of the filter window in units of ``time``
            window_length=self.window_length,
            edge_cutoff=self.edge_cutoff,
            break_tolerance=1,  # Split into segments at breaks longer than
            return_trend=True,
            cval=5.0,  # Tuning parameter for the robust estimators
        )
        if self.sigma_clip_flat is not None:
            msg = "Applying sigma clip on flattened lc with "
            msg += f"(lower,upper)=({self.sigma_clip_flat})"
            print(msg)
            idx = sigma_clip(
                wflat_lc,
                sigma_lower=self.sigma_clip_flat[0],
                sigma_upper=self.sigma_clip_flat[1],
            ).mask
        else:
            idx = np.zeros_like(wflat_lc, dtype=bool)
        flat_lc, trend_lc = self.raw_lc.flatten(return_trend=True)
        # replace flux values with that from wotan
        flat_lc = flat_lc[~idx]
        trend_lc = trend_lc[~idx]
        trend_lc.flux = wtrend_lc[~idx]
        flat_lc.flux = wflat_lc[~idx]
        return flat_lc, trend_lc

    def get_transit_mask(self):
        if np.all([self.tfop_epoch[0], self.tfop_period[0], self.tfop_dur[0]]):
            tmask = self.raw_lc.create_transit_mask(
                transit_time=self.tfop_epoch[0],
                period=self.tfop_period[0],
                duration=self.tfop_dur[0],
            )
            err_msg = "No masked transits"
            assert tmask.sum() > 0, err_msg
        else:
            tmask = np.zeros_like(self.raw_lc.time.value, dtype=bool)
        return tmask

    def make_summary_info(self):
        star_params = self.tfop_info["stellar_parameters"][1]
        params = {}
        param_names = ["srad", "mass", "teff", "logg", "dist"]
        for name in param_names:
            try:
                params[name] = float(star_params.get(name))
            except:
                params[name] = np.nan
            try:
                params[name + "_e"] = float(star_params.get(name + "_e"))
            except:
                params[name + "_e"] = np.nan
        meta = self.raw_lc.meta
        Rp = self.tls_results["rp_rs"] * params["srad"] * u.Rsun.to(u.Rearth)
        if self.pipeline in ["spoc", "tess-spoc"]:
            Rp_true = Rp * np.sqrt(1 + meta["CROWDSAP"])
        else:
            # FIXME: need to get dilution from other pipelines
            Rp_true = Rp
        msg = "\nCandidate Properties\n"
        msg += "-" * 30 + "\n"
        text = f"SDE={self.tls_results.SDE:.4f} (sector={self.sector}"
        text += f" in {self.all_sectors})"
        msg += "\n".join(textwrap.wrap(text, 60))
        msg += f"\nPeriod={self.tls_results.period:.4f}" + r"$\pm$"
        msg += f"{self.tls_results.period_uncertainty:.4f} d (TLS)"
        if self.tfop_period is not None:
            msg += f", {self.tfop_period[0]:.4f}" + r"$\pm$"
            msg += f"{self.tfop_period[1]:.4f} d (TFOP)\n"
        else:
            msg += "\n"
        msg += f"T0={self.tls_results.T0:.4f} "
        if self.tfop_period is not None:
            msg += f"(TLS), {self.tfop_epoch[0]:.4f}" + r"$\pm$"
            msg += f"{self.tfop_epoch[1]+TESS_TIME_OFFSET:.4f} "
            msg += f"BJD-{TESS_TIME_OFFSET} (TFOP)\n"
        else:
            msg += "BJD-{TESS_TIME_OFFSET} (TLS)\n"
        msg += f"Duration={self.tls_results.duration*24:.2f} hr (TLS)"
        if self.tfop_dur is not None:
            msg += f", {self.tfop_dur[0]*24:.2f}" + r"$\pm$"
            msg += f"{self.tfop_dur[1]*24:.2f} hr (TFOP)\n"
        else:
            msg += "\n"
        msg += f"Depth={(1-self.tls_results.depth)*1e3:.2f} ppt (TLS)"
        if self.tfop_depth is not None:
            msg += f", {self.tfop_depth[0]:.1f}" + r"$\pm$"
            msg += f"{self.tfop_depth[1]:.1f} ppt (TFOP)\n"
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
        msg += f"Gaia DR2 ID={self.gaiaid}"
        # msg += f"TIC ID={self.ticid}" + " " * 5
        msg += f", Tmag={meta['TESSMAG']:.2f}\n"
        msg += (
            f"Distance={params['dist']:.1f}"
            + r"$\pm$"
            + f"{params['dist_e']:.1f} pc\n"
        )
        # msg += f"GOF_AL={astrometric_gof_al:.2f} (hints binarity if >20)\n"
        # D = gp.astrometric_excess_noise_sig
        # msg += f"astro. excess noise sig={D:.2f} (hints binarity if >5)\n"
        msg += (
            f"Rstar={params['srad']:.2f}"
            + r"$\pm$"
            + f"{params['srad_e']:.2f} "
            + r"R$_{\odot}$"
            + " " * 5
        )
        msg += (
            f"Teff={int(params['teff'])}"
            + r"$\pm$"
            + f"{int(params['teff_e'])} K"
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
        msg += f"Rotation speed={per/t:.2f} km/s\n"
        if self.simbad_obj_type is not None:
            msg += f"Simbad Object: {self.simbad_obj_type}"
        # msg += f"met={feh:.2f}"+r"$\pm$"+f"{feh_err:.2f} dex " + " " * 6
        return msg

    def append_tls_results(self):
        self.tls_results["time_raw"] = self.raw_lc.time.value
        self.tls_results["flux_raw"] = self.raw_lc.flux.value
        self.tls_results["err_raw"] = self.raw_lc.flux_err.value
        self.tls_results["time_flat"] = self.flat_lc.time.value
        self.tls_results["flux_flat"] = self.flat_lc.flux.value
        self.tls_results["err_flat"] = self.flat_lc.flux_err.value
        self.tls_results["Porb_min"] = self.Porb_min
        self.tls_results["Porb_max"] = self.Porb_max
        self.tls_results["period_tfop"] = self.tfop_period
        self.tls_results["T0_tfop"] = self.tfop_epoch
        self.tls_results["duration_tfop"] = self.tfop_dur
        self.tls_results["depth_tfop"] = self.tfop_depth
        self.tls_results["gaiaid"] = self.gaiaid
        self.tls_results["ticid"] = self.ticid
        self.tls_results["toiid"] = self.toiid
        self.tls_results["sector"] = self.sector
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
        if self.simbad_obj_type is not None:
            self.tls_results["simbad_obj"] = self.simbad_obj_type

    def plot_tql(self, **kwargs: dict) -> pl.Figure:
        if kwargs.get("savefig"):
            self.savefig = kwargs.get("savefig")
        if kwargs.get("overwrite"):
            self.overwrite = kwargs.get("overwite")
        if kwargs.get("plot"):
            self.plot = kwargs.get("plot")

        fig, axes = pl.subplots(3, 3, figsize=(16, 12), tight_layout=True)

        # +++++++++++++++++++++ax: Raw + trend
        ax = axes.flatten()[0]
        self.raw_lc.scatter(ax=ax, label=f"raw (exp={self.exptime} s)")
        label = f"baseline trend\nmethod={self.flatten_method}"
        label += f"(window_size={self.window_length:.2f})"
        self.trend_lc.plot(ax=ax, color="r", lw=2, label=label)

        # +++++++++++++++++++++ax2 Lomb-scargle periodogram
        ax = axes.flatten()[1]

        if self.pg_method == "gls":
            self.gls = self.init_gls()
            ax = plot_gls_periodogram(
                self.gls,
                offset=0.1,
                N_peaks=3,
                relative_height=10,
                FAP_levels=[0.1, 0.01, 0.001],
                ax=ax,
            )
            # create a dummy pg class
            pg = self.raw_lc[~self.tmask].to_periodogram(method="lombscargle")
            self.Prot_ls = self.gls.best["P"]
            if math.isnan(self.gls.power.max()):
                ax.clear()
                self.pg_method = "lombscargle"
                pg = plot_periodogram(
                    self.raw_lc[~self.tmask], method="lombscargle", ax=ax
                )
                self.Prot_ls = pg.period_at_max_power.value
        else:
            self.gls = None
            pg = plot_periodogram(
                self.raw_lc[~self.tmask], method=self.pg_method, ax=ax
            )
            self.Prot_ls = pg.period_at_max_power.value

        # add gls to tls_results
        self.append_tls_results()

        # +++++++++++++++++++++ax phase-folded at Prot + sinusoidal model
        ax = axes.flatten()[2]
        # raw
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
                label="masked and folded at Prot",
                show_colorbar=False,  # colorbar_label="Time [BTJD]"
            )
        )
        _ = (
            pg.model(self.raw_lc[~self.tmask].time)
            .fold(
                self.Prot_ls,
                normalize_phase=False,
                wrap_phase=self.Prot_ls / 2,
            )
            .plot(label=f"{self.pg_method} model", color="r", lw=3, ax=ax)
        )
        ax.set_xlabel("Rotation Phase [days]")
        # ax.set_xlim(-0.5, 0.5)

        # +++++++++++++++++++++ax5: TLS periodogram
        ax = axes.flatten()[4]
        ax = plot_tls(
            self.tls_results,
            period_min=self.Porb_min,
            period_max=self.Porb_max,
            ax=ax,
        )

        # +++++++++++++++++++++ax4: flattened lc
        ax = axes.flatten()[3]
        self.flat_lc.scatter(ax=ax, label="flat")
        tmask2 = self.flat_lc.create_transit_mask(
            transit_time=self.tls_results.T0,
            period=self.tls_results.period,
            duration=self.tls_results.duration,
        )
        self.flat_lc[tmask2].scatter(ax=ax, color="r", label="transit")
        # +++++++++++++++++++++ax: tpf
        ax = axes.flatten()[5]
        if self.pipeline in [
            "cdips",
            "pathos",
            "tglc",
            "tasoc",
            "gsfc-eleanor-lite",
        ]:
            errmsg = "Pipeline to be added soon."
            raise NotImplementedError(errmsg)
        elif self.pipeline == "qlp":
            self.tpf = self.get_tpf_tesscut()
            self.sap_mask = "square"
        else:
            self.tpf = self.get_tpf(
                sector=self.sector,
                author=self.pipeline,  # exptime=self.exptime, #cadence=self.cadence
            )
            self.sap_mask = "pipeline"
        try:
            # query image to get projection
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
                # threshold_sigma=l.threshold_sigma,
                # percentile=l.percentile,
                survey=self.archival_survey,
                fov_rad=fov_rad,
                verbose=True,
                ax=ax,
            )
        except Exception as e:
            print(e)
            print("Querying archival image failed. Plotting tpf instead.")
            _ = plot_gaia_sources_on_tpf(
                tpf=self.tpf,
                target_gaiaid=self.gaiaid,
                gaia_sources=None,
                kmax=1,
                depth=1 - self.tls_results.depth,
                sap_mask=self.sap_mask,
                aper_radius=2,
                # threshold_sigma=l.threshold_sigma,
                # percentile=l.percentile,
                cmap="viridis",
                dmag_limit=8,
                verbose=True,
                ax=ax,
            )

        # +++++++++++++++++++++ax: odd-even
        ax = axes.flatten()[6]
        _ = plot_odd_even_transit(
            self.fold_lc, self.tls_results, bin_mins=10, ax=ax
        )

        # +++++++++++++++++++++ax: secondary eclipse
        ax = axes.flatten()[7]
        _ = plot_secondary_eclipse(
            self.flat_lc, self.tls_results, tmask2, bin_mins=10, ax=ax
        )

        # +++++++++++++++++++++ax: summary panel
        ax = axes.flatten()[8]
        ax.axis([0, 10, 0, 10])
        msg = self.make_summary_info()
        ax.text(-1, 11, msg, ha="left", va="top", fontsize=10, wrap=True)
        ax.axis("off")

        title = ""
        if self.toiid is not None:
            title = f"TOI {self.toiid} | "
        title += f"TIC {self.ticid} "
        title += f"| sector {self.sector} "
        title += f"| {self.pipeline.upper()} pipeline"
        # if toi_params["Comments"] or toi_params["Comments"] != "nan":
        #     comment = f"Comment: {toi_params['Comments']}"
        #     msg += "\n".join(textwrap.wrap(comment, 60))
        fig.suptitle(title, y=1.0, fontsize=20)

        if (self.outdir is not None) & (not Path(self.outdir).exists()):
            Path(self.outdir).mkdir()

        name = self.target_name.replace(" ", "")
        fp = Path(
            self.outdir,
            f"{name}_s{str(self.sector).zfill(2)}_{self.pipeline}_{self.cadence[0]}c",
        )
        fp = self.check_file_exists()
        png_file = fp.with_suffix(".png")
        if self.savefig:
            save_figure(fig, png_file, dpi=100, writepdf=False)

        if self.savetls:
            h5_file = Path(self.outdir, fp.name + "_tls").with_suffix(".h5")
            fk.save(h5_file, self.tls_results)
            print("Saved: ", h5_file)

        self.timer_end = timer()
        elapsed_time = self.timer_end - self.timer_start
        print(f"#----------Runtime: {elapsed_time:.2f} s----------#\n")
        if not self.plot:
            # del fig
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
        )
        fig = ql.plot_tql()

    except Exception:
        # Get current system exception
        ex_type, ex_value, ex_traceback = sys.exc_info()
        # Extract unformatter stack traces as tuples
        trace_back = traceback.extract_tb(ex_traceback)

        print(f"\nException type: {ex_type.__name__}")
        print(f"Exception message: {ex_value}")
        # pdb.post_mortem(ex_traceback)
        # Format stacktrace
        for trace in trace_back:
            print(f"Line : {trace[1]}")
            print(f"Func : {trace[2]}")
            # print(f"Message : {trace[3]}")
            print(f"File : {trace[0]}")
        # traceback.print_exc()
