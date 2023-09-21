import sys
import traceback
import textwrap
from pathlib import Path
from time import time as timer
import matplotlib.pyplot as pl
import numpy as np
from transitleastsquares import transitleastsquares as tls
from wotan import flatten
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clip
from astropy.wcs import WCS
import astropy.units as u
import lightkurve as lk
from astroquery.skyview import SkyView
from aesthetic.plot import set_style
from aesthetic.plot import savefig as save_figure
import flammkuchen as fk
from utils import get_tfop_info, TESS_TIME_OFFSET, TESS_pix_scale
from plot import (
    plot_gaia_sources_on_survey,
    plot_gaia_sources_on_tpf,
    plot_odd_even_transit,
    plot_secondary_eclipse,
    plot_tls,
    plot_periodogram,
)

set_style("science")


class TessQuickLook:
    def __init__(
        self,
        target_name: str,
        sector=-1,
        author: str = "SPOC",
        flux_type="pdcsap",
        cadence: str = "short",
        pg_method: str = "lombscargle",
        flatten_method: str = "biweight",
        gp_kernel: str = "matern",
        gp_kernel_size: float = 1,
        window_length: float = None,
        edge_cutoff: float = 0.1,
        sigma_clip_raw: tuple = (5, 5),
        sigma_clip_flat: tuple = None,
        ephem_mask: list = None,
        Porb_min: float = None,
        Porb_max: float = None,
        plot: bool = True,
        savefig: bool = False,
        savetls: bool = False,
        overwrite: bool = False,
        outdir: str = ".",
    ):
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
            idx = [i[:3].lower() == "toi" for i in self.names]
            toiid = int(self.names[idx][0].split("-")[-1])
        self.toiid = toiid
        self.ticid = int(self.tfop_info.get("basic_info")["tic_id"])
        self.flux_type = flux_type
        self.cadence = cadence
        self.pg_method = pg_method
        self.raw_lc = self.get_lc(
            author=author, sector=sector, cadence=cadence
        )
        self.flatten_method = flatten_method
        self.gp_kernel = (
            gp_kernel,
        )  # squared_exp, matern, periodic, periodic_auto
        self.gp_kernel_size = (gp_kernel_size,)
        self.sigma_clip_raw = (sigma_clip_raw,)
        self.sigma_clip_flat = (sigma_clip_flat,)
        self.edge_cutoff = edge_cutoff
        self.ephem_mask = ephem_mask
        _ = self.get_toi_ephem()
        if window_length is None:
            self.window_length = (
                self.tfop_dur[0] * 3 if self.tfop_dur[0] else 0.5
            )
        else:
            self.window_length = window_length
        # self.flat_lc, self.trend_lc = self.raw_lc.flatten(return_trend=True)
        self.flat_lc, self.trend_lc = self.flatten_raw_lc()
        self.tmask = self.get_transit_mask()
        self.Porb_min = 0.1 if Porb_min is None else Porb_min
        self.Porb_max = (
            (max(self.flat_lc.time.value) - min(self.flat_lc.time.value)) / 2
            if Porb_max is None
            else Porb_max
        )
        self.tls_results = tls(
            self.flat_lc.time.value,
            self.flat_lc.flux.value,
            self.flat_lc.flux_err.value,
        ).power(
            period_min=self.Porb_min,  # Roche limit default
            period_max=self.Porb_max,
        )
        self.tls_results["Porb_min"], self.tls_results["Porb_max"] = (
            self.Porb_min,
            self.Porb_max,
        )
        self.fold_lc = self.flat_lc.fold(
            period=self.tls_results.period, epoch_time=self.tls_results.T0
        )
        self.plot = plot
        self.savefig = savefig
        self.savetls = savetls
        self.overwrite = overwrite
        self.outdir = outdir
        _ = self.plot_tql()

    def __repr__(self):
        """Override to print a readable string representation of class"""

        included_args = [
            # ===target attributes===
            "name",
            "toiid",
            "ctoiid",
            "ticid",
            "epicid",
            "gaiaDR2id",
            "ra_deg",
            "dec_deg",
            "target_coord",
            "search_radius",
            "mission",
            "campaign",
            # "all_sectors",
            # "all_campaigns",
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

        if self.target_name[:3].lower() == "toi":
            query_name = f"TIC{self.ticid}"
        else:
            query_name = self.target_name.replace(" ", "")

        # get all info
        search_result = lk.search_lightcurve(query_name)
        all_sectors = sorted(
            set([int(i[-2:]) for i in search_result.table["mission"].tolist()])
        )
        if kwargs.get("sector") is None:
            print(f"Available sectors: {all_sectors}")
        else:
            err_msg = f"sector={kwargs.get('sector')} not in {all_sectors}"
            assert kwargs.get("sector") in all_sectors, err_msg

        if kwargs.get("author") is None:
            kwargs["author"] = "SPOC"
        else:
            all_authors = set(search_result.table["provenance_name"].tolist())
            err_msg = f"author={kwargs.get('author')} not in {all_authors}"
            assert kwargs.get("author").upper() in all_authors, err_msg
        self.pipeline = kwargs["author"]
        self.all_pipelines = all_authors

        # get specific lc in cache
        search_result = lk.search_lightcurve(query_name, **kwargs)

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
            print(
                f"Downloaded {kwargs.get('author')} lc in sector {lc.sector}."
            )
            self.sector = lc.sector
        if lc.meta["AUTHOR"].lower() == "spoc":
            lc = lc.select_flux(self.flux_type + "_flux")
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

        if self.target_name[:3].lower() == "toi":
            query_name = f"TIC{self.ticid}"
        else:
            query_name = self.target_name.replace(" ", "")
        search_result = lk.search_targetpixelfile(query_name, **kwargs)
        idx = sector_orig if sector_orig == -1 else 0
        tpf = search_result[idx].download()
        print(f"Downloaded {kwargs.get('author')} tpf in sector {tpf.sector}.")
        return tpf

    def get_toi_ephem(self, idx=1, params=["epoch", "per", "dur"]) -> list:
        print(f"Querying ephemeris for {self.target_name}:")
        planet_params = self.tfop_info.get("planet_parameters")[idx]
        vals = []
        for p in params:
            val = planet_params.get(p)
            val = float(val) if val else 0.1
            err = planet_params.get(p + "_e")
            err = float(err) if err else 0.1
            print(f"{p}: {val}, {err}")
            vals.append((val, err))
        self.tfop_epoch = np.array(vals[0]) - 2457000
        self.tfop_period = np.array(vals[1])
        self.tfop_dur = np.array(vals[2]) / 24
        return vals

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
        if self.sigma_clip_flat[0] is not None:
            print(
                f"Applying sigma clip on flattened lc with \
                 (lower,upper)=({self.sigma_clip_flat})"
            )
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
        Rstar = float(star_params["srad"])
        Rstar_err = float(star_params["srad_e"])
        Mstar = float(star_params["mass"])
        Mstar_err = float(star_params["mass_e"])
        Teff = float(star_params["teff"])
        Teff_err = float(star_params["teff_e"])
        logg = float(star_params["logg"])
        logg_err = float(star_params["logg_e"])
        # feh = star_params["met"]
        # feh_err = star_params["met_e"]
        dist = float(star_params["dist"])
        dist_err = float(star_params["dist_e"])
        meta = self.raw_lc.meta

        Rp = self.tls_results["rp_rs"] * Rstar * u.Rsun.to(u.Rearth)
        Rp_true = Rp * np.sqrt(1 + meta["CROWDSAP"])
        msg = "\nCandidate Properties\n"
        msg += "-" * 30 + "\n"
        text = f"SDE={self.tls_results.SDE:.4f} (sector={self.sector} \
                in {self.all_sectors}, {self.pipeline.upper()} pipeline)"
        msg += "\n".join(textwrap.wrap(text, 60))
        msg += (
            f"\nPeriod={self.tls_results.period:.4f}\
                +/-{self.tls_results.period_uncertainty:.4f} d"
            + " " * 5
        )
        msg += f"Duration={self.tls_results.duration*24:.2f} hr" + "\n"
        msg += f"T0={self.tls_results.T0+TESS_TIME_OFFSET:.4f} BJD" + " " * 11
        msg += f"Depth={(1-self.tls_results.depth)*100:.2f}%\n"

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
            msg += f"Rp={Rp:.2f} " + r"R$_{\oplus}$" + "(diluted)" + " " * 10
            msg += f"Rp={Rp_true:.2f} " + r"R$_{\oplus}$" + "(undiluted)\n"
        msg += (
            f"Odd-Even mismatch={self.tls_results.odd_even_mismatch:.2f}"
            + r"$\sigma$"
        )
        msg += "\n" * 2
        msg += "Stellar Properties\n"
        msg += "-" * 30 + "\n"
        msg += f"TIC ID={self.ticid}" + " " * 5
        msg += f"Tmag={meta['TESSMAG']:.2f}\n"
        # msg += f"Gaia DR2 ID={self.gaiaid}\n"
        msg += f"Distance={dist:.1f}+/-{dist_err:.1f}pc\n"
        # msg += f"GOF_AL={astrometric_gof_al:.2f} (hints binarity if >20)\n"
        # D = gp.astrometric_excess_noise_sig
        # msg += f"astro. excess noise sig={D:.2f} (hints binarity if >5)\n"
        msg += (
            f"Rstar={Rstar:.2f}+/-{Rstar_err:.2f} " + r"R$_{\odot}$" + " " * 5
        )
        msg += f"Teff={Teff}+/-{Teff_err} K" + "\n"
        msg += (
            f"Mstar={Mstar:.2f}+/-{Mstar_err:.2f} " + r"M$_{\odot}$" + " " * 5
        )
        msg += f"logg={logg:.2f}+/-{logg_err:.2f} cgs\n"
        # msg += f"met={feh:.2f}+/-{feh_err:.2f} dex " + " " * 6
        return msg

    def plot_tql(self, **kwargs: dict) -> pl.Figure:
        if kwargs.get("savefig"):
            self.savefig = kwargs.get("savefig")
        if kwargs.get("overwrite"):
            self.overwrite = kwargs.get("overwite")
        if kwargs.get("plot"):
            self.plot = kwargs.get("plot")

        start = timer()
        fig, axes = pl.subplots(3, 3, figsize=(16, 12), tight_layout=True)

        # +++++++++++++++++++++ax: Raw + trend
        ax = axes.flatten()[0]
        self.raw_lc.scatter(ax=ax, label="raw")
        label = f"baseline trend\nmethod={self.flatten_method}"
        label += f"(window_size={self.window_length:.2f})"
        self.trend_lc.plot(ax=ax, color="r", lw=2, label=label)

        # +++++++++++++++++++++ax2 Lomb-scargle periodogram
        ax = axes.flatten()[1]

        best_period, ls_model = plot_periodogram(
            self.raw_lc[~self.tmask], method=self.pg_method, ax=ax
        )
        # +++++++++++++++++++++ax phase-folded at Prot + sinusoidal model
        ax = axes.flatten()[2]
        # raw
        _ = (
            self.raw_lc[~self.tmask]
            .fold(best_period)
            .scatter(
                ax=ax,
                c=self.raw_lc.time.value[~self.tmask],
                cmap=pl.get_cmap("Blues_r"),
                label="masked and folded at Prot",
                show_colorbar=False,  # colorbar_label="Time [BTJD]"
            )
        )
        _ = (
            ls_model(self.raw_lc.time)
            .fold(best_period)
            .plot(label=f"{self.pg_method} model", color="r", lw=3, ax=ax)
        )

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
            duration=self.tls_results.duration * 24,
        )
        self.flat_lc[tmask2].scatter(ax=ax, color="r", label="transit")
        # +++++++++++++++++++++ax: tpf
        ax = axes.flatten()[5].remove()
        self.tpf = self.get_tpf(
            sector=self.sector, author=self.pipeline, cadence=self.cadence
        )
        try:
            survey = "DSS2 Red"
            # query image to get projection
            ny, nx = self.tpf.flux.shape[1:]
            diag = np.sqrt(nx**2 + ny**2)
            fov_rad = (0.4 * diag * TESS_pix_scale).to(u.arcmin)
            position = self.target_coord.icrs.to_string()
            results = SkyView.get_images(
                position=position,
                coordinates="icrs",
                survey=survey,
                radius=fov_rad,
                grid=True,
            )
            if len(results) > 0:
                hdu = results[0][0]
            else:
                errmsg = (
                    "SkyView returned empty result. Try a different survey."
                )
                raise ValueError(errmsg)
            ax = fig.add_subplot(3, 3, 6, projection=WCS(hdu.header))
            _ = plot_gaia_sources_on_survey(
                tpf=self.tpf,
                target_gaiaid=self.gaiaid,
                gaia_sources=None,
                kmax=1,
                depth=1 - self.tls_results.depth,
                # sap_mask='pipeline',
                # aper_radius=l.aper_radius,
                # threshold_sigma=l.threshold_sigma,
                # percentile=l.percentile,
                survey=survey,
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
                # gaia_sources=None,
                kmax=1,
                depth=1 - self.tls_results.depth,
                # sap_mask=l.sap_mask,
                # aper_radius=l.aper_radius,
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

        if self.toiid is None:
            title = f"TIC {self.ticid} (sector {self.sector})"
        else:
            title = (
                f"TOI {self.toiid} | TIC {self.ticid} (sector {self.sector})"
            )
            # if toi_params["Comments"] or toi_params["Comments"] != "nan":
            #     comment = f"Comment: {toi_params['Comments']}"
            #     msg += "\n".join(textwrap.wrap(comment, 60))
        fig.suptitle(title, y=1.0, fontsize=20)
        end = timer()

        if (self.outdir is not None) & (not Path(self.outdir).exists()):
            Path(self.outdir).mkdir()

        lctype = self.flux_type if self.pipeline == "spoc" else self.pipeline
        name = self.target_name.replace(" ", "")
        fp = Path(
            self.outdir,
            f"{name}_s{str(self.sector).zfill(2)}_{lctype}_{self.cadence[0]}c",
        )
        png_file = fp.with_suffix(".png")
        if png_file.exists() and not self.overwrite:
            raise FileExistsError(
                f"{png_file} already exists! Set overwrite=True."
            )

        if self.savefig:
            save_figure(fig, png_file, dpi=100, writepdf=False)

        if self.savetls:
            h5_file = Path(self.outdir, fp.name + "_tls").with_suffix(".h5")
            # self.tls_results["gaiaid"] = self.gaiaid
            self.tls_results["ticid"] = self.ticid
            fk.save(h5_file, self.tls_results)
            print("Saved: ", h5_file)
        print(f"#----------Runtime: {end-start:.2f} s----------#\n")
        if not self.plot:
            # del fig
            pl.clf()
        return fig


if __name__ == "__main__":
    try:
        tql = TessQuickLook(
            "TOI-837",
            sector=-1,
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
