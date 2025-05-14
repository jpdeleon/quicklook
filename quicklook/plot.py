import importlib.resources as pkg_resources
import matplotlib.pyplot as pl
import numpy as np
from matplotlib.patches import Circle
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.visualization import ZScaleInterval
import astropy.units as u

# from astroplan.plots import plot_finder_image
from astroquery.mast import Catalogs
from scipy.ndimage import zoom
import pandas as pd
from quicklook.utils import (
    compute_secthresh,
    is_point_inside_mask,
    PadWithZeros,
    parse_aperture_mask,
    TESS_pix_scale,
)
from quicklook.measure import find_contours

__all__ = [
    "use_style",
    "plot_odd_even_transit",
    "plot_secondary_eclipse",
    "plot_periodogram",
    "plot_gaia_sources_on_tpf",
    "plot_gaia_sources_on_survey",
    "plot_archival_images",
    "plot_tls",
]

# http://gsss.stsci.edu/SkySurveys/Surveys.htm
dss_description = {
    "dss1": "POSS1 Red in the north; POSS2/UKSTU Blue in the south",
    "poss2ukstu_red": "POSS2/UKSTU Red",
    "poss2ukstu_ir": "POSS2/UKSTU Infrared",
    "poss2ukstu_blue": "POSS2/UKSTU Blue",
    "poss1_blue": "POSS1 Blue",
    "poss1_red": "POSS1 Red",
    "all": "best among all plates",
    "quickv": "Quick-V Survey",
    "phase2_gsc2": "HST Phase 2 Target Positioning (GSC 2)",
    "phase2_gsc1": "HST Phase 2 Target Positioning (GSC 1)",
}


def use_style(name="science"):
    with pkg_resources.path(__package__, f"{name}.mplstyle") as style_path:
        pl.style.use(str(style_path))


def plot_odd_even_transit(fold_lc, tls_results, bin_mins=10, ax=None):
    if ax is None:
        _, ax = pl.subplots()
    yline = tls_results.depth
    fold_lc.scatter(ax=ax, c="k", alpha=0.5, label="_nolegend_", zorder=1)
    fold_lc[fold_lc.even_mask].bin(time_bin_size=bin_mins * u.minute).errorbar(
        label="even transit",
        c="#1f77b4",
        marker="o",
        lw=2,
        markersize=8,
        ax=ax,
        zorder=2,
    )
    fold_lc[fold_lc.odd_mask].bin(time_bin_size=bin_mins * u.minute).errorbar(
        label="odd transit",
        c="#d62728",
        marker="o",
        lw=2,
        markersize=8,
        ax=ax,
        zorder=2,
    )
    ax.plot(
        (tls_results.model_folded_phase - 0.5) * tls_results.period,
        tls_results.model_folded_model,
        "-k",
        lw=3,
        zorder=3,
        label="TLS model",
    )
    ax.axhline(yline, 0, 1, lw=2, ls="--", c="k")
    # transit duration in phase
    t14 = tls_results.duration
    ax.axvline(-t14 / 2, 0, 1, label="__nolegend__", c="k", ls="--")
    ax.axvline(t14 / 2, 0, 1, label="__nolegend__", c="k", ls="--")
    ax.set_xlabel("Orbital Phase")
    ax.set_xlim(-t14 * 2, t14 * 2)
    # y1, y2 = ax.get_ylim()
    ax.legend()
    return ax


def plot_secondary_eclipse(flat_lc, tls_results, tmask, bin_mins=10, ax=None):
    if ax is None:
        _, ax = pl.subplots()
    # mask transit and shift phase
    fold_lc2 = flat_lc[~tmask].fold(
        period=tls_results.period,
        epoch_time=tls_results.T0 + tls_results.period / 2,
        # normalize_phase=False,
        # wrap_phase=tls_results.period
    )
    half_phase = 0.5  # tls_results.period/2
    fold_lc2.time = fold_lc2.time + half_phase * u.day
    fold_lc2.scatter(ax=ax, c="k", alpha=0.5, label="_nolegend_", zorder=1)
    yline = tls_results.depth
    ax.axhline(
        yline,
        0,
        1,
        lw=2,
        label=f"TLS depth={(1-yline)*1e3:.2f} ppt",
        c="k",
        ls="--",
    )
    t14 = tls_results.duration
    ax.axvline(
        half_phase - t14 / 2, 0, 1, label="__nolegend__", c="k", ls="--"
    )
    ax.axvline(
        half_phase + t14 / 2, 0, 1, label="__nolegend__", c="k", ls="--"
    )
    try:
        secthresh = compute_secthresh(fold_lc2, t14)
    except Exception as e:
        print(e)
        secthresh = np.nan
    fold_lc2.scatter(ax=ax, c="k", alpha=0.5, label="_nolegend_", zorder=1)
    try:
        fold_lc2.bin(time_bin_size=bin_mins * u.minute).errorbar(
            ax=ax,
            marker="o",
            markersize=8,
            lw=2,
            label=f"sec_eclipse_thresh={secthresh*1e3:.2f} ppt",
            zorder=2,
        )
    except Exception as e:
        print(e)
    ax.set_xlabel("Orbital Phase")
    ax.set_xlim(half_phase - t14 * 2, half_phase + t14 * 2)
    ax.legend()
    return ax


def plot_periodogram(lc, method="lombscargle", verbose=True, ax=None) -> tuple:
    if ax is None:
        _, ax = pl.subplots()
    baseline = int(lc.time.value[-1] - lc.time.value[0])
    Prot_max = baseline / 2.0
    pg = lc.to_periodogram(
        method=method,
        # minimum_period=0.5,
        maximum_period=Prot_max,
        # minimum_frequency = 2.0,
        # maximum_frequency = 1/Prot_max
    )
    best_period = pg.period_at_max_power.value
    pg.plot(ax=ax, view="period", lw=2, color="k", label="__nolegend__")
    ax.axvline(
        best_period,
        0,
        1,
        color="r",
        ls="--",
        lw=2,
        label=f"peak={best_period:.2f}",
    )
    ax.set_xscale("log")
    ax.set_ylabel("Lomb-Scargle Power")
    ax.set_xlabel("Rotation period [days]")
    ax.legend(title="Rotation period [d]")
    xmin, _ = ax.get_xlim()
    ax.set_xlim(xmin, Prot_max)
    if verbose:
        print(pg.show_properties())
    return pg


def plot_gls_periodogram(
    gls,
    offset=0.1,
    N_peaks=3,
    relative_height=10,
    FAP_levels=[0.1, 0.01, 0.001],
    ax=None,
    verbose=True,
) -> tuple:
    """
    Based on:
    https://github.com/SLSkrzypinski/TESS_diagnosis
    """
    from scipy.signal import find_peaks

    linestyles = [":", "dotted", "solid"]
    if ax is None:
        _, ax = pl.subplots()
    x = 1 / gls.freq
    # idx = np.argsort(x)
    # x = x[idx]
    # y = gls.power[idx]
    y = gls.power
    power_levels = [[gls.powerLevel(i)] * len(x) for i in FAP_levels]
    best_period = gls.best["P"]
    # Find peaks
    max_power = y.max()
    peaks = find_peaks(y, height=max_power / relative_height)
    # print(peaks)
    peak_pos = peaks[0]
    peak_pos = peak_pos[
        (x[peak_pos] < best_period - offset)
        | (x[peak_pos] > best_period + offset)
    ]
    peak_pos = peak_pos[
        (x[peak_pos] < best_period / 2 - offset)
        | (x[peak_pos] > best_period / 2 + offset)
    ]
    peak_pos = peak_pos[
        (x[peak_pos] < 2 * best_period - offset)
        | (x[peak_pos] > 2 * best_period + offset)
    ]
    while len(peak_pos) > N_peaks:
        peak_pos = np.delete(peak_pos, peak_pos.argmin())
    peaks = x[peak_pos]
    heights = y[peak_pos]

    ax.plot(x, y, c="C0", linewidth=0.8)
    # mark best period and multiples
    ax.axvline(
        best_period,
        0,
        1,
        c="k",
        lw=2,
        alpha=0.5,
        label=f"best={best_period:.3}",
    )
    # mark best period and multiples
    if best_period * 2 <= max(x):
        ax.axvline(
            x=best_period * 2, color="k", ls="--", linewidth=2, alpha=0.5
        )
    ax.axvline(
        x=best_period / 2,
        color="k",
        ls="--",
        linewidth=2,
        alpha=0.5,
        label="best alias",
    )
    for i in range(len(peaks)):
        ax.scatter(
            peaks[i],
            heights[i],
            c="k",
            s=20,
            label=f"P$_{i+1}$={peaks[i]:.3f}",
        )
    # mark FAP levels
    for i in range(len(FAP_levels)):
        label = (
            f"FAP={max(FAP_levels)}" if i == np.argmax(FAP_levels) else None
        )
        ax.plot(
            x,
            power_levels[i],
            linestyle=linestyles[i],
            linewidth=0.8,
            c="red",
            label=label,
        )
    ax.minorticks_on()
    ax.set_ylabel("Gen. Lomb-Scargle Power")
    ax.set_xlabel("Rotation period [days]")
    ax.legend(title="Prot peaks [d]")
    return ax


def plot_tpf(tpf, aperture_mask="pipeline", ax=None) -> pl.axis:
    if ax is None:
        _, ax = pl.subplots()
    tpf.plot(
        ax=ax,
        aperture_mask=aperture_mask,
        mask_color="red",
        show_colorbar=False,
    )
    ax.set_title("")
    return ax


def plot_dss_image(
    hdu, cmap="gray", contrast=0.5, coord_format="dd:mm:ss", ax=None
):
    """
    Plot output of get_dss_data:
    hdu = get_dss_data(ra, dec)
    """
    data, header = hdu.data, hdu.header
    interval = ZScaleInterval(contrast=contrast)
    zmin, zmax = interval.get_limits(data)

    if ax is None:
        fig = pl.figure(constrained_layout=True)
        ax = fig.add_subplot(projection=WCS(header))
    ax.imshow(data, vmin=zmin, vmax=zmax, cmap=cmap)
    ax.set_xlabel("RA")
    ax.set_ylabel("DEC", y=0.9)
    title = f"{header['SURVEY']} ({header['FILTER']})\n"
    title += f"{header['DATE-OBS'][:10]}"
    ax.set_title(title)
    # set RA from hourangle to degree
    if hasattr(ax, "coords"):
        ax.coords[1].set_major_formatter(coord_format)
        ax.coords[0].set_major_formatter(coord_format)
    return ax


def plot_aperture_outline(
    img,
    mask,
    ax=None,
    imgwcs=None,
    cmap="viridis",
    color_aper="C6",
    figsize=None,
):
    """
    see https://github.com/rodluger/everest/blob/master/everest/standalone.py
    """
    interval = ZScaleInterval(contrast=0.5)
    ny, nx = mask.shape
    contour = np.zeros((ny, nx))
    contour[np.where(mask)] = 1
    contour = np.lib.pad(contour, 1, PadWithZeros)
    highres = zoom(contour, 100, order=0, mode="nearest")
    extent = np.array([-1, nx, -1, ny])

    if ax is None:
        fig, ax = pl.subplots(
            subplot_kw={"projection": imgwcs}, figsize=figsize
        )
        ax.set_xlabel("RA")
        ax.set_ylabel("Dec")
    _ = ax.contour(
        highres,
        levels=[0.5],
        linewidths=[3],
        extent=extent,
        origin="lower",
        colors=color_aper,
    )
    zmin, zmax = interval.get_limits(img)
    ax.matshow(
        img, origin="lower", cmap=cmap, vmin=zmin, vmax=zmax, extent=extent
    )
    # verts = cs.allsegs[0][0]
    return ax


def plot_gaia_sources_on_tpf(
    tpf,
    target_gaiaid,
    gaia_sources=None,
    sap_mask="pipeline",
    depth=None,
    kmax=1,
    dmag_limit=8,
    fov_rad=None,
    cmap="viridis",
    figsize=None,
    ax=None,
    invert_xaxis=False,
    invert_yaxis=False,
    pix_scale=TESS_pix_scale,
    verbose=True,
    **mask_kwargs,
):
    """
    plot Gaia sources brighter than dmag_limit; only annotated with starids
    are those that are bright enough to cause reproduce the transit depth;
    starids are in increasing separation

    dmag_limit : float
        maximum delta mag to consider; computed based on depth if None

    TODO: correct for proper motion difference between
    survey image and Gaia DR2 positions
    """
    if verbose:
        print("Plotting nearby Gaia sources on tpf.")
    assert target_gaiaid is not None
    img = np.nanmedian(tpf.flux, axis=0)
    # make aperture mask
    mask = parse_aperture_mask(tpf, sap_mask=sap_mask, **mask_kwargs)

    ax = plot_aperture_outline(
        img, mask=mask, imgwcs=tpf.wcs, figsize=figsize, cmap=cmap, ax=ax
    )
    if fov_rad is None:
        nx, ny = tpf.shape[1:]
        diag = np.sqrt(nx**2 + ny**2)
        fov_rad = (0.4 * diag * pix_scale).to(u.arcmin).round(0).value

    if gaia_sources is None:
        if verbose:
            print("Querying Gaia sources around the target.")
        target_coord = SkyCoord(
            ra=tpf.get_header()["RA_OBJ"],
            dec=tpf.get_header()["DEC_OBJ"],
            unit="deg",
        )
        gaia_sources = Catalogs.query_region(
            target_coord,
            radius=fov_rad,
            catalog="gaiadr3",  # version=3
        ).to_pandas()
    assert len(gaia_sources) > 1, "gaia_sources contains single entry"
    # find sources within mask
    # target is assumed to be the first row
    idx = gaia_sources["source_id"].astype(int).isin([target_gaiaid])
    target_gmag = gaia_sources.loc[idx, "phot_g_mean_mag"].values[0]
    # sources_inside_aperture = []
    if depth is not None:
        # compute delta mag limit given transit depth
        dmag_limit = (
            np.log10(kmax / depth - 1) if dmag_limit is None else dmag_limit
        )

        # get min_gmag inside mask
        ra, dec = gaia_sources[["ra", "dec"]].values.T
        pix_coords = tpf.wcs.all_world2pix(np.c_[ra, dec], 0)
        contour_points = find_contours(mask, level=0.1)[0]
        isinside = [
            is_point_inside_mask(contour_points, pix) for pix in pix_coords
        ]
        # sources_inside_aperture.append(isinside)
        min_gmag = gaia_sources.loc[isinside, "phot_g_mean_mag"].min()
        if (target_gmag - min_gmag) != 0:
            print(
                f"target Gmag={target_gmag:.2f} is not the brightest \
                within aperture (Gmag={min_gmag:.2f})"
            )
    else:
        min_gmag = gaia_sources.phot_g_mean_mag.min()  # brightest
        dmag_limit = (
            gaia_sources.phot_g_mean_mag.max()
            if dmag_limit is None
            else dmag_limit
        )
    base_ms = 128.0  # base marker size
    starid = 1
    # if very crowded, plot only top N
    gmags = gaia_sources.phot_g_mean_mag
    dmags = gmags - target_gmag
    rank = np.argsort(dmags.values)
    for _, row in gaia_sources.iterrows():
        # FIXME: why some indexes are missing?
        ra, dec, gmag, id = row[["ra", "dec", "phot_g_mean_mag", "source_id"]]
        dmag = gmag - target_gmag
        pix = tpf.wcs.all_world2pix(np.c_[ra, dec], 0)[0]
        contour_points = find_contours(mask, level=0.1)[0]

        color, alpha = "red", 1.0
        # change marker color and transparency based on location & dmag
        if is_point_inside_mask(contour_points, pix):
            if int(id) == int(target_gaiaid):
                # plot x on target
                ax.plot(
                    pix[1],
                    pix[0],
                    marker="x",
                    ms=base_ms / 16,
                    c="k",
                    zorder=3,
                )
            if depth is not None:
                # compute flux ratio with respect to brightest star
                gamma = 1 + 10 ** (0.4 * (min_gmag - gmag))
                if depth > kmax / gamma:
                    # orange if flux is insignificant
                    color = "C1"
        else:
            # outside aperture
            color, alpha = "C1", 0.5

        ax.scatter(
            pix[1],
            pix[0],
            s=base_ms / 2**dmag,  # fainter -> smaller
            c=color,
            alpha=alpha,
            zorder=2,
            edgecolor=None,
        )
        # choose which star to annotate
        if len(gmags) < 20:
            # sparse: annotate all
            ax.text(pix[1], pix[0], str(starid), color="white", zorder=100)
        elif len(gmags) > 50:
            # crowded: annotate only 15 smallest dmag ones
            if rank[starid - 1] < 15:
                ax.text(pix[1], pix[0], str(starid), color="white", zorder=100)
            elif (color == "red") & (dmag < dmag_limit):
                # plot if within aperture and significant source of dilution
                ax.text(pix[1], pix[0], str(starid), color="white", zorder=100)
        elif color == "red":
            # neither sparse nor crowded
            # annotate if inside aperture
            ax.text(pix[1], pix[0], str(starid), color="white", zorder=100)
        starid += 1
    # Make legend with 4 sizes representative of delta mags
    dmags = dmags[dmags < dmag_limit]
    _, dmags = pd.cut(dmags, 3, retbins=True)
    for dmag in dmags:
        size = base_ms / 2**dmag
        # -1, -1 is outside the fov
        # dmag = 0 if float(dmag)==0 else 0
        ax.scatter(
            -1,
            -1,
            s=size,
            c="red",
            alpha=0.6,
            edgecolor=None,
            zorder=10,
            clip_on=True,
            label=r"$\Delta m= $" + f"{dmag:.1f}",
        )
    ax.legend(fancybox=True, framealpha=0.5)
    # set img limits
    xdeg = (nx * pix_scale).to(u.arcmin)
    ydeg = (ny * pix_scale).to(u.arcmin)
    # orient such that north is up; east is left
    if invert_yaxis:
        # ax.invert_yaxis()  # increasing upward
        raise NotImplementedError()
    if invert_xaxis:
        # ax.invert_xaxis() #decresing rightward
        raise NotImplementedError()
    if hasattr(ax, "coords"):
        ax.coords[0].set_major_formatter("dd:mm")
        ax.coords[1].set_major_formatter("dd:mm")
    pl.setp(
        ax, xlim=(0, nx), ylim=(0, ny), xlabel=f"({xdeg:.2f} x {ydeg:.2f})"
    )
    return ax


def plot_gaia_sources_on_survey(
    tpf,
    target_gaiaid,
    hdu=None,
    gaia_sources=None,
    fov_rad=None,
    depth=0.0,
    kmax=1.0,
    sap_mask="pipeline",
    survey="dss1",
    ax=None,
    color_aper="C0",  # pink
    figsize=None,
    invert_xaxis=False,
    invert_yaxis=False,
    pix_scale=TESS_pix_scale,
    verbose=True,
    **mask_kwargs,
):
    """Plot (superpose) Gaia sources on archival image

    Parameters
    ----------
    target_coord : astropy.coordinates
        target coordinate
    gaia_sources : pd.DataFrame
        gaia sources table
    fov_rad : astropy.unit
        FOV radius
    survey : str
        image survey
    verbose : bool
        print texts
    ax : axis
        subplot axis
    color_aper : str
        aperture outline color (default=C6)
    kwargs : dict
        keyword arguments for aper_radius, percentile
    Returns
    -------
    ax : axis
        subplot axis

    TODO: correct for proper motion difference between
    survey image and Gaia DR2 positions
    """
    errmsg = f"{survey} not in {list(dss_description.keys())}"
    assert survey in list(dss_description.keys()), errmsg
    if verbose:
        print("Plotting nearby Gaia sources on survey image.")
    assert target_gaiaid is not None
    ny, nx = tpf.flux.shape[1:]
    if fov_rad is None:
        diag = np.sqrt(nx**2 + ny**2)
        fov_rad = (0.4 * diag * pix_scale).to(u.arcmin).round(2)
    target_coord = SkyCoord(ra=tpf.ra * u.deg, dec=tpf.dec * u.deg)
    if gaia_sources is None:
        if verbose:
            print("Querying Gaia sources around the target.")
        gaia_sources = Catalogs.query_region(
            target_coord, radius=fov_rad, catalog="Gaia", version=2
        ).to_pandas()
    assert len(gaia_sources) > 1, "gaia_sources contains single entry"
    # make aperture mask
    mask = parse_aperture_mask(tpf, sap_mask=sap_mask, **mask_kwargs)
    maskhdr = tpf.hdu[2].header
    # make aperture mask outline
    contour = np.zeros((ny, nx))
    contour[np.where(mask)] = 1
    contour = np.lib.pad(contour, 1, PadWithZeros)
    highres = zoom(contour, 100, order=0, mode="nearest")
    extent = np.array([-1, nx, -1, ny])

    if verbose:
        print(
            f"Querying {survey} ({fov_rad:.2f} x {fov_rad:.2f}) archival image..."
        )
    # -----------create figure---------------#
    if (ax is None) or (hdu is None):
        # get img hdu for subplot projection
        try:
            hdu = get_dss_data(
                ra=target_coord.ra.deg,
                dec=target_coord.dec.deg,
                survey=survey,
                width=fov_rad.value,
                height=fov_rad.value,
            )
        except Exception:
            errmsg = "survey image not available"
            raise FileNotFoundError(errmsg)
        fig = pl.figure(figsize=figsize)
        # define scaling in projection
        ax = fig.add_subplot(111, projection=WCS(hdu.header))
    # plot survey img
    ax.imshow(hdu.data, cmap="Greys", origin="lower")
    ax.set(xlabel="RA", ylabel="DEC")
    imgwcs = WCS(hdu.header)
    mx, my = hdu.data.shape
    # plot mask
    _ = ax.contour(
        highres,
        levels=[0.5],
        extent=extent,
        origin="lower",
        linewidths=[3],
        colors=color_aper,
        transform=ax.get_transform(WCS(maskhdr)),
    )
    idx = gaia_sources["source_id"].astype(int).isin([target_gaiaid])
    target_gmag = gaia_sources.loc[idx, "phot_g_mean_mag"].values[0]

    for _, row in gaia_sources.iterrows():
        marker, s = "o", 100
        r, d, mag, id = row[["ra", "dec", "phot_g_mean_mag", "source_id"]]
        pix = imgwcs.all_world2pix(np.c_[r, d], 1)[0]
        if int(id) != int(target_gaiaid):
            gamma = 1 + 10 ** (0.4 * (mag - target_gmag))
            if depth > kmax / gamma:
                # too deep to have originated from secondary star
                edgecolor = "C1"
                alpha = 1  # 0.5
            else:
                # possible NEBs
                edgecolor = "C3"
                alpha = 1
        else:
            s = 200
            edgecolor = "C2"
            marker = "s"
            alpha = 1
        ax.scatter(
            pix[0],
            pix[1],
            marker=marker,
            s=s,
            edgecolor=edgecolor,
            alpha=alpha,
            facecolor="none",
        )
    # orient such that north is up; left is east
    if invert_yaxis:
        # ax.invert_yaxis()
        raise NotImplementedError()
    if invert_xaxis:
        # ax.invert_xaxis()
        raise NotImplementedError()
    if hasattr(ax, "coords"):
        ax.coords[0].set_major_formatter("dd:mm")
        ax.coords[1].set_major_formatter("dd:mm")
    # set img limits
    pl.setp(
        ax,
        xlim=(0, mx),
        ylim=(0, my),
    )
    ax.set_title(
        f"{survey.upper()} survey (FOV={fov_rad.value:.2f}' x {fov_rad.value:.2f}')",
        y=0.99,
    )
    return ax


def plot_tls(
    tls_results: dict,
    period_min: float = None,
    period_max: float = None,
    ax=None,
) -> pl.axis:
    if ax is None:
        _, ax = pl.subplots()

    label = f"best={tls_results.period:.3f}"
    ax.axvline(tls_results.period, alpha=0.4, lw=3, label=label)
    ax.set_xlim(np.min(tls_results.periods), np.max(tls_results.periods))

    for i in range(2, 10):
        higher_harmonics = i * tls_results.period
        if period_min <= higher_harmonics <= period_max:
            ax.axvline(higher_harmonics, alpha=0.4, lw=1, linestyle="dashed")
        lower_harmonics = tls_results.period / i
        if period_min <= lower_harmonics <= period_max:
            ax.axvline(lower_harmonics, alpha=0.4, lw=1, linestyle="dashed")
    ax.set_ylabel("Transit Least Squares SDE")
    ax.set_xlabel("Orbital Period [days]")
    ax.plot(tls_results.periods, tls_results.power, color="black", lw=0.5)
    ax.set_xlim(tls_results["Porb_min"], tls_results["Porb_max"])
    # do not show negative SDE
    y1, y2 = ax.get_ylim()
    y1 = 0 if y1 < 0 else y1
    ax.set_ylim(y1, y2)
    ax.legend(title="Porb peaks [d]")
    return ax


def get_dss_data(
    ra,
    dec,
    survey="poss2ukstu_red",
    plot=False,
    height=1,
    width=1,
    epoch="J2000",
):
    """
    Digitized Sky Survey (DSS)
    http://archive.stsci.edu/cgi-bin/dss_form
    Parameters
    ----------
    survey : str
        (default=poss2ukstu_red) see `dss_description`
    height, width : float
        image cutout height and width [arcmin]
    epoch : str
        default=J2000
    Returns
    -------
    hdu
    """
    survey_list = list(dss_description.keys())
    if survey not in survey_list:
        raise ValueError(f"{survey} not in:\n{survey_list}")
    base_url = "http://archive.stsci.edu/cgi-bin/dss_search?v="
    url = f"{base_url}{survey}&r={ra}&d={dec}&e={epoch}&h={height}&w={width}"
    url += "&f=fits&c=none&s=on&fov=NONE&v3"
    try:
        hdulist = fits.open(url)
        # hdulist.info()

        hdu = hdulist[0]
        # data = hdu.data
        # header = hdu.header
        if plot:
            _ = plot_dss_image(hdu)
        return hdu
    except Exception as e:
        if isinstance(e, OSError):
            print(f"Error: {e}\nsurvey={survey} image is likely unavailable.")
        else:
            raise Exception(f"Error: {e}")


def plot_archival_images(
    ra,
    dec,
    survey1="dss1",
    survey2="ps1",  # "poss2ukstu_red",
    filter="i",
    fp1=None,
    fp2=None,
    height=1,
    width=1,
    cmap="gray",
    reticle=True,
    grid=True,
    color="red",
    contrast=0.5,
    fontsize=14,
    coord_format="dd:mm:ss",
    return_baseline=False,
):
    """
    Plot two archival images
    See e.g.
    https://s3.amazonaws.com/aasie/images/1538-3881/159/3/100/ajab5f15f2_hr.jpg
    Uses reproject to have identical fov:
    https://reproject.readthedocs.io/en/stable/

    Parameters
    ----------
    ra, dec : float
        target coordinates in degrees
    survey1, survey2 : str
        survey from which the images will come from
    fp1, fp2 : path
        filepaths if the images were downloaded locally
    height, width
        fov of view in arcmin (default=1')
    filter : str
        (g,r,i,z,y) filter if survey = PS1
    cmap : str
        colormap (default='gray')
    reticle : bool
        plot circle to mark the original position of target in survey1
    color : str
        default='red'
    contrast : float
        ZScale contrast
    Notes:
    ------
    Account for space motion:
    https://docs.astropy.org/en/stable/coordinates/apply_space_motion.html

    The position offset can be computed as:
    ```
    import numpy as np
    pm = np.hypot(pmra, pmdec) #mas/yr
    offset = pm*baseline_year/1e3
    ```
    """
    pl.rcParams["font.size"] = fontsize
    pl.rcParams["xtick.labelsize"] = fontsize

    if (survey1 == "ps1") or (survey2 == "ps1"):
        try:
            import panstarrs3 as p3

            fov = np.hypot(width, height) * u.arcmin
            ps = p3.Panstarrs(
                ra=ra,
                dec=dec,
                fov=fov.to(u.arcsec),
                format="fits",
                color=False,
            )
            img, hdr = ps.get_fits(filter=filter, verbose=False)
        except Exception:
            raise ModuleNotFoundError(
                "pip install git+https://github.com/jpdeleon/panstarrs3.git"
            )

    # poss1
    if fp1 is not None and fp2 is not None:
        hdu1 = fits.open(fp1)[0]
        hdu2 = fits.open(fp2)[0]
    else:
        if survey1 == "ps1":
            hdu1 = fits.open(ps.get_url()[0])[0]
            hdu1.header["DATE-OBS"] = Time(
                hdu1.header["MJD-OBS"], format="mjd"
            ).strftime("%Y-%m-%d")
            hdu1.header["FILTER"] = hdu1.header["FPA.FILTER"].split(".")[0]
            hdu1.header["SURVEY"] = "Panstarrs1"
        else:
            hdu1 = get_dss_data(
                ra, dec, height=height, width=width, survey=survey1
            )
        if survey2 == "ps1":
            hdu2 = fits.open(ps.get_url()[0])[0]
            hdu2.header["DATE-OBS"] = Time(
                hdu2.header["MJD-OBS"], format="mjd"
            ).strftime("%Y-%m-%d")
            hdu2.header["FILTER"] = hdu2.header["FPA.FILTER"].split(".")[0]
            hdu2.header["SURVEY"] = "Panstarrs1"
        else:
            hdu2 = get_dss_data(
                ra, dec, height=height, width=width, survey=survey2
            )
    try:
        from reproject import reproject_interp
    except Exception:
        cmd = "pip install reproject"
        raise ModuleNotFoundError(cmd)

    projected_img, footprint = reproject_interp(hdu2, hdu1.header)

    fig = pl.figure(figsize=(10, 5), constrained_layout=False)
    interval = ZScaleInterval(contrast=contrast)

    # data1 = hdu1.data
    header1 = hdu1.header
    ax1 = fig.add_subplot("121", projection=WCS(header1))
    _ = plot_dss_image(
        hdu1, cmap=cmap, contrast=contrast, coord_format=coord_format, ax=ax1
    )
    if reticle:
        c = Circle(
            (ra, dec),
            0.001,
            edgecolor=color,
            facecolor="none",
            lw=2,
            transform=ax1.get_transform("fk5"),
        )
        ax1.add_patch(c)
    filt1 = (
        hdu1.header["FILTER"]
        if hdu1.header["FILTER"] is not None
        else survey1.split("_")[1]
    )
    title = f"{header1['SURVEY']} ({filt1})\n"
    title += f"{header1['DATE-OBS'][:10]}"
    ax1.set_title(title)
    # set RA from hourangle to degree
    if hasattr(ax1, "coords"):
        ax1.coords[0].set_major_formatter(coord_format)
        ax1.coords[1].set_major_formatter(coord_format)

    # recent
    data2, header2 = hdu2.data, hdu2.header
    ax2 = fig.add_subplot("122", projection=WCS(header1))
    # _ = plot_dss_image(hdu2, ax=ax2)
    zmin, zmax = interval.get_limits(data2)
    ax2.imshow(projected_img, origin="lower", vmin=zmin, vmax=zmax, cmap=cmap)
    if reticle:
        c = Circle(
            (ra, dec),
            0.001,
            edgecolor=color,
            facecolor="none",
            lw=2,
            transform=ax2.get_transform("fk5"),
        )
        ax2.add_patch(c)
        # ax2.scatter(ra, dec, 'r+')
    filt2 = (
        hdu2.header["FILTER"]
        if hdu2.header["FILTER"] is not None
        else survey2.split("_")[1]
    )
    ax2.coords["dec"].set_axislabel_position("r")
    ax2.coords["dec"].set_ticklabel_position("r")
    ax2.coords["dec"].set_axislabel("DEC")
    ax2.set_xlabel("RA")
    title = f"{header2['SURVEY']} ({filt2})\n"
    title += f"{header2['DATE-OBS'][:10]}"
    ax2.set_title(title)
    # set RA from hourangle to degree
    if hasattr(ax2, "coords"):
        ax2.coords[0].set_major_formatter(coord_format)
        ax2.coords[1].set_major_formatter(coord_format)

    if grid:
        [ax.grid(True) for ax in fig.axes]
    fig.tight_layout(rect=[0, 0.03, 0.5, 0.9])
    fig.suptitle(".", y=0.995)
    fig.tight_layout()
    if return_baseline:
        baseline = int(header2["DATE-OBS"][:4]) - int(header1["DATE-OBS"][:4])
        return fig, baseline
    else:
        return fig
