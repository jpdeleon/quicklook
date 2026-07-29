import json
import os
import sys
from pathlib import Path
from urllib.request import urlopen
import numpy as np
import pandas as pd
import astropy.units as u
from loguru import logger

from importlib.resources import files

DATA_PATH = files("quicklook").joinpath("data")
TESS_TIME_OFFSET = 2_457_000
TESS_pix_scale = 21 * u.arcsec  # / u.pixel
# K2_TIME_OFFSET = 2_454_833  # BKJD
# Kepler_pix_scale = 3.98 * u.arcsec  # /pix

__all__ = [
    "get_available_sectors",
    "get_exofop_json",
    "get_params_from_exofop",
    "parse_aperture_mask",
    "compute_secthresh",
    "mag_to_flux",
    "is_point_inside_mask",
    "PadWithZeros",
]


def _get_tglc_sectors_via_tesscut(target_name, tic_id=None):
    """Enumerate TESS sectors covering a target using the TESScut service.

    Avoids the MAST Observations portal (`/api/v0/invoke`) so it survives
    outages of that subsystem. TGLC runs locally on FFI cutouts, so sector
    availability is driven purely by sky coverage, not by any pipeline
    products existing on MAST.
    """
    from astropy.coordinates import SkyCoord
    from astroquery.mast import Tesscut

    if tic_id is not None:
        resolve_name = f"TIC {tic_id}"
    else:
        resolve_name = target_name.replace("_", " ").replace("-", " ")

    coord = SkyCoord.from_name(resolve_name)
    sector_table = Tesscut.get_sectors(coordinates=coord)
    if len(sector_table) == 0:
        return []
    # MAST sometimes returns synthetic/HLSP composite sectorNames whose
    # ``sector`` value is not a real observing sector (e.g. ``tess-s1751-1-3``
    # for TOI-5169). Filter those out so they don't pollute each-sector
    # queues with downstream-incompatible sector numbers.
    MAX_REAL_SECTOR = 300
    return sorted({int(s) for s in sector_table["sector"] if 1 <= int(s) <= MAX_REAL_SECTOR})


def get_available_sectors(target_name, pipeline="SPOC", tic_id=None):
    """Query available sectors for target+pipeline.

    Parameters
    ----------
    target_name : str
        Target name (e.g., "TOI-1234" or "TIC123456")
    pipeline : str
        Pipeline name (e.g., "SPOC", "QLP", "TGLC")
    tic_id : int, optional
        Pre-resolved TIC ID to skip the ExoFOP lookup.

    Returns
    -------
    list
        Sorted list of unique sector numbers available for the target and pipeline.
    """
    import lightkurve as lk

    if pipeline.upper() == "TGLC":
        return _get_tglc_sectors_via_tesscut(target_name, tic_id=tic_id)

    # Resolve to TIC ID via ExoFOP, matching TessQuickLook behavior
    if tic_id is not None:
        query_name = f"TIC{tic_id}"
    else:
        try:
            tic_id = get_tic_id(target_name)
            query_name = f"TIC{tic_id}"
        except (ValueError, KeyError, OSError):
            query_name = target_name

    search = lk.search_lightcurve(query_name)
    if len(search) == 0:
        return []

    df = search.table.to_pandas()

    if pipeline.upper() == "SPOC":
        sectors = []
        for mission in df["mission"].tolist():
            x = mission.split()
            if len(x) == 3:
                sectors.append(int(x[-1]))
        return sorted(set(sectors))

    sectors = []
    for mission, author in zip(df["mission"].tolist(), df["provenance_name"].tolist()):
        if author.lower() == pipeline.lower():
            x = mission.split()
            if len(x) == 3:
                sectors.append(int(x[-1]))

    return sorted(set(sectors))


def get_available_pipelines(target_name, tic_id=None):
    """Query available pipelines (HLSP / mission products) for a target.

    Parameters
    ----------
    target_name : str
        Target name (e.g. "TOI-1234" or "TIC123456").
    tic_id : int, optional
        Pre-resolved TIC ID to skip the ExoFOP lookup.

    Returns
    -------
    list[str]
        Sorted list of unique pipeline names (lowercase) available for
        the target. Always includes ``"tglc"`` because the local ePSF
        fallback in ``TessQuickLook`` can produce a TGLC light curve
        even when MAST lists no TGLC HLSP product.
    """
    import lightkurve as lk

    if tic_id is not None:
        query_name = f"TIC{tic_id}"
    else:
        try:
            tic_id = get_tic_id(target_name)
            query_name = f"TIC{tic_id}"
        except (ValueError, KeyError, OSError):
            query_name = target_name

    search = lk.search_lightcurve(query_name)
    if len(search) == 0:
        # Local TGLC ePSF fallback is always an option even when MAST
        # has no products for the target.
        return ["tglc"]

    df = search.table.to_pandas()
    pipelines = sorted({p.lower() for p in df["provenance_name"].tolist()})
    if "tglc" not in pipelines:
        pipelines.append("tglc")
        pipelines.sort()
    return pipelines


def get_tois(
    clobber=False,
    outdir=DATA_PATH,
    verbose=False,
    remove_FP=True,
    remove_known_planets=False,
    add_FPP=False,
):
    """Download TOI list from TESS Alert/TOI Release.

    Parameters
    ----------
    clobber : bool
        re-download table and save as csv file
    outdir : str
        download directory location
    verbose : bool
        print texts

    Returns
    -------
    d : pandas.DataFrame
        TOI table as dataframe
    """
    dl_link = "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv"
    fp = Path(outdir, "TOIs.csv")
    if not Path(outdir).exists():
        Path(outdir).mkdir()

    if not Path(fp).exists() or clobber:
        d = pd.read_csv(dl_link)  # , dtype={'RA': float, 'Dec': float})
        msg = f"Downloading {dl_link}\n"
        if add_FPP:
            fp2 = Path(outdir, "Giacalone2020/tab4.txt")
            classified = ascii.read(fp2).to_pandas()
            fp3 = Path(outdir, "Giacalone2020/tab5.txt")
            unclassified = ascii.read(fp3).to_pandas()
            fpp = pd.concat(
                [
                    classified[["TOI", "FPP-2m", "FPP-30m"]],
                    unclassified[["TOI", "FPP"]],
                ],
                sort=True,
            )
            d = pd.merge(d, fpp, how="outer").drop_duplicates()
        d.to_csv(fp, index=False)
    else:
        d = pd.read_csv(fp).drop_duplicates()
        msg = f"Loaded: {fp}\n"
    assert len(d) > 1000, f"{fp} likely has been overwritten!"

    # remove False Positives
    if remove_FP:
        d = d[d["TFOPWG Disposition"] != "FP"]
        msg += "TOIs with TFPWG disposition==FP are removed.\n"
    if remove_known_planets:
        planet_keys = [
            "HD",
            "GJ",
            "LHS",
            "XO",
            "WASP",
            "SWASP",
            "HAT",
            "HATS",
            "KELT",
            "TrES",
            "QATAR",
            "CoRoT",
            "K2",  # , "EPIC"
            "Kepler",  # "KOI"
        ]
        keys = []
        for key in planet_keys:
            idx = ~np.array(d["Comments"].str.contains(key).tolist(), dtype=bool)
            d = d[idx]
            if idx.sum() > 0:
                keys.append(key)
        msg += f"{keys} planets are removed.\n"
    msg += f"Saved: {fp}\n"
    if verbose:
        logger.info(msg)
    return d.sort_values("TOI")


def get_exofop_json(target_name: str) -> dict:
    base_url = "https://exofop.ipac.caltech.edu/tess"
    if target_name.lower()[:4] == "gaia":
        query_name = target_name.replace(" ", "_")
    else:
        query_name = target_name.replace(" ", "")
    url = f"{base_url}/target.php?id={query_name}&json"
    response = urlopen(url, timeout=30)
    if response.code != 200:
        raise ConnectionError(f"ExoFOP-TESS returned HTTP {response.code} for {target_name}")
    try:
        data_json = json.loads(response.read())
        return data_json
    except (json.JSONDecodeError, UnicodeDecodeError) as e:
        raise ValueError(f"No TIC data found for {target_name}") from e


def get_params_from_exofop(exofop_info, name="planet_parameters", idx=None):
    params_dict = exofop_info.get(name)
    if not params_dict:
        return None
    try:
        if idx is None:
            key = "pdate" if name == "planet_parameters" else "sdate"
            # get the latest parameter based on upload date
            dates = [d.get(key) for d in params_dict]
            df = pd.DataFrame({"date": dates})
            df["date"] = pd.to_datetime(df["date"], format="mixed", errors="coerce")
            if df["date"].isna().all():
                idx = 0  # no parseable dates; fall back to first entry
            else:
                idx = df["date"].idxmax()
        return params_dict[idx]
    except (KeyError, TypeError, IndexError, ValueError) as e:
        logger.warning(f"Could not parse ExoFOP params ({name}): {e}")
        return None


def get_tic_id(target_name: str) -> int:
    return int(get_exofop_json(target_name)["basic_info"]["tic_id"])


def get_toi_ephem(target_name: str, idx=1, params=["epoch", "per", "dur"]) -> list:
    logger.info(f"Querying ephemeris for {target_name}:")
    r = get_exofop_json(target_name)
    planet_params = r["planet_parameters"][idx]
    vals = []
    for p in params:
        val = planet_params.get(p)
        val = float(val) if val else 0.1
        err = planet_params.get(p + "_e")
        err = float(err) if err else 0.1
        logger.info(f"     {p}: {val}, {err}")
        vals.append((val, err))
    return vals


def mag_to_flux(mag: np.ndarray, mag_err: np.ndarray = None):
    """
    Convert magnitudes and magnitude uncertainties to relative flux and flux uncertainties.

    Parameters
    ----------
    mag : np.ndarray
        Magnitude values.
    mag_err : np.ndarray
        Magnitude uncertainties.

    Returns
    -------
    flux : np.ndarray
        Relative flux values (normalized to 1 at mag=0 reference).
    flux_err : np.ndarray
        Corresponding flux uncertainties.
    """
    # Convert magnitude to relative flux (zero-point cancels in relative space)
    flux = 10 ** (-0.4 * mag)

    if mag_err is not None:
        # Propagate uncertainty analytically
        flux_err = 0.921034 * flux * mag_err
        return flux, flux_err
    else:
        return flux


def parse_aperture_mask(
    tpf,
    sap_mask="pipeline",
    aper_radius=None,
    percentile=None,
    verbose=False,
    threshold_sigma=None,
):
    """Parse and make aperture mask"""
    if verbose:
        if sap_mask == "round":
            logger.info(f"aperture photometry mask: {sap_mask} (r={aper_radius} pix)")
        elif sap_mask == "square":
            logger.info(f"aperture photometry mask: {sap_mask} ({aper_radius}x{aper_radius} pix)")
        elif sap_mask == "percentile":
            logger.info(f"aperture photometry mask: {sap_mask} ({percentile}%)")
        else:
            logger.info(f"aperture photometry mask: {sap_mask}")

    median_img = np.nanmedian(tpf.flux, axis=0).value
    if (sap_mask == "pipeline") or (sap_mask is None):
        errmsg = "tpf does not have pipeline mask"
        assert tpf.pipeline_mask is not None, errmsg
        mask = tpf.pipeline_mask  # default
    elif sap_mask == "all":
        mask = np.ones((tpf.shape[1], tpf.shape[2]), dtype=bool)
    elif sap_mask == "round":
        assert aper_radius is not None, "supply aper_radius"
        mask = make_round_mask(median_img, radius=aper_radius)
    elif sap_mask == "square":
        assert aper_radius is not None, "supply aper_radius/size"
        mask = make_square_mask(median_img, size=aper_radius)
    elif sap_mask == "threshold":
        assert threshold_sigma is not None, "supply threshold_sigma"
        # FIXME: make sure aperture is contiguous
        mask = tpf.create_threshold_mask(threshold_sigma)
    elif sap_mask == "percentile":
        assert percentile is not None, "supply percentile"
        mask = median_img > np.nanpercentile(median_img, percentile)
    else:
        raise ValueError("Unknown aperture mask")
    return mask


def make_round_mask(img, radius, xy_center=None):
    """Make round mask in units of pixels

    Parameters
    ----------
    img : numpy ndarray
        image
    radius : int
        aperture mask radius or size
    xy_center : tuple
        aperture mask center position

    Returns
    -------
    mask : np.ma.masked_array
        aperture mask
    """
    offset = 2  # from center
    xcen, ycen = img.shape[0] // 2, img.shape[1] // 2
    if xy_center is None:  # use the middle of the image
        y, x = np.unravel_index(np.argmax(img), img.shape)
        xy_center = [x, y]
        # check if near edge
        if np.any([abs(x - xcen) > offset, abs(y - ycen) > offset]):
            logger.info("Brightest star is detected far from the center.")
            logger.info("Aperture mask is placed at the center instead.")
            xy_center = [xcen, ycen]

    Y, X = np.ogrid[: img.shape[0], : img.shape[1]]
    dist_from_center = np.sqrt((X - xy_center[0]) ** 2 + (Y - xy_center[1]) ** 2)

    mask = dist_from_center <= radius
    return np.ma.masked_array(img, mask=mask).mask


def make_square_mask(img, size, xy_center=None):
    """Make rectangular mask with optional rotation

    Parameters
    ----------
    img : numpy ndarray
        image
    size : int
        aperture mask size
    xy_center : tuple
        aperture mask center position
    angle : int
        rotation

    Returns
    -------
    mask : np.ma.masked_array
        aperture mask
    """
    offset = 2  # from center
    xcen, ycen = img.shape[0] // 2, img.shape[1] // 2
    if xy_center is None:  # use the middle of the image
        y, x = np.unravel_index(np.argmax(img), img.shape)
        xy_center = [x, y]
        # check if near edge
        if np.any([abs(x - xcen) > offset, abs(y - ycen) > offset]):
            logger.info("Brightest star detected is far from the center.")
            logger.info("Aperture mask is placed at the center instead.")
            xy_center = [xcen, ycen]
    mask = np.zeros_like(img, dtype=bool)
    mask[ycen - size : ycen + size + 1, xcen - size : xcen + size + 1] = True  # noqa
    # if angle:
    #    #rotate mask
    #    mask = rotate(mask, angle, axes=(1, 0),
    #                  reshape=True, output=bool, order=0)
    return mask


def compute_secthresh(fold_lc, t14):
    """
    Similar to Mayo+2018, compute `secthresh` by binning the phase-folded
    lightcurves by measuring the transit duration and taking thrice the value
    of the standard deviation of the mean in each bin.
    """
    means = []
    start, end = -0.5, 0.5
    chunks = np.arange(start, end, t14)
    for n, x in enumerate(chunks):
        if n == 0:
            x1 = start
            x2 = x
        elif n == len(chunks):
            x1 = x
            x2 = end
        else:
            x1 = chunks[n - 1]
            x2 = x
        idx = (fold_lc.phase.value > x1) & (fold_lc.phase.value < x2)
        if sum(idx) > 3:
            mean = np.nanmean(fold_lc.flux[idx].value)
            # print(mean)
            means.append(mean)
    return 3 * np.nanstd(means)


def get_cartersian_distance(x1, y1, x2, y2):
    return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)


def is_point_inside_mask(border, target):
    """determine if target coordinate is within polygon border"""
    degree = 0
    for i in range(len(border) - 1):
        a = border[i]
        b = border[i + 1]

        # calculate distance of vector
        A = get_cartersian_distance(a[0], a[1], b[0], b[1])
        B = get_cartersian_distance(target[0], target[1], a[0], a[1])
        C = get_cartersian_distance(target[0], target[1], b[0], b[1])

        # calculate direction of vector
        ta_x = a[0] - target[0]
        ta_y = a[1] - target[1]
        tb_x = b[0] - target[0]
        tb_y = b[1] - target[1]

        cross = tb_y * ta_x - tb_x * ta_y
        clockwise = cross < 0

        # calculate sum of angles
        if clockwise:
            degree = degree + np.rad2deg(np.arccos((B * B + C * C - A * A) / (2.0 * B * C)))
        else:
            degree = degree - np.rad2deg(np.arccos((B * B + C * C - A * A) / (2.0 * B * C)))

    if abs(round(degree) - 360) <= 3:
        return True
    return False


def PadWithZeros(vector, pad_width, iaxis, kwargs):
    vector[: pad_width[0]] = 0  # noqa
    vector[-pad_width[1] :] = 0  # noqa
    return vector
