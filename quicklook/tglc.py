"""
This code was copied verbatim from:
https://github.com/TeHanHunter/TESS_Gaia_Light_Curve

TGLC workflow


"""

import json
import sys
import numpy as np
import warnings
import pickle
import os
import requests
from os.path import exists
from urllib.parse import quote as urlencode
import time

from astropy.wcs import WCS
import astropy.units as u
from astropy.table import Table, hstack, vstack, unique, Column
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astroquery.mast import Catalogs, Tesscut
from wotan import flatten
from tqdm import trange

Gaia.ROW_LIMIT = -1
Gaia.MAIN_GAIA_TABLE = (
    "gaiadr3.gaia_source"  # TODO: dr3 MJD = 2457388.5, TBJD = 388.5
)


def convert_gaia_id(catalogdata_tic):
    query = """
            SELECT dr2_source_id, dr3_source_id
            FROM gaiadr3.dr2_neighbourhood
            WHERE dr2_source_id IN {gaia_ids}
            """
    gaia_array = np.array(
        [str(item) for item in catalogdata_tic["GAIA"]], dtype=object
    )
    gaia_array = gaia_array[gaia_array != "None"]
    # np.save('gaia_array.npy', gaia_array)
    segment = (len(gaia_array) - 1) // 10000
    gaia_tuple = tuple(gaia_array[:10000])
    results = Gaia.launch_job_async(
        query.format(gaia_ids=gaia_tuple)
    ).get_results()
    # np.save('result.npy', np.array(results))
    for i in range(segment):
        gaia_array_cut = gaia_array[((i + 1) * 10000) : ((i + 2) * 10000)]
        gaia_tuple_cut = tuple(gaia_array_cut)
        results = vstack(
            [
                results,
                Gaia.launch_job_async(
                    query.format(gaia_ids=gaia_tuple_cut)
                ).get_results(),
            ]
        )
    tic_ids = []
    for j in range(len(results)):
        tic_ids.append(
            int(
                catalogdata_tic["ID"][
                    np.where(
                        catalogdata_tic["GAIA"]
                        == str(results["dr2_source_id"][j])
                    )
                ][0]
            )
        )
    tic_ids = Column(np.array(tic_ids), name="TIC")
    results.add_column(tic_ids)
    return results


def mast_query(request):
    """Perform a MAST query.

    Parameters
    ----------
    request (dictionary): The MAST request json object

    Returns head,content where head is the response HTTP headers, and content is the returned data
    """

    # Base API url
    request_url = "https://mast.stsci.edu/api/v0/invoke"
    # Grab Python Version
    version = ".".join(map(str, sys.version_info[:3]))
    # Create Http Header Variables
    headers = {
        "Content-type": "application/x-www-form-urlencoded",
        "Accept": "text/plain",
        "User-agent": "python-requests/" + version,
    }
    # Encoding the request as a json string
    req_string = json.dumps(request)
    req_string = urlencode(req_string)
    # Perform the HTTP request
    resp = requests.post(
        request_url, data="request=" + req_string, headers=headers
    )
    # Pull out the headers and response content
    head = resp.headers
    content = resp.content.decode("utf-8")
    return head, content


def mast_json2table(json_obj):
    data_table = Table()
    for col, atype in [(x["name"], x["type"]) for x in json_obj["fields"]]:
        if atype == "string":
            atype = "str"
        if atype == "boolean":
            atype = "bool"
        data_table[col] = np.array(
            [x.get(col, None) for x in json_obj["data"]], dtype=atype
        )
    return data_table


def tic_advanced_search_position_rows(
    ra=1.0, dec=1.0, radius=0.5, limit_mag=16
):
    request = {
        "service": "Mast.Catalogs.Filtered.Tic.Position.Rows",
        "format": "json",
        "params": {
            "columns": "ID, GAIA",
            "filters": [
                {
                    "paramName": "Tmag",
                    "values": [{"min": -10.0, "max": (limit_mag + 0.5)}],
                }
            ],
            "ra": ra,
            "dec": dec,
            "radius": radius,
        },
    }

    headers, out_string = mast_query(request)
    out_data = json.loads(out_string)
    return mast_json2table(out_data)


class Source(object):
    def __init__(
        self,
        x=0,
        y=0,
        flux=None,
        time=None,
        wcs=None,
        quality=None,
        mask=None,
        exposure=1800,
        sector=0,
        size=150,
        camera=1,
        ccd=1,
        cadence=None,
    ):
        """
        Source object that includes all data from TESS and Gaia DR2
        :param x: int, required
        starting horizontal pixel coordinate
        :param y: int, required
        starting vertical pixel coordinate
        :param flux: np.ndarray, required
        3d data cube, the time series of a all FFI from a CCD
        :param time: np.ndarray, required
        1d array of time
        :param wcs: astropy.wcs.wcs.WCS, required
        WCS Keywords of the TESS FFI
        :param sector: int, required
        TESS sector number
        :param size: int, optional
        the side length in pixel  of TESScut image
        :param camera: int, optional
        camera number
        :param ccd: int, optional
        CCD number
        :param cadence: list, required
        list of cadences of TESS FFI
        """
        super(Source, self).__init__()
        if cadence is None:
            cadence = []
        if quality is None:
            quality = []
        if wcs is None:
            wcs = []
        if time is None:
            time = []
        if flux is None:
            flux = []

        self.size = size
        self.sector = sector
        self.camera = camera
        self.ccd = ccd
        self.cadence = cadence
        self.quality = quality
        self.exposure = exposure
        self.wcs = wcs
        co1 = 38.5
        co2 = 116.5
        catalog_1 = self.search_gaia(x, y, co1, co1)
        catalog_2 = self.search_gaia(x, y, co1, co2)
        catalog_3 = self.search_gaia(x, y, co2, co1)
        catalog_4 = self.search_gaia(x, y, co2, co2)
        catalogdata = vstack(
            [catalog_1, catalog_2, catalog_3, catalog_4], join_type="exact"
        )
        catalogdata = unique(catalogdata, keys="DESIGNATION")
        coord = wcs.pixel_to_world(
            [x + (size - 1) / 2 + 44], [y + (size - 1) / 2]
        )[0].to_string()
        ra = float(coord.split()[0])
        dec = float(coord.split()[1])
        catalogdata_tic = tic_advanced_search_position_rows(
            ra=ra, dec=dec, radius=(self.size + 2) * 21 * 0.707 / 3600
        )
        # print(f'no_of_stars={len(catalogdata_tic)}, camera={camera}, ccd={ccd}: ra={ra}, dec={dec}, radius={(self.size + 2) * 21 * 0.707 / 3600}')
        self.tic = convert_gaia_id(catalogdata_tic)
        self.flux = flux[:, y : y + size, x : x + size]
        self.mask = mask[y : y + size, x : x + size]
        self.time = np.array(time)
        median_time = np.median(self.time)
        interval = (median_time - 388.5) / 365.25
        # Julian Day Number:	2457000.0 (TBJD=0)
        # Calendar Date/Time:	2014-12-08 12:00:00 388.5 days before J2016

        num_gaia = len(catalogdata)
        tic_id = np.zeros(num_gaia)
        x_gaia = np.zeros(num_gaia)
        y_gaia = np.zeros(num_gaia)
        tess_mag = np.zeros(num_gaia)
        in_frame = [True] * num_gaia
        for i, designation in enumerate(catalogdata["DESIGNATION"]):
            ra = catalogdata["ra"][i]
            dec = catalogdata["dec"][i]
            if not np.isnan(catalogdata["pmra"].mask[i]):  # masked?
                ra += (
                    catalogdata["pmra"][i]
                    * np.cos(np.deg2rad(dec))
                    * interval
                    / 1000
                    / 3600
                )
            if not np.isnan(catalogdata["pmdec"].mask[i]):
                dec += catalogdata["pmdec"][i] * interval / 1000 / 3600
            pixel = self.wcs.all_world2pix(
                np.array(
                    [catalogdata["ra"][i], catalogdata["dec"][i]]
                ).reshape((1, 2)),
                0,
                quiet=True,
            )
            x_gaia[i] = pixel[0][0] - x - 44
            y_gaia[i] = pixel[0][1] - y
            try:
                tic_id[i] = catalogdata_tic["ID"][
                    np.where(
                        catalogdata_tic["GAIA"] == designation.split()[2]
                    )[0][0]
                ]
            except Exception as e:
                print(e)
                tic_id[i] = np.nan
            if np.isnan(catalogdata["phot_g_mean_mag"][i]):
                in_frame[i] = False
            elif catalogdata["phot_g_mean_mag"][i] >= 25:
                in_frame[i] = False
            elif (
                -4 < x_gaia[i] < self.size + 3
                and -4 < y_gaia[i] < self.size + 3
            ):
                dif = (
                    catalogdata["phot_bp_mean_mag"][i]
                    - catalogdata["phot_rp_mean_mag"][i]
                )
                tess_mag[i] = (
                    catalogdata["phot_g_mean_mag"][i]
                    - 0.00522555 * dif**3
                    + 0.0891337 * dif**2
                    - 0.633923 * dif
                    + 0.0324473
                )
                if np.isnan(tess_mag[i]):
                    tess_mag[i] = catalogdata["phot_g_mean_mag"][i] - 0.430
                if np.isnan(tess_mag[i]):
                    in_frame[i] = False
            else:
                in_frame[i] = False

        tess_flux = 10 ** (-tess_mag / 2.5)
        t = Table()
        t["tess_mag"] = tess_mag[in_frame]
        t["tess_flux"] = tess_flux[in_frame]
        t["tess_flux_ratio"] = tess_flux[in_frame] / np.nanmax(
            tess_flux[in_frame]
        )
        t["sector_{self.sector}_x"] = x_gaia[in_frame]
        t["sector_{self.sector}_y"] = y_gaia[in_frame]
        catalogdata = hstack([catalogdata[in_frame], t])
        catalogdata.sort("tess_mag")
        self.gaia = catalogdata

    def search_gaia(self, x, y, co1, co2):
        coord = self.wcs.pixel_to_world([x + co1 + 44], [y + co2])[
            0
        ].to_string()
        radius = u.Quantity((self.size / 2 + 4) * 21 * 0.707 / 3600, u.deg)
        attempt = 0
        while attempt < 5:
            try:
                catalogdata = Gaia.cone_search_async(
                    coord,
                    radius=radius,
                    columns=[
                        "DESIGNATION",
                        "phot_g_mean_mag",
                        "phot_bp_mean_mag",
                        "phot_rp_mean_mag",
                        "ra",
                        "dec",
                        "pmra",
                        "pmdec",
                    ],
                ).get_results()
                return catalogdata
            except Exception as e:
                print(e)
                attempt += 1
                time.sleep(10)
                print(
                    f"Trying Gaia search again. Coord = {coord}, radius = {radius}"
                )


class Source_cut(object):
    def __init__(
        self,
        name,
        size=50,
        sector=None,
        cadence=None,
        limit_mag=None,
        transient=None,
    ):
        """
        Source_cut object that includes all data from TESS and Gaia DR3
        :param name: str, required
        Target identifier (e.g. "NGC 7654" or "M31"),
        or coordinate in the format of ra dec (e.g. '351.40691 61.646657')
        :param size: int, optional
        The side length in pixel  of TESScut image
        :param cadence: list, required
        list of cadences of TESS FFI
        """
        super(Source_cut, self).__init__()
        if cadence is None:
            cadence = []
        if size < 25:
            warnings.warn("FFI cut size too small, try at least 25*25")
        self.name = name
        self.size = size
        self.sector = 0
        self.wcs = []
        self.time = []
        self.flux = []
        self.flux_err = []
        self.gaia = []
        self.cadence = cadence
        self.quality = []
        self.mask = []
        self.transient = transient

        target = Catalogs.query_object(
            self.name, radius=21 * 0.707 / 3600, catalog="Gaia", version=2
        )
        if len(target) == 0:
            target = Catalogs.query_object(
                self.name,
                radius=5 * 21 * 0.707 / 3600,
                catalog="Gaia",
                version=2,
            )
        ra = target[0]["ra"]
        dec = target[0]["dec"]
        coord = SkyCoord(
            ra=ra, dec=dec, unit=(u.degree, u.degree), frame="icrs"
        )
        radius = u.Quantity((self.size + 6) * 21 * 0.707 / 3600, u.deg)
        print(f"Target Gaia: {target[0]['designation']}")
        catalogdata = Gaia.cone_search_async(
            coord,
            radius=radius,
            columns=[
                "DESIGNATION",
                "phot_g_mean_mag",
                "phot_bp_mean_mag",
                "phot_rp_mean_mag",
                "ra",
                "dec",
                "pmra",
                "pmdec",
            ],
        ).get_results()
        print(f"Found {len(catalogdata)} Gaia DR3 objects.")
        catalogdata_tic = tic_advanced_search_position_rows(
            ra=ra,
            dec=dec,
            radius=(self.size + 2) * 21 * 0.707 / 3600,
            limit_mag=limit_mag,
        )
        print(f"Found {len(catalogdata_tic)} TIC objects.")
        self.tic = convert_gaia_id(catalogdata_tic)
        sector_table = Tesscut.get_sectors(coordinates=coord)
        if len(sector_table) == 0:
            warnings.warn("TESS has not observed this position yet :(")
        if sector is None:
            hdulist = Tesscut.get_cutouts(coordinates=coord, size=self.size)
        elif sector == "first":
            hdulist = Tesscut.get_cutouts(
                coordinates=coord,
                size=self.size,
                sector=sector_table["sector"][0],
            )
            sector = sector_table["sector"][0]
        elif sector == "last":
            hdulist = Tesscut.get_cutouts(
                coordinates=coord,
                size=self.size,
                sector=sector_table["sector"][-1],
            )
            sector = sector_table["sector"][-1]
        else:
            hdulist = Tesscut.get_cutouts(
                coordinates=coord, size=self.size, sector=sector
            )
        self.catalogdata = catalogdata
        self.sector_table = sector_table
        self.camera = int(sector_table[0]["camera"])
        self.ccd = int(sector_table[0]["ccd"])
        self.hdulist = hdulist
        sector_list = []
        for i in range(len(hdulist)):
            sector_list.append(hdulist[i][0].header["SECTOR"])
        self.sector_list = sector_list
        if sector is None:
            self.select_sector(sector=sector_table["sector"][0])
        else:
            self.select_sector(sector=sector)

    def select_sector(self, sector=1):
        """
        select sector to use if target is in multi-sectors
        :param sector: int, required
        TESS sector number
        """
        if self.sector == sector:
            print(f"Already in sector {sector}.")
            return
        elif sector not in self.sector_table["sector"]:
            print(
                f"Sector {sector} does not cover this region. Please refer to sector table."
            )
            return

        index = self.sector_list.index(sector)
        self.sector = sector
        hdu = self.hdulist[index]
        self.camera = int(hdu[0].header["CAMERA"])
        self.ccd = int(hdu[0].header["CCD"])
        wcs = WCS(hdu[2].header)
        data_time = hdu[1].data["TIME"]
        data_flux = hdu[1].data["FLUX"]
        data_flux_err = hdu[1].data["FLUX_ERR"]
        data_quality = hdu[1].data["QUALITY"]
        # data_time = data_time[np.where(data_quality == 0)]
        # data_flux = data_flux[np.where(data_quality == 0), :, :][0]
        # data_flux_err = data_flux_err[np.where(data_quality == 0), :, :][0]
        self.wcs = wcs
        self.time = data_time
        self.flux = data_flux
        self.flux_err = data_flux_err
        self.quality = data_quality
        median_time = np.median(data_time)
        interval = (median_time - 388.5) / 365.25

        mask = np.ones(np.shape(data_flux[0]))
        bad_pixels = np.zeros(np.shape(data_flux[0]))
        med_flux = np.median(data_flux, axis=0)
        bad_pixels[med_flux > 0.8 * np.nanmax(med_flux)] = 1
        bad_pixels[med_flux < 0.2 * np.nanmedian(med_flux)] = 1
        bad_pixels[np.isnan(med_flux)] = 1
        mask = np.ma.masked_array(mask, mask=bad_pixels)
        self.mask = mask

        gaia_targets = self.catalogdata[
            "DESIGNATION",
            "phot_g_mean_mag",
            "phot_bp_mean_mag",
            "phot_rp_mean_mag",
            "ra",
            "dec",
            "pmra",
            "pmdec",
        ]

        # inject transients
        if self.transient is not None:
            gaia_targets.add_row(
                [
                    self.transient[0],
                    20,
                    20,
                    20,
                    self.transient[1],
                    self.transient[2],
                    0,
                    0,
                ]
            )

        gaia_targets["phot_bp_mean_mag"].fill_value = np.nan
        gaia_targets["phot_rp_mean_mag"].fill_value = np.nan
        gaia_targets["pmra"].fill_value = np.nan
        gaia_targets["pmdec"].fill_value = np.nan
        gaia_targets = gaia_targets.filled()
        num_gaia = len(gaia_targets)
        # tic_id = np.zeros(num_gaia)
        x_gaia = np.zeros(num_gaia)
        y_gaia = np.zeros(num_gaia)
        tess_mag = np.zeros(num_gaia)
        in_frame = [True] * num_gaia
        for i, designation in enumerate(gaia_targets["DESIGNATION"]):
            ra = gaia_targets["ra"][i]
            dec = gaia_targets["dec"][i]
            if not np.isnan(gaia_targets["pmra"][i]):
                ra += (
                    gaia_targets["pmra"][i]
                    * np.cos(np.deg2rad(dec))
                    * interval
                    / 1000
                    / 3600
                )
            if not np.isnan(gaia_targets["pmdec"][i]):
                dec += gaia_targets["pmdec"][i] * interval / 1000 / 3600
            pixel = self.wcs.all_world2pix(
                np.array([ra, dec]).reshape((1, 2)), 0
            )
            x_gaia[i] = pixel[0][0]
            y_gaia[i] = pixel[0][1]
            if np.isnan(gaia_targets["phot_g_mean_mag"][i]):
                in_frame[i] = False
            elif gaia_targets["phot_g_mean_mag"][i] >= 25:
                in_frame[i] = False
            elif (
                -4 < x_gaia[i] < self.size + 3
                and -4 < y_gaia[i] < self.size + 3
            ):
                dif = (
                    gaia_targets["phot_bp_mean_mag"][i]
                    - gaia_targets["phot_rp_mean_mag"][i]
                )
                tess_mag[i] = (
                    gaia_targets["phot_g_mean_mag"][i]
                    - 0.00522555 * dif**3
                    + 0.0891337 * dif**2
                    - 0.633923 * dif
                    + 0.0324473
                )
                if np.isnan(tess_mag[i]):
                    tess_mag[i] = gaia_targets["phot_g_mean_mag"][i] - 0.430
            else:
                in_frame[i] = False
        tess_flux = 10 ** (-tess_mag / 2.5)
        t = Table()
        t["tess_mag"] = tess_mag[in_frame]
        t["tess_flux"] = tess_flux[in_frame]
        t["tess_flux_ratio"] = tess_flux[in_frame] / np.max(
            tess_flux[in_frame]
        )
        t["sector_{self.sector}_x"] = x_gaia[in_frame]
        t["sector_{self.sector}_y"] = y_gaia[in_frame]
        gaia_targets = hstack([gaia_targets[in_frame], t])
        if self.transient is not None:
            gaia_targets["tess_flux"][
                np.where(gaia_targets["DESIGNATION"] == self.transient[0])[0][
                    0
                ]
            ] = 0
        gaia_targets.sort("tess_mag")
        self.gaia = gaia_targets


def ffi_cut(
    target="",
    local_directory="",
    size=90,
    sector=None,
    limit_mag=None,
    transient=None,
):
    """
    Function to generate Source_cut objects
    :param target: string, required
    target name
    :param local_directory: string, required
    output directory
    :param size: int, required
    FFI cut side length
    :param sector: int, required
    TESS sector number
    :return: tglc.ffi_cut.Source_cut
    """
    if sector is None:
        source_name = f"source_{target}"
    elif sector == "first":
        source_name = f"source_{target}_earliest_sector"
    elif sector == "last":
        source_name = f"source_{target}_last_sector"
    else:
        source_name = f"source_{target}_sector_{sector}"
    source_exists = exists(f"{local_directory}source/{source_name}.pkl")
    if (
        source_exists
        and os.path.getsize(f"{local_directory}source/{source_name}.pkl") > 0
    ):
        with open(
            f"{local_directory}source/{source_name}.pkl", "rb"
        ) as input_:
            source = pickle.load(input_)
        print(source.sector_table)
        print("Loaded ffi_cut from directory. ")
    else:
        with open(
            f"{local_directory}source/{source_name}.pkl", "wb"
        ) as output:
            source = Source_cut(
                target,
                size=size,
                sector=sector,
                limit_mag=limit_mag,
                transient=transient,
            )
            pickle.dump(source, output, pickle.HIGHEST_PROTOCOL)
    return source


def lc_output(
    source,
    local_directory="",
    index=0,
    time=None,
    psf_lc=None,
    cal_psf_lc=None,
    aper_lc=None,
    cal_aper_lc=None,
    bg=None,
    tess_flag=None,
    tglc_flag=None,
    cadence=None,
    aperture=None,
    cut_x=None,
    cut_y=None,
    star_x=2,
    star_y=2,
    x_aperture=None,
    y_aperture=None,
    near_edge=False,
    local_bg=None,
    save_aper=False,
    portion=1,
    prior=None,
    transient=None,
    target_5x5=None,
    field_stars_5x5=None,
):
    """
    lc output to .FITS file in MAST HLSP standards
    :param tglc_flag: np.array(), required
    TGLC quality flags
    :param source: tglc.ffi_cut.Source or tglc.ffi_cut.Source_cut, required
    Source or Source_cut object
    :param local_directory: string, required
    output directory
    :param index: int, required
    star index
    :param time: list, required
    epochs of FFI
    :param lc: list, required
    ePSF light curve fluxes
    :param cal_lc: list, required
    ePSF light curve fluxes, detrended
    :param cadence: list, required
    list of cadences of TESS FFI
    :return:
    """
    if transient is None:
        objid = [
            int(s)
            for s in (source.gaia[index]["DESIGNATION"]).split()
            if s.isdigit()
        ][0]
    else:
        objid = transient[0]
    # source_path = f"{local_directory}hlsp_tglc_tess_ffi_gaiaid-{objid}-s{source.sector:04d}-cam{source.camera}-ccd{source.ccd}_tess_v1_llc.fits"
    # source_exists = exists(source_path)
    # if source_exists and (os.path.getsize(source_path) > 0):
    #     print('LC exists, please (re)move the file if you wish to overwrite.')
    #     return
    if np.isnan(source.gaia[index]["phot_bp_mean_mag"]) or np.ma.is_masked(
        source.gaia[index]["phot_bp_mean_mag"]
    ):
        gaia_bp = "NaN"
    else:
        gaia_bp = source.gaia[index]["phot_bp_mean_mag"]
    if np.isnan(source.gaia[index]["phot_rp_mean_mag"]) or np.ma.is_masked(
        source.gaia[index]["phot_rp_mean_mag"]
    ):
        gaia_rp = "NaN"
    else:
        gaia_rp = source.gaia[index]["phot_rp_mean_mag"]
    psf_err = 1.4826 * np.nanmedian(np.abs(psf_lc - np.nanmedian(psf_lc)))
    if np.isnan(psf_err):
        psf_err = "NaN"
    aper_err = 1.4826 * np.nanmedian(np.abs(aper_lc - np.nanmedian(aper_lc)))
    if np.isnan(aper_err):
        aper_err = "NaN"
    cal_psf_err = 1.4826 * np.nanmedian(
        np.abs(cal_psf_lc - np.nanmedian(cal_psf_lc))
    )
    if np.isnan(cal_psf_err):
        cal_psf_err = "NaN"
    cal_aper_err = 1.4826 * np.nanmedian(
        np.abs(cal_aper_lc - np.nanmedian(cal_aper_lc))
    )
    if np.isnan(cal_aper_err):
        cal_aper_err = "NaN"
    try:
        ticid = str(
            source.tic["TIC"][np.where(source.tic["dr3_source_id"] == objid)][
                0
            ]
        )
    except Exception as e:
        print(e)
        ticid = ""
    try:
        raw_flux = np.nanmedian(source.flux[:, star_y, star_x])
    except Exception as e:
        print(e)
        raw_flux = None
    if save_aper:
        primary_hdu = fits.PrimaryHDU(aperture)
    else:
        primary_hdu = fits.PrimaryHDU()
    # Simulated star images based on ePSF, used to estimate contamination ratio and others
    image_data = np.zeros((3, 5, 5))
    image_data[0] = target_5x5
    image_data[1] = field_stars_5x5
    # This is the pixel-wise contamination ratio
    image_data[2] = field_stars_5x5 / target_5x5
    image_hdu = fits.ImageHDU(data=image_data)

    primary_hdu.header = fits.Header(
        cards=[
            fits.Card("SIMPLE", True, "conforms to FITS standard"),
            fits.Card("EXTEND", True),
            fits.Card("NEXTEND", 1, "number of standard extensions"),
            fits.Card("EXTNAME", "PRIMARY", "name of extension"),
            fits.Card(
                "EXTDATA",
                "aperture",
                "decontaminated FFI cut for aperture photometry",
            ),
            fits.Card("EXTVER", 1, "extension version"),
            fits.Card("TIMESYS", "TDB", "TESS Barycentric Dynamical Time"),
            fits.Card("BUNIT", "e-/s", "flux unit"),
            fits.Card("STAR_X", x_aperture, "star x position in cut"),
            fits.Card("STAR_Y", y_aperture, "star y position in cut"),
            fits.Card("COMMENT", "hdul[0].data[:,star_y,star_x]=lc"),
            fits.Card(
                "ORIGIN",
                "UCSB/TGLC",
                "institution responsible for creating this file",
            ),
            fits.Card("TELESCOP", "TESS", "telescope"),
            fits.Card("INSTRUME", "TESS Photometer", "detector type"),
            fits.Card(
                "FILTER", "TESS", "the filter used for the observations"
            ),
            fits.Card(
                "OBJECT",
                source.gaia[index]["DESIGNATION"],
                "string version of Gaia DR3 ID",
            ),
            fits.Card("GAIADR3", objid, "integer version of Gaia DR3 ID"),
            fits.Card("TICID", ticid, "TESS Input Catalog ID"),
            fits.Card("SECTOR", source.sector, "observation sector"),
            fits.Card("CAMERA", source.camera, "camera No."),
            fits.Card("CCD", source.ccd, "CCD No."),
            fits.Card("CUT_x", cut_x, "FFI cut x index"),
            fits.Card("CUT_y", cut_y, "FFI cut y index"),
            fits.Card("CUTSIZE", source.size, "FFI cut size"),
            fits.Card(
                "RADESYS", "ICRS", "reference frame of celestial coordinates"
            ),
            fits.Card(
                "RA_OBJ",
                source.gaia[index]["ra"],
                "[deg] right ascension, J2000",
            ),
            fits.Card(
                "DEC_OBJ",
                source.gaia[index]["dec"],
                "[deg] declination, J2000",
            ),
            fits.Card(
                "TESSMAG",
                source.gaia[index]["tess_mag"],
                "TESS magnitude, fitted by Gaia DR3 bands",
            ),
            fits.Card(
                "GAIA_G",
                source.gaia[index]["phot_g_mean_mag"],
                "Gaia DR3 g band magnitude",
            ),
            fits.Card("GAIA_bp", gaia_bp, "Gaia DR3 bp band magnitude"),
            fits.Card("GAIA_rp", gaia_rp, "Gaia DR3 rp band magnitude"),
            fits.Card("RAWFLUX", raw_flux, "median flux of raw FFI"),
            fits.Card(
                "CONTAMRT",
                round(
                    np.nansum(field_stars_5x5[1:4, 1:4])
                    / np.nansum(target_5x5[1:4, 1:4]),
                    9,
                ),
                "contamination ratio of default 3*3 aperture",
            ),
            fits.Card("CALIB", "TGLC", "pipeline used for image calibration"),
        ]
    )
    if save_aper:
        primary_hdu.header.comments["NAXIS1"] = "Time (hdul[1].data['time'])"
        primary_hdu.header.comments["NAXIS2"] = "x size of cut"
        primary_hdu.header.comments["NAXIS3"] = "y size of cut"

    t_start = source.time[0]
    t_stop = source.time[-1]
    if source.sector < 27:  # primary
        exposure_time = 1800
    elif source.sector < 56:  # first extended
        exposure_time = 600
    else:  # second extended
        exposure_time = 200
    c1 = fits.Column(name="time", array=np.array(time), format="D")
    c2 = fits.Column(
        name="psf_flux", array=np.array(psf_lc), format="E"
    )  # psf factor
    # c3 = fits.Column(name='psf_flux_err',
    #                  array=1.4826 * np.median(np.abs(psf_lc - np.median(psf_lc))) * np.ones(len(psf_lc)), format='E')
    c4 = fits.Column(name="aperture_flux", array=aper_lc / portion, format="E")
    # c5 = fits.Column(name='aperture_flux_err',
    #                  array=1.4826 * np.median(np.abs(aper_lc - np.median(aper_lc))) * np.ones(len(aper_lc)), format='E')
    c6 = fits.Column(
        name="cal_psf_flux", array=np.array(cal_psf_lc), format="E"
    )
    # c7 = fits.Column(name='cal_psf_flux_err',
    #                  array=1.4826 * np.median(np.abs(cal_psf_lc - np.median(cal_psf_lc))) * np.ones(len(cal_psf_lc)),
    #                  format='E')
    c8 = fits.Column(
        name="cal_aper_flux", array=np.array(cal_aper_lc), format="E"
    )
    # c9 = fits.Column(name='cal_aper_flux_err',
    #                  array=1.4826 * np.median(np.abs(cal_aper_lc - np.median(cal_aper_lc))) * np.ones(len(cal_aper_lc)),
    #                  format='E')
    c10 = fits.Column(name="background", array=bg, format="E")  # add tilt
    c11 = fits.Column(
        name="cadence_num", array=np.array(cadence), format="J"
    )  # 32 bit int
    c12 = fits.Column(
        name="TESS_flags", array=np.array(tess_flag), format="I"
    )  # 16 bit int
    c13 = fits.Column(name="TGLC_flags", array=tglc_flag, format="I")
    table_hdu = fits.BinTableHDU.from_columns(
        [c1, c2, c4, c6, c8, c10, c11, c12, c13]
    )
    table_hdu.header.append(
        ("INHERIT", "T", "inherit the primary header"), end=True
    )
    table_hdu.header.append(
        ("EXTNAME", "LIGHTCURVE", "name of extension"), end=True
    )
    table_hdu.header.append(
        ("EXTVER", 1, "extension version"), end=True  # TODO: version?
    )
    table_hdu.header.append(("TELESCOP", "TESS", "telescope"), end=True)
    table_hdu.header.append(
        ("INSTRUME", "TESS Photometer", "detector type"), end=True
    )
    table_hdu.header.append(
        ("FILTER", "TESS", "the filter used for the observations"), end=True
    )
    table_hdu.header.append(
        (
            "OBJECT",
            source.gaia[index]["DESIGNATION"],
            "string version of Gaia DR3 ID",
        ),
        end=True,
    )
    table_hdu.header.append(
        ("GAIADR3", objid, "integer version of GaiaDR3 designation"), end=True
    )
    table_hdu.header.append(
        ("RADESYS", "ICRS", "reference frame of celestial coordinates"),
        end=True,
    )
    table_hdu.header.append(
        ("RA_OBJ", source.gaia[index]["ra"], "[deg] right ascension, J2000"),
        end=True,
    )
    table_hdu.header.append(
        ("DEC_OBJ", source.gaia[index]["dec"], "[deg] declination, J2000"),
        end=True,
    )
    table_hdu.header.append(
        ("TIMEREF", "SOLARSYSTEM", "barycentric correction applied to times"),
        end=True,
    )
    table_hdu.header.append(
        ("TASSIGN", "SPACECRAFT", "where time is assigned"), end=True
    )
    table_hdu.header.append(
        ("BJDREFI", 2457000, "integer part of BJD reference date"), end=True
    )
    table_hdu.header.append(
        ("BJDREFR", 0.0, "fraction of the day in BJD reference date"), end=True
    )
    table_hdu.header.append(
        ("TIMESYS", "TDB", "TESS Barycentric Dynamical Time"), end=True
    )
    table_hdu.header.append(("TIMEUNIT", "d", "time unit for TIME"), end=True)
    # table_hdu.header.append(('BUNIT', 'e-/s', 'psf_flux unit'), end=True)
    table_hdu.header.append(
        ("TELAPS", t_stop - t_start, "[d] TSTOP-TSTART"), end=True
    )
    table_hdu.header.append(
        ("TSTART", t_start, "[d] observation start time in TBJD"), end=True
    )
    table_hdu.header.append(
        ("TSTOP", t_stop, "[d] observation end time in TBJD"), end=True
    )
    table_hdu.header.append(
        ("MJD_BEG", t_start + 56999.5, "[d] start time in barycentric MJD"),
        end=True,
    )
    table_hdu.header.append(
        ("MJD_END", t_stop + 56999.5, "[d] end time in barycentric MJD"),
        end=True,
    )
    table_hdu.header.append(
        (
            "TIMEDEL",
            (t_stop - t_start) / len(source.time),
            "[d] time resolution of data",
        ),
        end=True,
    )
    table_hdu.header.append(
        ("XPTIME", exposure_time, "[s] exposure time"), end=True
    )
    table_hdu.header.append(
        ("PSF_ERR", psf_err, "[e-/s] PSF flux error"), end=True
    )
    table_hdu.header.append(
        ("APER_ERR", aper_err, "[e-/s] aperture flux error"), end=True
    )
    table_hdu.header.append(
        ("CPSF_ERR", cal_psf_err, "[e-/s] calibrated PSF flux error"), end=True
    )
    table_hdu.header.append(
        ("CAPE_ERR", cal_aper_err, "[e-/s] calibrated aperture flux error"),
        end=True,
    )
    table_hdu.header.append(
        ("NEAREDGE", near_edge, "distance to edges of FFI <= 2"), end=True
    )
    table_hdu.header.append(
        ("LOC_BG", local_bg, "[e-/s] locally modified background"), end=True
    )
    table_hdu.header.append(
        ("COMMENT", "TRUE_BG = hdul[1].data['background'] + LOC_BG"), end=True
    )
    table_hdu.header.append(
        ("WOTAN_WL", 1, "wotan detrending window length"), end=True
    )
    table_hdu.header.append(
        ("WOTAN_MT", "biweight", "wotan detrending method"), end=True
    )
    if isinstance(float, prior):
        table_hdu.header.append(
            ("PRIOR", prior, "prior of field stars"), end=True
        )

    hdul = fits.HDUList([primary_hdu, table_hdu, image_hdu])
    hdul.writeto(
        f"{local_directory}hlsp_tglc_tess_ffi_gaiaid-{objid}-s{source.sector:04d}-cam{source.camera}-ccd{source.ccd}_tess_v2_llc.fits",
        overwrite=True,
    )
    return


def bilinear(x, y, repeat=1):
    """
    A bilinear formula
    np.array([1 - x - y + x * y, x - x * y, y - x * y, x * y] * repeat)
    b, d = array[1]
    a, c = array[0]
    :param x: x
    :param y: y
    :param repeat: side length of epsf
    :return: bilinear interpolation
    """
    return np.array([1 - x - y + x * y, x - x * y, y - x * y, x * y] * repeat)


def get_psf(
    source, factor=2, psf_size=11, edge_compression=1e-4, c=np.array([0, 0, 0])
):
    """
    Generate matrix for PSF fitting
    :param source: tglc.ffi_cut.Source or tglc.ffi_cut.Source_cut, required
    Source or Source_cut object
    :param factor: int, optional
    effective PSF oversampling factor
    :param psf_size: int, optional
    effective PSF side length
    :param edge_compression: float, optional
    parameter for edge compression
    :param c: np.ndarray, optional
    manual modification of Gaia positions in the format of [x, y, theta]
    :return: A, star_info, over_size, x_round, y_round
    A: 2d matrix for least_square
    star_info: star parameters
    over_size: size of oversampled grid of ePSF
    x_round: star horizontal pixel coordinates rounded
    y_round: star vertical pixel coordinates rounded
    """
    # even only
    if factor % 2 != 0:
        raise ValueError("Factor must be even.")
    psf_size = psf_size
    half_size = int((psf_size - 1) / 2)
    over_size = psf_size * factor + 1
    size = source.size  # TODO: must be even?
    flux_ratio = np.array(source.gaia["tess_flux_ratio"])
    # flux_ratio = 0.9998 * flux_ratio + 0.0002
    # x_shift = np.array(source.gaia[f'sector_{source.sector}_x'])
    # y_shift = np.array(source.gaia[f'sector_{source.sector}_y'])

    x_shift = np.array(source.gaia[f"sector_{source.sector}_x"])
    y_shift = np.array(source.gaia[f"sector_{source.sector}_y"])

    # x_shift = (x_ - c[0]) * np.cos(c[2]) - (y_ - c[1]) * np.sin(c[2])
    # y_shift = (x_ - c[0]) * np.sin(c[2]) + (y_ - c[1]) * np.cos(c[2])

    x_round = np.round(x_shift).astype(int)
    y_round = np.round(y_shift).astype(int)

    left = np.maximum(0, x_round - half_size)
    right = np.minimum(size, x_round + half_size) + 1
    down = np.maximum(0, y_round - half_size)
    up = np.minimum(size, y_round + half_size) + 1
    x_residual = x_shift % (1 / factor) * factor
    y_residual = y_shift % (1 / factor) * factor

    x_p = np.arange(size)
    y_p = np.arange(size)
    coord = np.arange(size**2).reshape(size, size)
    xx, yy = np.meshgrid(
        (np.arange(size) - (size - 1) / 2), (np.arange(size) - (size - 1) / 2)
    )

    if isinstance(Source, source):
        bg_dof = 6
        A = np.zeros((size**2, over_size**2 + bg_dof))
        A[:, -1] = np.ones(size**2)
        A[:, -2] = yy.flatten()
        A[:, -3] = xx.flatten()
        A[:, -4] = source.mask.data.flatten()
        A[:, -5] = (source.mask.data * xx).flatten()
        A[:, -6] = (source.mask.data * yy).flatten()
    else:
        bg_dof = 3
        A = np.zeros((size**2, over_size**2 + bg_dof))
        A[:, -1] = np.ones(size**2)
        A[:, -2] = yy.flatten()
        A[:, -3] = xx.flatten()
    star_info = []
    for i in range(len(source.gaia)):
        #     if i == 8:
        #         continue
        x_psf = factor * (x_p[left[i] : right[i]] - x_round[i] + half_size) + (
            x_shift[i] % 1
        ) // (1 / factor)
        y_psf = factor * (y_p[down[i] : up[i]] - y_round[i] + half_size) + (
            y_shift[i] % 1
        ) // (1 / factor)
        x_psf, y_psf = np.meshgrid(x_psf, y_psf)  # super slow here
        a = np.array(x_psf + y_psf * over_size, dtype=np.int64).flatten()
        index = coord[down[i] : up[i], left[i] : right[i]]
        A[
            np.repeat(index, 4),
            np.array([a, a + 1, a + over_size, a + over_size + 1]).flatten(
                order="F"
            ),
        ] += flux_ratio[i] * bilinear(
            x_residual[i], y_residual[i], repeat=len(a)
        )
        # star_info.append(
        #     (np.repeat(index, 4), np.array([a, a + 1, a + over_size, a + over_size + 1]).flatten(order='F'),
        #      flux_ratio[i] * bilinear(x_residual[i], y_residual[i], repeat=len(a))))
        star_info.append(
            (index, a, flux_ratio[i] * bilinear(x_residual[i], y_residual[i]))
        )
    coord_ = np.arange(-psf_size * factor / 2 + 1, psf_size * factor / 2 + 2)
    x_coord, y_coord = np.meshgrid(coord_, coord_)
    variance = psf_size
    dist = (
        1 - np.exp(-0.5 * (x_coord**4 + y_coord**4) / variance**4)
    ) * edge_compression  # 1e-3
    A_mod = np.diag(dist.flatten())
    A_mod = np.concatenate(
        (A_mod, (np.zeros((over_size**2, bg_dof)))), axis=-1
    )
    A = np.append(A, A_mod, axis=0)
    return A, star_info, over_size, x_round, y_round


def fit_lc(
    A,
    source,
    star_info=None,
    x=0.0,
    y=0.0,
    star_num=0,
    factor=2,
    psf_size=11,
    e_psf=None,
    near_edge=False,
):
    """
    Produce matrix for least_square fitting without a certain target
    :param A: np.ndarray, required
    2d matrix for least_square
    :param source: tglc.ffi_cut.Source or tglc.ffi_cut.Source_cut, required
    Source or Source_cut object
    :param star_info: np.ndarray, required
    star parameters
    :param x: float, required
    target horizontal pixel coordinate
    :param y: float, required
    target vertical pixel coordinate
    :param star_num: int, required
    target star index
    :param factor: int, optional
    effective PSF oversampling factor
    :param psf_size: int, optional
    effective PSF side length
    :param e_psf: np.ndarray, required
    effective PSF as a 3d array as a timeseries
    :param near_edge: boolean, required
    whether the star is 2 pixels or closer to the edge of a CCD
    :return: aperture lightcurve, PSF lightcurve, vertical pixel coord, horizontal pixel coord, portion of light in aperture
    """
    over_size = psf_size * factor + 1
    a = star_info[star_num][1]
    star_info_num = (
        np.repeat(star_info[star_num][0], 4),
        np.array([a, a + 1, a + over_size, a + over_size + 1]).flatten(
            order="F"
        ),
        np.tile(star_info[star_num][2], len(a)),
    )
    size = source.size  # TODO: must be even?
    # star_position = int(x + source.size * y - 5 * size - 5)
    # aper_lc
    cut_size = 5
    in_frame = np.where(np.invert(np.isnan(source.flux[0])))
    left = np.maximum(np.min(in_frame[1]), x - cut_size // 2)
    right = np.minimum(np.max(in_frame[1]), x + cut_size // 2 + 1)
    down = np.maximum(np.min(in_frame[0]), y - cut_size // 2)
    up = np.minimum(np.max(in_frame[0]), y + cut_size // 2 + 1)
    coord = np.arange(size**2).reshape(size, size)
    index = np.array(coord[down:up, left:right]).flatten()
    A_cut = np.zeros((len(index), np.shape(A)[1]))
    A_target = np.zeros((len(index), np.shape(A)[1]))
    for i in range(len(index)):
        A_ = np.zeros(np.shape(A)[-1])
        star_pos = np.where(star_info_num[0] == index[i])[0]
        A_[star_info_num[1][star_pos]] = star_info_num[2][star_pos]
        A_target[i] = A_
        A_cut[i] = A[index[i], :] - A_
    aperture = np.zeros((len(source.time), len(index)))
    for j in range(len(source.time)):
        aperture[j] = np.array(
            source.flux[j][down:up, left:right]
        ).flatten() - np.dot(A_cut, e_psf[j])
    aperture = aperture.reshape((len(source.time), up - down, right - left))
    target_5x5 = np.dot(A_target, np.nanmedian(e_psf, axis=0)).reshape(
        cut_size, cut_size
    )
    field_stars_5x5 = np.dot(A_cut, np.nanmedian(e_psf, axis=0)).reshape(
        cut_size, cut_size
    )
    if target_5x5.shape != (cut_size, cut_size):
        # Pad with nans to get to 5x5 shape
        # Pad amount in a direction is (expected_num_pix) - (actual_num_pix)
        pad_left = (cut_size // 2) - (x - left)
        pad_right = (cut_size // 2 + 1) - (right - x)
        pad_down = (cut_size // 2) - (y - down)
        pad_up = (cut_size // 2 + 1) - (up - y)
        target_5x5 = np.pad(
            target_5x5,
            [(pad_down, pad_up), (pad_left, pad_right)],
            constant_values=np.nan,
        )
        field_stars_5x5 = np.pad(
            field_stars_5x5,
            [(pad_down, pad_up), (pad_left, pad_right)],
            constant_values=np.nan,
        )

    # psf_lc
    over_size = psf_size * factor + 1
    if near_edge:  # TODO: near_edge
        psf_lc = np.zeros(len(source.time))
        psf_lc[:] = np.NaN
        e_psf_1d = np.nanmedian(e_psf[:, : over_size**2], axis=0).reshape(
            over_size, over_size
        )
        portion = (
            (36 / 49) * np.nansum(e_psf_1d[8:15, 8:15]) / np.nansum(e_psf_1d)
        )  # only valid for factor = 2
        return (
            aperture,
            psf_lc,
            y - down,
            x - left,
            portion,
            target_5x5,
            field_stars_5x5,
        )
    left_ = left - x + 5
    right_ = right - x + 5
    down_ = down - y + 5
    up_ = up - y + 5

    left_11 = np.maximum(-x + 5, 0)
    right_11 = np.minimum(size - x + 5, 11)
    down_11 = np.maximum(-y + 5, 0)
    up_11 = np.minimum(size - y + 5, 11)
    coord = np.arange(psf_size**2).reshape(psf_size, psf_size)
    index = coord[down_11:up_11, left_11:right_11]
    if isinstance(Source, source):
        bg_dof = 6
    else:
        bg_dof = 3
    A = np.zeros((psf_size**2, over_size**2 + bg_dof))
    A[np.repeat(index, 4), star_info_num[1]] = star_info_num[2]
    psf_shape = np.dot(e_psf, A.T).reshape(
        len(source.time), psf_size, psf_size
    )
    psf_sim = psf_shape[:, down_:up_, left_:right_]
    # psf_sim = np.transpose(psf_shape[:, down_:up_, left_: right_], (0, 2, 1))

    psf_lc = np.zeros(len(source.time))
    A_ = np.zeros((cut_size**2, 4))
    xx, yy = np.meshgrid(
        (np.arange(cut_size) - (cut_size - 1) / 2),
        (np.arange(cut_size) - (cut_size - 1) / 2),
    )
    A_[:, -1] = np.ones(cut_size**2)
    A_[:, -2] = yy.flatten()
    A_[:, -3] = xx.flatten()
    edge_pixel = np.array(
        [0, 1, 2, 3, 4, 5, 9, 10, 14, 15, 19, 20, 21, 22, 23, 24]
    )
    # edge_pixel = np.array([0, 1, 2, 3, 4, 5, 6,
    #                        7, 8, 9, 10, 11, 12, 13,
    #                        14, 15, 19, 20,
    #                        21, 22, 26, 27,
    #                        28, 29, 33, 34,
    #                        35, 36, 37, 38, 39, 40, 41,
    #                        42, 43, 44, 45, 46, 47, 48])
    med_aperture = np.median(aperture, axis=0).flatten()
    outliers = np.abs(
        med_aperture[edge_pixel] - np.nanmedian(med_aperture[edge_pixel])
    ) > 1 * np.std(med_aperture[edge_pixel])
    epsf_sum = np.sum(np.nanmedian(psf_shape, axis=0))
    for j in range(len(source.time)):
        if np.isnan(psf_sim[j, :, :]).any():
            psf_lc[j] = np.nan
        else:
            aper_flat = aperture[j, :, :].flatten()
            A_[:, 0] = psf_sim[j, :, :].flatten() / epsf_sum
            a = np.delete(A_, edge_pixel[outliers], 0)
            aper_flat = np.delete(aper_flat, edge_pixel[outliers])
            psf_lc[j] = np.linalg.lstsq(a, aper_flat)[0][0]
    portion = np.nansum(psf_shape[:, 4:7, 4:7]) / np.nansum(psf_shape)
    # print(np.nansum(psf_shape[:, 5, 5]) / np.nansum(psf_shape))
    # np.save(f'toi-5344_psf_{source.sector}.npy', psf_shape)
    return (
        aperture,
        psf_lc,
        y - down,
        x - left,
        portion,
        target_5x5,
        field_stars_5x5,
    )


def fit_lc_float_field(
    A,
    source,
    star_info=None,
    x=np.array([]),
    y=np.array([]),
    star_num=0,
    factor=2,
    psf_size=11,
    e_psf=None,
    near_edge=False,
    prior=0.001,
):
    """
    Produce matrix for least_square fitting without a certain target
    :param A: np.ndarray, required
    2d matrix for least_square
    :param source: tglc.ffi_cut.Source or tglc.ffi_cut.Source_cut, required
    Source or Source_cut object
    :param star_info: np.ndarray, required
    star parameters
    :param x: float, required
    target horizontal pixel coordinate
    :param y: float, required
    target vertical pixel coordinate
    :param star_num: int, required
    target star index
    :param factor: int, optional
    effective PSF oversampling factor
    :param psf_size: int, optional
    effective PSF side length
    :param e_psf: np.ndarray, required
    effective PSF as a 3d array as a timeseries
    :param near_edge: boolean, required
    whether the star is 2 pixels or closer to the edge of a CCD
    :return: aperture lightcurve, PSF lightcurve, vertical pixel coord, horizontal pixel coord, portion of light in aperture
    """
    over_size = psf_size * factor + 1
    a = star_info[star_num][1]
    star_info_num = (
        np.repeat(star_info[star_num][0], 4),
        np.array([a, a + 1, a + over_size, a + over_size + 1]).flatten(
            order="F"
        ),
        np.tile(star_info[star_num][2], len(a)),
    )
    size = source.size  # TODO: must be even?
    # star_position = int(x + source.size * y - 5 * size - 5)
    # aper_lc
    cut_size = 5
    in_frame = np.where(np.invert(np.isnan(source.flux[0])))
    left = np.maximum(np.min(in_frame[1]), x[star_num] - cut_size // 2)
    right = np.minimum(np.max(in_frame[1]), x[star_num] + cut_size // 2 + 1)
    down = np.maximum(np.min(in_frame[0]), y[star_num] - cut_size // 2)
    up = np.minimum(np.max(in_frame[0]), y[star_num] + cut_size // 2 + 1)
    coord = np.arange(size**2).reshape(size, size)
    index = np.array(coord[down:up, left:right]).flatten()
    A_cut = np.zeros((len(index), np.shape(A)[1]))
    for i in range(len(index)):
        A_ = np.zeros(np.shape(A)[-1])
        star_pos = np.where(star_info_num[0] == index[i])[0]
        A_[star_info_num[1][star_pos]] = star_info_num[2][star_pos]
        A_cut[i] = A[index[i], :] - A_
    aperture = np.zeros((len(source.time), len(index)))
    for j in range(len(source.time)):
        aperture[j] = np.array(
            source.flux[j][down:up, left:right]
        ).flatten() - np.dot(A_cut, e_psf[j])
    aperture = aperture.reshape((len(source.time), up - down, right - left))

    # psf_lc
    over_size = psf_size * factor + 1
    if near_edge:  # TODO: near_edge
        psf_lc = np.zeros(len(source.time))
        psf_lc[:] = np.NaN
        e_psf_1d = np.nanmedian(e_psf[:, : over_size**2], axis=0).reshape(
            over_size, over_size
        )
        portion = (
            (36 / 49) * np.nansum(e_psf_1d[8:15, 8:15]) / np.nansum(e_psf_1d)
        )  # only valid for factor = 2
        return (
            aperture,
            psf_lc,
            y[star_num] - down,
            x[star_num] - left,
            portion,
        )
    # left_ = left - x[star_num] + 5
    # right_ = right - x[star_num] + 5
    # down_ = down - y[star_num] + 5
    # up_ = up - y[star_num] + 5
    if isinstance(Source, source):
        bg_dof = 6
    else:
        bg_dof = 3
    field_star_num = []
    for j in range(len(source.gaia)):
        if np.abs(x[j] - x[star_num]) < 5 and np.abs(y[j] - y[star_num]) < 5:
            field_star_num.append(j)

    psf_lc = np.zeros(len(source.time))
    A_ = np.zeros((cut_size**2 + len(field_star_num), len(field_star_num) + 3))
    xx, yy = np.meshgrid(
        (np.arange(cut_size) - (cut_size - 1) / 2),
        (np.arange(cut_size) - (cut_size - 1) / 2),
    )
    A_[: (cut_size**2), -1] = np.ones(cut_size**2)
    A_[: (cut_size**2), -2] = yy.flatten()
    A_[: (cut_size**2), -3] = xx.flatten()
    psf_sim = np.zeros(
        (len(source.time), 11**2 + len(field_star_num), len(field_star_num))
    )
    for j, star in enumerate(field_star_num):
        a = star_info[star][1]
        star_info_star = (
            np.repeat(star_info[star][0], 4),
            np.array([a, a + 1, a + over_size, a + over_size + 1]).flatten(
                order="F"
            ),
            np.tile(star_info[star][2], len(a)),
        )
        delta_x = x[star_num] - x[star]
        delta_y = y[star_num] - y[star]
        # for psf_sim
        left_shift = np.maximum(delta_x, 0)
        right_shift = np.minimum(11 + delta_x, 11)
        down_shift = np.maximum(delta_y, 0)
        up_shift = np.minimum(11 + delta_y, 11)
        # for psf_shape
        left_shift_ = np.maximum(-delta_x, 0)
        right_shift_ = np.minimum(11 - delta_x, 11)
        down_shift_ = np.maximum(-delta_y, 0)
        up_shift_ = np.minimum(11 - delta_y, 11)

        left_11 = np.maximum(-x[star] + 5, 0)
        right_11 = np.minimum(size - x[star] + 5, 11)
        down_11 = np.maximum(-y[star] + 5, 0)
        up_11 = np.minimum(size - y[star] + 5, 11)

        coord = np.arange(psf_size**2).reshape(psf_size, psf_size)
        index = coord[down_11:up_11, left_11:right_11]
        A = np.zeros((psf_size**2, over_size**2 + bg_dof))
        A[np.repeat(index, 4), star_info_star[1]] = star_info_star[2]
        psf_shape = np.dot(e_psf, A.T).reshape(
            len(source.time), psf_size, psf_size
        )
        epsf_sum = np.sum(np.nanmedian(psf_shape, axis=0))
        psf_sim_index = coord[
            down_shift:up_shift, left_shift:right_shift
        ].flatten()
        psf_sim[:, psf_sim_index, j] = (
            psf_shape[
                :, down_shift_:up_shift_, left_shift_:right_shift_
            ].reshape(len(source.time), -1)
            / epsf_sum
        )
        if star != star_num:
            psf_sim[:, 11**2 + j, j] = np.ones(len(source.time)) / (
                prior
                * 1.5e4
                * 10 ** ((10 - source.gaia[star]["tess_mag"]) / 2.5)
            )
        else:
            portion = np.nansum(psf_shape[:, 4:7, 4:7]) / np.nansum(psf_shape)

    star_index = np.where(np.array(field_star_num) == star_num)[0]
    field_star = (
        psf_sim[0, np.arange(11**2).reshape(11, 11)[3:8, 3:8], :].reshape(
            cut_size**2, len(field_star_num)
        )
        * source.gaia["tess_flux_ratio"][field_star_num]
    )
    field_star[:, star_index] = 0
    for j in range(len(source.time)):
        if np.isnan(psf_sim[j, :, :]).any():
            psf_lc[j] = np.nan
        else:
            aper_flat = aperture[j, :, :].flatten()
            aper_flat = np.append(
                aper_flat, np.zeros(len(field_star_num) - 1)
            )  # / prior
            aper_flat[cut_size**2 + star_index] = 0
            postcards = psf_sim[
                j, np.arange(11**2).reshape(11, 11)[3:8, 3:8], :
            ].reshape(cut_size**2, len(field_star_num))
            A_[: cut_size**2, : len(field_star_num)] = postcards
            field_star = (
                postcards * source.gaia["tess_flux_ratio"][field_star_num]
            )
            field_star[:, star_index] = 0
            # A_[:(cut_size ** 2), -4] = np.sum(field_star, axis=1)
            A_[cut_size**2 :, : len(field_star_num)] = psf_sim[
                j, 11**2 :, :
            ].reshape(len(field_star_num), len(field_star_num))
            a = np.delete(A_, cut_size**2 + star_index, 0)
            psf_lc[j] = np.linalg.lstsq(a, aper_flat)[0][star_index]
    return aperture, psf_lc, y[star_num] - down, x[star_num] - left, portion


def fit_psf(A, source, over_size, power=0.8, time=0):
    """
    fit_psf using least_square (improved performance by changing to np.linalg.solve)
    :param A: np.ndarray, required
    2d matrix for least_square
    :param source: tglc.ffi_cut.Source or tglc.ffi_cut.Source_cut, required
    Source or Source_cut object
    :param over_size: int, required
    size of oversampled grid of ePSF
    :param power: float, optional
    power for weighting bright stars' contribution to the fit. 1 means same contribution from all stars,
    <1 means emphasizing dimmer stars
    :param time: int, required
    time index of this ePSF fit
    :return: fit result
    """
    saturated_index = source.mask.mask.flatten()
    flux = source.flux[time].flatten()
    saturated_index[flux < 0.8 * np.nanmedian(flux)] = True

    b = np.delete(flux, saturated_index)
    scaler = np.abs(np.delete(flux, saturated_index)) ** power
    b = np.append(b, np.zeros(over_size**2))
    scaler = np.append(scaler, np.ones(over_size**2))

    # fit = np.linalg.lstsq(A / scaler[:, np.newaxis], b / scaler, rcond=None)[0]
    a = np.delete(A, np.where(saturated_index), 0) / scaler[:, np.newaxis]
    b = b / scaler
    alpha = np.dot(a.T, a)
    beta = np.dot(a.T, b)
    try:
        fit = np.linalg.solve(alpha, beta)
    except np.linalg.LinAlgError:
        fit = np.full(np.shape(a)[1], np.nan)
        # fit = np.linalg.lstsq(a, b, rcond=None)[0]
    return fit


def bg_mod(
    source,
    q=None,
    aper_lc=None,
    psf_lc=None,
    portion=None,
    star_num=0,
    near_edge=False,
):
    """
    background modification
    :param source: tglc.ffi_cut.Source or tglc.ffi_cut.Source_cut, required
    Source or Source_cut object
    :param q: list, optional
    list of booleans that filter the data points
    :param aper_lc: np.ndarray, required
    aperture light curve
    :param psf_lc: np.ndarray, required
    PSF light curve
    :param portion: float, required
    portion of light in aperture
    :param star_num: int, required,
    star index
    :param near_edge: boolean, required
    whether the star is 2 pixels or closer to the edge of a CCD
    :return: local background, modified aperture light curve, modified PSF light curve
    """
    bar = 15000 * 10 ** ((source.gaia["tess_mag"][star_num] - 10) / -2.5)
    # print(bar)
    # med_epsf = np.nanmedian(e_psf[:, :23 ** 2].reshape(len(source.time), 23, 23), axis=0)
    # centroid_to_aper_ratio = 4/9 * np.sum(med_epsf[10:13, 10:13]) / np.sum(med_epsf)
    # centroid_to_aper_ratio = np.nanmedian(ratio)
    # flux_bar = aperture_bar * centroid_to_aper_ratio
    # lightcurve = lightcurve + (flux_bar - np.nanmedian(lightcurve[q]))
    aperture_bar = bar * portion
    # print(bar)
    local_bg = np.nanmedian(aper_lc[q]) - aperture_bar
    if np.isnan(local_bg):
        local_bg = 0
    aper_lc = aper_lc - local_bg
    psf_bar = bar
    local_bg = np.nanmedian(psf_lc[q]) - psf_bar
    if np.isnan(local_bg):
        local_bg = 0
    psf_lc = psf_lc - local_bg
    negative_arg_aper = np.where(aper_lc <= 0)  # Negative frames
    aper_lc[negative_arg_aper] = np.nan
    negative_arg_psf = np.where(psf_lc <= 0)
    psf_lc[negative_arg_psf] = np.nan
    # removes very large outliers to prevent wotan to freeze
    cal_aper_lc = aper_lc / np.nanmedian(aper_lc)
    cal_aper_lc[np.where(cal_aper_lc > 100)] = np.nan
    if np.isnan(cal_aper_lc).all():
        print(
            "Calibrated aperture flux are not accessible or processed incorrectly. "
        )
    else:
        _, trend = flatten(
            source.time,
            cal_aper_lc - np.nanmin(cal_aper_lc) + 1000,
            window_length=1,
            method="biweight",
            return_trend=True,
        )
        cal_aper_lc = (
            cal_aper_lc - np.nanmin(cal_aper_lc) + 1000 - trend
        ) / np.nanmedian(cal_aper_lc) + 1
        # cal_aper_lc = flatten(source.time, cal_aper_lc, window_length=1, method='biweight',
        #                       return_trend=False)
    if near_edge:
        cal_psf_lc = psf_lc
        return local_bg, aper_lc, psf_lc, cal_aper_lc, cal_psf_lc
    else:
        cal_psf_lc = psf_lc / np.nanmedian(psf_lc)
        cal_psf_lc[np.where(cal_psf_lc > 100)] = np.nan
        if np.isnan(cal_psf_lc).all():
            print(
                "Calibrated PSF flux are not accessible or processed incorrectly. "
            )
        else:
            _, trend = flatten(
                source.time,
                cal_psf_lc - np.nanmin(cal_psf_lc) + 1000,
                window_length=1,
                method="biweight",
                return_trend=True,
            )
            cal_psf_lc = (
                cal_psf_lc - np.nanmin(cal_psf_lc) + 1000 - trend
            ) / np.nanmedian(cal_psf_lc) + 1
            # cal_psf_lc = flatten(source.time, cal_psf_lc, window_length=1, method='biweight',
            #                      return_trend=False)
    # aper_mad = 1.4826 * np.nanmedian(np.abs(cal_aper_lc - 1))
    # if aper_mad > 0.02:
    #     psf_mad = 1.4826 * np.nanmedian(np.abs(cal_psf_lc - 1))
    #     cal_psf_lc /= psf_mad / aper_mad
    #     cal_psf_lc += 1 - np.median(cal_psf_lc)
    return local_bg, aper_lc, psf_lc, cal_aper_lc, cal_psf_lc


def epsf(
    source,
    psf_size=11,
    factor=2,
    local_directory="",
    target=None,
    cut_x=0,
    cut_y=0,
    sector=0,
    limit_mag=16,
    edge_compression=1e-4,
    power=1.4,
    name=None,
    save_aper=False,
    no_progress_bar=False,
    prior=None,
):
    """
    User function that unites all necessary steps
    :param source: TGLC.ffi_cut.Source or TGLC.ffi_cut.Source_cut, required
    Source or Source_cut object
    :param factor: int, optional
    effective PSF oversampling factor
    :param local_directory: string, required
    output directory
    :param sector: int, required
    TESS sector number
    :param limit_mag: int, required
    upper limiting magnitude of the lightcurves that are outputted.
    :param edge_compression: float, optional
    parameter for edge compression
    :param power: float, optional
    power for weighting bright stars' contribution to the fit. 1 means same contribution from all stars,
    <1 means emphasizing dimmer stars
    :return:
    """
    if target is None:
        target = f"{cut_x:02d}_{cut_y:02d}"
    A, star_info, over_size, x_round, y_round = get_psf(
        source,
        psf_size=psf_size,
        factor=factor,
        edge_compression=edge_compression,
    )
    lc_directory = f"{local_directory}lc/{source.camera}-{source.ccd}/"
    epsf_loc = f"{local_directory}epsf/{source.camera}-{source.ccd}/epsf_{target}_sector_{sector}_{source.camera}-{source.ccd}.npy"
    if isinstance(Source_cut, source):
        bg_dof = 3
        lc_directory = f"{local_directory}lc/"
        epsf_loc = f"{local_directory}epsf/epsf_{target}_sector_{sector}.npy"
    else:
        bg_dof = 6
    os.makedirs(lc_directory, exist_ok=True)
    # sim_image = np.dot(A[:source.size ** 2, :], fit_psf(A, source, over_size, power=power, time=2817).T)
    # # residual = np.abs(source.flux[2817].flatten() - sim_image)
    # residual = source.flux[2817].flatten() - sim_image
    # return residual.reshape((source.size, source.size))

    epsf_exists = exists(epsf_loc)
    if epsf_exists:
        e_psf = np.load(epsf_loc)
        print(f"Loaded ePSF {target} from directory. ")
    else:
        e_psf = np.zeros((len(source.time), over_size**2 + bg_dof))
        for i in trange(
            len(source.time), desc="Fitting ePSF", disable=no_progress_bar
        ):
            e_psf[i] = fit_psf(A, source, over_size, power=power, time=i)
        if np.isnan(e_psf).any():
            warnings.warn(
                f"TESS FFI cut includes Nan values. Please shift the center of the cutout to remove Nan near edge. Target: {target}"
            )
        np.save(epsf_loc, e_psf)
    # contamination_8 = np.dot(A[:source.size ** 2, :], e_psf[0].T)
    # np.save('/mnt/c/users/tehan/desktop/7654/contamination_8_.npy', contamination_8)
    # TODO: quality use which background?
    background = np.dot(A[: source.size**2, -bg_dof:], e_psf[:, -bg_dof:].T)
    quality_raw = np.zeros(len(source.time), dtype=np.int16)
    sigma = 1.4826 * np.nanmedian(
        np.abs(e_psf[:, -1] - np.nanmedian(e_psf[:, -1]))
    )
    quality_raw[
        abs(e_psf[:, -1] - np.nanmedian(e_psf[:, -1])) >= 3 * sigma
    ] += 1
    index_1 = np.where(np.array(source.quality) == 0)[0]
    index_2 = np.where(quality_raw == 0)[0]
    index = np.intersect1d(index_1, index_2)
    if isinstance(Source_cut, source):
        in_frame = np.where(np.invert(np.isnan(source.flux[0])))
        x_left = np.min(in_frame[1]) - 0.5
        x_right = source.size - np.max(in_frame[1]) + 0.5
        y_left = np.min(in_frame[0]) - 0.5
        y_right = source.size - np.max(in_frame[0]) + 0.5
    else:
        x_left = 1.5 if cut_x != 0 else -0.5
        x_right = 2.5 if cut_x != 13 else 0.5
        y_left = 1.5 if cut_y != 0 else -0.5
        y_right = 2.5 if cut_y != 13 else 0.5

    num_stars = np.array(source.gaia["tess_mag"]).searchsorted(
        limit_mag, "right"
    )
    x_aperture = source.gaia[f"sector_{source.sector}_x"] - np.maximum(
        0, x_round - 2
    )
    y_aperture = source.gaia[f"sector_{source.sector}_y"] - np.maximum(
        0, y_round - 2
    )

    start = 0
    end = num_stars
    if name is not None:
        try:
            start = int(
                np.where(
                    source.gaia["DESIGNATION"]
                    == "Gaia DR3 "
                    + str(
                        source.tic["dr3_source_id"][
                            np.where(source.tic["TIC"] == name)
                        ][0]
                    )
                )[0][0]
            )
            end = start + 1
        except Exception as e:
            print(e)
            try:
                start = int(np.where(source.gaia["DESIGNATION"] == name)[0][0])
                end = start + 1
            except IndexError:
                print(
                    f"Target not found in the requested sector (Sector {sector}). This can be caused by a lack of Gaia "
                    "ID or an incomplete TESS to Gaia crossmatch table. Please check whether the output light curve Gaia"
                    " DR3 ID agrees with your target."
                )
                start = 0
                end = 1
    for i in trange(start, end, desc="Fitting lc", disable=no_progress_bar):
        if (
            x_left <= x_round[i] < source.size - x_right
            and y_left <= y_round[i] < source.size - y_right
        ):
            if isinstance(Source, source):
                x_left = 1.5
                x_right = 2.5
                y_left = 1.5
                y_right = 2.5
            if x_left + 2 <= x_round[i] < source.size - (
                x_right + 2
            ) and y_left + 2 <= y_round[i] < source.size - (y_right + 2):
                near_edge = False
            else:
                near_edge = True

            if isinstance(float, prior) or isinstance(np.float64, prior):
                aperture, psf_lc, star_y, star_x, portion = fit_lc_float_field(
                    A,
                    source,
                    star_info=star_info,
                    x=x_round,
                    y=y_round,
                    star_num=i,
                    e_psf=e_psf,
                    near_edge=near_edge,
                    prior=prior,
                )
            else:
                (
                    aperture,
                    psf_lc,
                    star_y,
                    star_x,
                    portion,
                    target_5x5,
                    field_stars_5x5,
                ) = fit_lc(
                    A,
                    source,
                    star_info=star_info,
                    x=x_round[i],
                    y=y_round[i],
                    star_num=i,
                    e_psf=e_psf,
                    near_edge=near_edge,
                )
            aper_lc = np.sum(
                aperture[
                    :,
                    max(0, star_y - 1) : min(5, star_y + 2),
                    max(0, star_x - 1) : min(5, star_x + 2),
                ],
                axis=(1, 2),
            )
            if source.sector < 27:  # primary
                exposure_time = 1800
            elif source.sector < 56:  # first extended
                exposure_time = 600
            else:  # second extended
                exposure_time = 200
            saturated_arg_aper = np.where(
                aper_lc > 1e5 * 9 * 2e5 / exposure_time
            )  # saturation is 2e5 e-
            aper_lc[saturated_arg_aper] = np.nan
            saturated_arg_psf = np.where(
                psf_lc > 1e5 * 9 * 2e5 / exposure_time
            )
            psf_lc[saturated_arg_psf] = np.nan

            local_bg, aper_lc, psf_lc, cal_aper_lc, cal_psf_lc = bg_mod(
                source,
                q=index,
                portion=portion,
                psf_lc=psf_lc,
                aper_lc=aper_lc,
                near_edge=near_edge,
                star_num=i,
            )
            background_ = background[x_round[i] + source.size * y_round[i], :]
            quality = np.zeros(len(source.time), dtype=np.int16)
            sigma = 1.4826 * np.nanmedian(
                np.abs(background_ - np.nanmedian(background_))
            )
            quality[
                abs(background_ - np.nanmedian(background_)) >= 5 * sigma
            ] += 1
            # if cut_x == 7:
            #     lc_directory = f'{local_directory}lc/{source.camera}-{source.ccd}_extra/'
            #     os.makedirs(lc_directory, exist_ok=True)
            if np.isnan(aper_lc).all():
                continue
            else:
                if isinstance(Source, source):
                    # if cut_x >= 7:
                    #     lc_directory = f'{local_directory}lc/{source.camera}-{source.ccd}_extra/'
                    lc_output(
                        source,
                        local_directory=lc_directory,
                        index=i,
                        tess_flag=source.quality,
                        cut_x=cut_x,
                        cut_y=cut_y,
                        cadence=source.cadence,
                        aperture=aperture.astype(np.float32),
                        star_y=y_round[i],
                        star_x=x_round[i],
                        tglc_flag=quality,
                        bg=background_,
                        time=source.time,
                        psf_lc=psf_lc,
                        cal_psf_lc=cal_psf_lc,
                        aper_lc=aper_lc,
                        cal_aper_lc=cal_aper_lc,
                        local_bg=local_bg,
                        x_aperture=x_aperture[i],
                        y_aperture=y_aperture[i],
                        near_edge=near_edge,
                        save_aper=save_aper,
                        portion=portion,
                        prior=prior,
                        transient=source.transient,
                        target_5x5=target_5x5,
                        field_stars_5x5=field_stars_5x5,
                    )
                else:
                    lc_output(
                        source,
                        local_directory=lc_directory,
                        index=i,
                        tess_flag=source.quality,
                        cut_x=cut_x,
                        cut_y=cut_y,
                        cadence=source.cadence,
                        aperture=aperture.astype(np.float32),
                        star_y=y_round[i],
                        star_x=x_round[i],
                        tglc_flag=quality,
                        bg=background_,
                        time=source.time,
                        psf_lc=psf_lc,
                        cal_psf_lc=cal_psf_lc,
                        aper_lc=aper_lc,
                        cal_aper_lc=cal_aper_lc,
                        local_bg=local_bg,
                        x_aperture=x_aperture[i],
                        y_aperture=y_aperture[i],
                        near_edge=near_edge,
                        save_aper=save_aper,
                        portion=portion,
                        prior=prior,
                        transient=source.transient,
                        target_5x5=target_5x5,
                        field_stars_5x5=field_stars_5x5,
                    )


if __name__ == "__main__":
    target = "TIC 264468702"
    size = 90
    limit_mag = 16
    sector = None
    prior = None
    transient = None
    if not exists():
        os.makedirs("~/.tglc")
    data_path = f"~/.tglc/{target.replace(' ','')}"
    os.makedirs(data_path + "logs/", exist_ok=True)
    os.makedirs(data_path + "lc/", exist_ok=True)
    os.makedirs(data_path + "epsf/", exist_ok=True)
    os.makedirs(data_path + "source/", exist_ok=True)

    target_ = Catalogs.query_object(
        target, radius=42 * 0.707 / 3600, catalog="Gaia", version=2
    )
    ra = target_[0]["ra"]
    dec = target_[0]["dec"]
    coord = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame="icrs")
    sector_table = Tesscut.get_sectors(coordinates=coord)
    sector = sector_table["sector"][0]

    catalogdata = Catalogs.query_object(
        str(target), radius=0.02, catalog="TIC"
    )
    if target[0:3] == "TIC":
        name = int(target[4:])
    elif transient is not None:
        name = transient[0]
    else:
        name = int(np.array(catalogdata["ID"])[0])

    source = ffi_cut(
        target=target,
        size=size,
        local_directory=data_path,
        sector=sector,
        limit_mag=limit_mag,
        transient=transient,
    )  # sector
    source.select_sector(sector=source.sector_table["sector"][0])
    epsf(
        source,
        factor=2,
        sector=source.sector,
        target=target,
        local_directory=data_path,
        name=name,
        limit_mag=limit_mag,
        save_aper=False,
        prior=prior,
    )
