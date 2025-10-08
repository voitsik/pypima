"""Some useful tools for PIMA data analysis."""

import io
from collections.abc import Iterable

import pycurl
from astropy.io import fits


def fits_has_calib_tables(file_name: str) -> bool:
    """Check FITS-IDI file for existence of calibration tables."""
    has_gain_curves = False
    has_sys_temperatures = False

    with fits.open(file_name, mode="denywrite") as hdulist:
        for hdu in hdulist:
            if isinstance(hdu, fits.BinTableHDU):
                if hdu.name == "SYSTEM_TEMPERATURE":
                    has_sys_temperatures = True
                elif hdu.name == "GAIN_CURVE":
                    has_gain_curves = True

    return has_gain_curves and has_sys_temperatures


def fits_has_orbiting_antenna(file_name_list: Iterable[str]) -> bool:
    """Check FITS-IDI files for existence of orbiting antenna."""
    has_orbit_table = False
    has_ra_in_array_table = False
    has_ra_in_antenna_table = False

    for file_name in file_name_list:
        with fits.open(file_name, mode="denywrite") as hdulist:
            for hdu in hdulist:
                if isinstance(hdu, fits.BinTableHDU):
                    if hdu.name == "ARRAY_GEOMETRY" and "RA" in hdu.data["ANNAME"]:
                        has_ra_in_array_table = True
                    elif hdu.name == "ANTENNA" and "RA" in hdu.data["ANNAME"]:
                        has_ra_in_antenna_table = True
                    elif hdu.name == "SPACECRAFT_ORBIT":
                        has_orbit_table = True

    return has_orbit_table or has_ra_in_array_table or has_ra_in_antenna_table


def download_file(
    url: str,
    buffer: io.BufferedIOBase,
    max_retries: int = 0,
    ftp_user: str | None = None,
):
    """
    Download data from `url` and write it to `buffer` using pycurl.

    Parameters
    ----------
    url : str
        URL
    buffer : object
        Object with `write` function. For inctance, BytesIO or file descriptor.
    max_retries : int
        Number of attempts to download.

    """
    done = False
    retries = 0

    curl = pycurl.Curl()
    curl.setopt(pycurl.URL, url)
    curl.setopt(pycurl.CONNECTTIMEOUT, 30)
    curl.setopt(pycurl.LOW_SPEED_LIMIT, 10000)
    curl.setopt(pycurl.LOW_SPEED_TIME, 60)
    curl.setopt(pycurl.NETRC, 1)
    if ftp_user:
        curl.setopt(pycurl.USERNAME, ftp_user)
    curl.setopt(pycurl.WRITEDATA, buffer)

    while not done:
        try:
            curl.perform()
        except pycurl.error as err:
            errno, errstr = err.args

            if errno == pycurl.E_OPERATION_TIMEDOUT and retries < max_retries:
                retries += 1

                # Try to continue data transfer
                curl.setopt(pycurl.RESUME_FROM_LARGE, buffer.tell())
            else:
                raise
        else:
            done = True

    curl.close()
