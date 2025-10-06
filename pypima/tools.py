"""Some useful tools for PIMA data analysis."""

import io

import pycurl
from astropy.io import fits


def check_fits_has_calib_tables(file_name: str) -> bool:
    """Check FITS file for existence of calibration tables."""
    has_gain_curves = False
    has_sys_temperatures = False

    with fits.open(file_name, mode="denywrite") as hdulist:
        for hdu in hdulist:
            if hdu.name == "SYSTEM_TEMPERATURE":
                has_sys_temperatures = True
            elif hdu.name == "GAIN_CURVE":
                has_gain_curves = True

    return has_gain_curves and has_sys_temperatures


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
