"""Some useful tools for PIMA data analysis."""

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
