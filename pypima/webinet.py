"""Functions to work with Webinet data archive."""

import os
from datetime import datetime

from .db import DataBase


def get_orbit_path(
    exper: str, db: DataBase, webinet_dir: str | os.PathLike | None = None
) -> str | None:
    """Return the path to the orbit file for a given experiment.

    Parameters
    ----------
        exper : str
            The experiment name.
        db : DataBase
            An instance of the DataBase class.
        webinet_dir : str, optional
            The base directory of the local copy of Webinet data archive.
            If None, download the file via FTP.

    Returns
    -------
        orbit path : str
            The full path to the orbit file or ``None`` if orbit is not found.
    """
    # Relative path to the reconstructed orbit files
    orbit_dir = "radioastron/oddata/reconstr"

    orbit_name = db.get_orbit_file(exper)
    if orbit_name is None:
        return None

    if webinet_dir is not None:
        if not os.path.isdir(webinet_dir):
            msg = f"Directory {webinet_dir} does not exist."
            raise ValueError(msg)

        orbit_path = os.path.join(webinet_dir, orbit_dir, orbit_name)
    else:
        orbit_path = f"ftp://webinet.asc.rssi.ru/{orbit_dir}/{orbit_name}"

    return orbit_path


def get_antab_path(
    exper: str,
    band: str,
    exper_date: datetime,
    webinet_dir: str | os.PathLike | None = None,
) -> str:
    """Return the path to the ANTAB file for a given experiment.

    Parameters
    ----------
        exper : str
            The experiment name.
        band : str
            Frequency band.
        exper_date : datetime
            The date of the experiment.
        webinet_dir : str, optional
            The base directory of the local copy of Webinet data archive.
            If None, download the file via FTP.

    Returns
    -------
        antab path : str
            The full path or URL to the ANTAB file.
    """
    date_str1 = exper_date.strftime("%Y_%m")
    date_str2 = exper_date.strftime("%Y_%m_%d")
    antab_path = "{0}/{1}/{2}_{3}/{3}{4}.antab2".format(
        "radioastron/ampcal", date_str1, date_str2, exper, band
    )

    if webinet_dir is not None:
        if not os.path.isdir(webinet_dir):
            msg = f"Directory {webinet_dir} does not exist."
            raise ValueError(msg)

        antab_path = os.path.join(webinet_dir, antab_path)
    else:
        antab_path = f"ftp://webinet.asc.rssi.ru/{antab_path}"

    return antab_path
