"""Module for handling UVFITS files and UV data.

Created on Thu Feb 11 14:00:15 2016

@author: Petr Voytsik
"""

import os
from dataclasses import InitVar, dataclass, field

import numpy as np
import numpy.typing as npt
from astropy.io import fits
from astropy.time import Time


class UVFitsError(Exception):
    """Base class for exceptions in this module."""

    def __init__(self, message, file_name):
        self.file_name = file_name
        self.message = message

    def __str__(self):
        return f"{self.message}: {self.file_name}"


def get_float(header: fits.Header, key: str) -> float:
    """Get float value from FITS header."""
    val = header[key]  # Will raise KeyError if not found

    if not isinstance(val, float):
        raise TypeError(f"Invalid {key} value: {val}")

    return val


def get_int(header: fits.Header, key: str) -> int:
    """Get int value from FITS header."""
    val = header[key]

    if not isinstance(val, int):
        raise TypeError(f"Invalid {key} value: {val}")

    return val


def baseline_decode(
    baseline: npt.NDArray,
) -> tuple[npt.NDArray, npt.NDArray, npt.NDArray]:
    """Decode `baseline` index and return subarray and antennas numbers."""
    ant1 = np.asarray(baseline, dtype=int) // 256
    ant2 = np.mod(np.asarray(baseline, dtype=int), 256)
    subarr = np.asarray(np.round(100 * (baseline - np.round(baseline))), dtype=int)

    return subarr, ant1, ant2


NUM2STOKES = {
    1: "I",
    2: "Q",
    3: "U",
    4: "V",
    -1: "RR",
    -2: "LL",
    -3: "RL",
    -4: "LR",
    -5: "XX",
    -6: "YY",
    -7: "XY",
    -8: "YX",
}


@dataclass
class UVFits:
    """Class for reading and handling UVFITS data files."""

    hdulist: fits.HDUList = field(init=False, compare=False, repr=False)
    antenna_table: list[dict] = field(init=False, default_factory=list)
    freq_table: list[dict] = field(init=False, default_factory=list)

    # Header
    freq: float = field(init=False, default=0.0)
    no_if: int = field(init=False, default=0)
    gcount: int = field(init=False, default=0)
    source: str = field(init=False, default="")
    exper_name: str = field(init=False, default="")
    band: str = field(init=False, default="")
    stokes: str = field(init=False, default="")

    # UV-data
    u_raw: npt.NDArray = field(init=False)
    v_raw: npt.NDArray = field(init=False)
    baselines: npt.NDArray = field(init=False)
    ant1_inds: npt.NDArray = field(init=False)
    ant2_inds: npt.NDArray = field(init=False)
    subarrays: npt.NDArray = field(init=False)
    # self.dates1 = None
    # self.dates2 = None
    dates: npt.NDArray = field(init=False)
    amplitudes: npt.NDArray = field(init=False)
    phases: npt.NDArray = field(init=False)
    weights: npt.NDArray = field(init=False)
    inttime: npt.NDArray = field(init=False)

    # Input
    path: InitVar[str | os.PathLike]

    def __post_init__(self, path: str | os.PathLike) -> None:
        self.hdulist = fits.open(path)

        self._read_header()
        self._read_an_table()
        self._read_fq_table()
        self._read_uv_data()

    # Support the 'with' statement
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def _read_header(self) -> None:
        """Read UVFITS file header."""
        primary_hdu = self.hdulist[0]

        if not isinstance(primary_hdu, fits.GroupsHDU):
            raise UVFitsError(
                "Invalid UVFITS file: Primary HDU is not GroupsHDU", self.file_name
            )

        header = primary_hdu.header

        if "GROUPS" not in header or not header["GROUPS"]:
            raise UVFitsError(
                "Invalid UVFITS file: no GROUPS in header", self.file_name
            )

        if header["NAXIS1"] != 0:
            raise UVFitsError("Invalid UVFITS file: NAXIS1 is not zero", self.file_name)

        toks = str(header["OBSERVER"]).split("_")
        if len(toks) == 2:  # Exper code in RA AGN survey format
            self.exper_name = toks[0]
            self.band = toks[1]
        else:  # General case
            self.exper_name = str(header["OBSERVER"])

        self.source = str(header["OBJECT"])
        self.gcount = get_int(header, "GCOUNT")
        if self.gcount <= 0:
            raise UVFitsError(f"Invalid GCOUNT value: {self.gcount}", self.file_name)

        # Go through the array axis
        for ind in range(2, 1 + get_int(header, "NAXIS")):
            naxis = get_int(header, f"NAXIS{ind}")
            ctype = header[f"CTYPE{ind}"]
            crval = get_float(header, f"CRVAL{ind}")
            crpix = get_float(header, f"CRPIX{ind}")
            cdelt = get_float(header, f"CDELT{ind}")

            if ctype == "FREQ":
                if naxis != 1:
                    raise UVFitsError(
                        "Multiple frequency setup not supported", self.file_name
                    )
                self.freq = crval
            elif ctype == "RA" or ctype == "DEC":
                if naxis != 1:
                    raise UVFitsError("Multi source file not supported", self.file_name)
            elif ctype == "STOKES":
                if naxis != 1:
                    raise UVFitsError("Multiple Stokes not supported", self.file_name)
                val = int(crval + (1 - crpix) * cdelt)
                self.stokes = NUM2STOKES[val]

    def _read_an_table(self) -> None:
        """Read antenna table(s)."""
        # Look for all AIPS AN tables
        for hdu in self.hdulist:
            if isinstance(hdu, fits.BinTableHDU):
                if "EXTNAME" in hdu.header and hdu.header["EXTNAME"] == "AIPS AN":
                    anname = hdu.data.field("ANNAME")
                    nosta = hdu.data.field("NOSTA").tolist()
                    self.antenna_table.append(dict(zip(nosta, anname)))

        if not self.antenna_table:
            raise UVFitsError("Could not find AN table", self.file_name)

    def _read_fq_table(self) -> None:
        """Read frequency table."""
        fq_table = self.hdulist["AIPS FQ"]

        if not isinstance(fq_table, fits.BinTableHDU):
            raise UVFitsError("Invalid FQ table: not BinTableHDU", self.file_name)

        self.no_if = get_int(fq_table.header, "NO_IF")
        if self.no_if <= 0:
            raise UVFitsError(
                f"Invalid NO_IF value in FQ table: {self.no_if}", self.file_name
            )

        # Assume we have only one frequency setup
        freq_data = fq_table.data[0]

        if self.no_if == 1:
            self.freq_table.append(
                {
                    "if_freq": freq_data["IF FREQ"],
                    "ch_width": freq_data["CH WIDTH"],
                    "sideband": freq_data["SIDEBAND"],
                }
            )
        elif self.no_if >= 2:
            for if_num in range(self.no_if):
                self.freq_table.append(
                    {
                        "if_freq": freq_data["IF FREQ"][if_num],
                        "ch_width": freq_data["CH WIDTH"][if_num],
                        "sideband": freq_data["SIDEBAND"][if_num],
                    }
                )

    def _read_uv_data(self):
        """Read UV data."""
        primary_hdu = self.hdulist[0]
        assert isinstance(primary_hdu, fits.GroupsHDU), "Primary HDU is not GroupsHDU"
        data = primary_hdu.data
        assert isinstance(data, fits.GroupData), "Primary HDU data is not GroupsData"

        # Assume we have only files produced by PIMA
        self.u_raw = data.par("UU")
        self.v_raw = data.par("VV")
        # self.dates1 = data.par(4)
        # self.dates2 = data.par(5)
        self.inttime = data.par("INTTIM")
        self.dates = np.asarray(
            Time(data.par(4), data.par(5), scale="utc", format="jd").datetime
        )
        self.baselines = data.par("BASELINE")
        self.subarrays, self.ant1_inds, self.ant2_inds = baseline_decode(self.baselines)

        vis = data.data[:, 0, 0, :, 0, 0, 0] + data.data[:, 0, 0, :, 0, 0, 1] * 1j
        self.amplitudes = np.absolute(vis)
        self.phases = np.angle(vis)
        self.weights = data.data[:, 0, 0, :, 0, 0, 2]

    @property
    def file_name(self) -> str:
        """Return file name."""
        filename = self.hdulist.filename()
        if not isinstance(filename, str):
            raise RuntimeError("Invalid filename type")

        return filename

    def get_ant_by_ind(self, ind: int) -> tuple[str, str]:
        """Return antenna names for given baseline."""
        subarr = self.subarrays[ind]
        ant1 = self.ant1_inds[ind]
        ant2 = self.ant2_inds[ind]
        ant1_name = self.antenna_table[subarr][ant1]
        ant2_name = self.antenna_table[subarr][ant2]

        return ant1_name, ant2_name

    def close(self):
        """Close HDU list."""
        self.hdulist.close()
