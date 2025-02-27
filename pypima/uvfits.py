"""Module for handling UVFITS files and UV data.

Created on Thu Feb 11 14:00:15 2016

@author: Petr Voytsik
"""

import numpy as np
from astropy.io import fits
from astropy.time import Time


class UVFitsError(Exception):
    """Base class for exceptions in this module."""

    def __init__(self, message, file_name):
        self.file_name = file_name
        self.message = message

    def __str__(self):
        return f"{self.message}: {self.file_name}"


def baseline_decode(baseline):
    """Decode `baseline` index and return subarray and antennas numbers."""
    ant1 = np.asarray(baseline, dtype=int) // 256
    ant2 = np.mod(np.asarray(baseline, dtype=int), 256)
    subarr = np.asarray(np.round(100 * (baseline - np.round(baseline))), dtype=int)

    return subarr, ant1, ant2


class UVFits:
    """Class for reading and handling UVFITS data files."""

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

    def __init__(self, file_name: str):
        # Init class variables
        self.file_name = file_name
        self.antenna_table = []
        self.freq_table = []

        # Header
        self.freq = 0.0
        self.no_if = 0
        self.gcount = 0
        self.source = None
        self.exper_name = None
        self.band = None

        # UV-data
        self.u_raw = None
        self.v_raw = None
        self.baselines = None
        self.ant1_inds = None
        self.ant2_inds = None
        self.subarrays = None
        self.dates1 = None
        self.dates2 = None
        self.dates = None
        self.amplitudes = None
        self.phases = None
        self.weights = None
        self.inttime = None

        self.hdulist = fits.open(file_name)

        self._read_header()
        self._read_an_table()
        self._read_fq_table()
        self._read_uv_data()

    # Support the 'with' statement
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def _read_header(self):
        """Read FITS file head."""
        header = self.hdulist[0].header

        if not header["groups"] or header["NAXIS1"] != 0:
            raise UVFitsError("", self.file_name)

        toks = header["OBSERVER"].split("_")
        if len(toks) == 2:  # Exper code in RA AGN survey format
            self.exper_name = toks[0]
            self.band = toks[1]
        else:  # General case
            self.exper_name = header["OBSERVER"]

        self.source = header["OBJECT"]
        self.gcount = header["GCOUNT"]
        # Go through the array axis
        for ind in range(2, 1 + header["NAXIS"]):
            naxis = header[f"NAXIS{ind}"]
            ctype = header[f"CTYPE{ind}"]
            crval = header[f"CRVAL{ind}"]
            crpix = header[f"CRPIX{ind}"]
            cdelt = header[f"CDELT{ind}"]

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
                val = crval + (1 - crpix) * cdelt
                self.stokes = self.NUM2STOKES[val]

    def _read_an_table(self):
        """Read antenna table."""
        subarray = 0
        for hdu in self.hdulist:
            if hdu.name == "AIPS AN":
                self.antenna_table.append({})
                for row in hdu.data:
                    anname = row["ANNAME"]
                    nosta = row["NOSTA"]
                    self.antenna_table[subarray][nosta] = anname
                subarray += 1

        if not self.antenna_table:
            raise UVFitsError("Could not find AN table", self.file_name)

    def _read_fq_table(self):
        """Read frequency table."""
        self.no_if = self.hdulist["AIPS FQ"].header["NO_IF"]

        # Assume we have only one frequency setup
        freq_data = self.hdulist["AIPS FQ"].data[0]

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
        else:
            raise UVFitsError(f"Invalid NO_IF value: {self.no_if}", self.file_name)

    def _read_uv_data(self):
        """Read UV data."""
        data = self.hdulist[0].data
        self.u_raw = data.par("UU")
        self.v_raw = data.par("VV")
        self.dates1 = data.par(4)
        self.dates2 = data.par(5)
        self.inttime = data.par("INTTIM")
        self.dates = Time(self.dates1, self.dates2, scale="utc", format="jd").datetime
        self.baselines = data.par("BASELINE")
        self.subarrays, self.ant1_inds, self.ant2_inds = baseline_decode(self.baselines)

        vis = data.data[:, 0, 0, :, 0, 0, 0] + data.data[:, 0, 0, :, 0, 0, 1] * 1j
        self.amplitudes = np.absolute(vis)
        self.phases = np.angle(vis)
        self.weights = data.data[:, 0, 0, :, 0, 0, 2]

    def get_ant_by_ind(self, ind):
        """Return antenna names for given baseline."""
        subarr = self.subarrays[ind]
        ant1 = self.ant1_inds[ind]
        ant2 = self.ant2_inds[ind]
        ant1_name = self.antenna_table[subarr][ant1]
        ant2_name = self.antenna_table[subarr][ant2]

        return ant1_name, ant2_name

    def print_info(self):
        """Print FITS file info to the standard output."""
        self.hdulist.info()

    def print_uv(self):
        """Print UV data to the standard output."""
        for ind in range(self.gcount):
            ant1_name, ant2_name = self.get_ant_by_ind(ind)

            for if_num in range(self.no_if):
                uu = self.u_raw[ind] * (self.freq + self.freq_table[if_num]["if_freq"])
                vv = self.v_raw[ind] * (self.freq + self.freq_table[if_num]["if_freq"])

                amp = self.amplitudes[ind, if_num]
                pha = self.phases[ind, if_num]
                weight = self.weights[ind, if_num]

                print(
                    ind + 1,
                    if_num + 1,
                    ant1_name,
                    ant2_name,
                    self.dates[ind],
                    f"{amp:10.6f} {pha:8.5f} {weight:10.1f}",
                    f"{uu:15.1f} {vv:15.1f}",
                )

    def useful_antennas(self):
        """Return set of antenna names."""
        antennas = set()

        for ind in range(self.gcount):
            ant1_name, ant2_name = self.get_ant_by_ind(ind)
            for if_num in range(self.no_if):
                if self.amplitudes[ind, if_num] > 0 and self.weights[ind, if_num] > 0:
                    antennas.add(ant1_name)
                    antennas.add(ant2_name)

        return antennas

    def close(self):
        """Close HDU list."""
        self.hdulist.close()
