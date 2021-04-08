"""
Read and analyze PIMA fri-files.

@author: Petr Voytsik
"""

import math
from datetime import datetime
from typing import NamedTuple

import numpy as np
from scipy.stats import rv_continuous

ED = 2.0 * 6371e5  # Earth diameter
C = 2.99792458e10  # Speed of light


class pfd_dist(rv_continuous):
    """Theoretical distribution of SNR values in no-signal case."""

    def _cdf(self, x, n_eff, sigma_eff):
        exp_arg = -(np.pi * x ** 2) / (4 * sigma_eff ** 2)

        return np.power(1 - np.exp(exp_arg), n_eff)

    def _pdf(self, x, n_eff, sigma_eff):
        exp_arg = -(np.pi * x ** 2) / (4 * sigma_eff ** 2)

        x1 = np.log((np.pi / 2) * n_eff * x / sigma_eff ** 2)
        x2 = exp_arg
        x3 = (n_eff - 1) * np.log(1 - np.exp(exp_arg))

        return np.exp(x1 + x2 + x3)

    def _ppf(self, q, n_eff, sigma_eff):
        return 2 * sigma_eff * np.sqrt(-np.log(1 - (1 - q) ** (1 / n_eff)) / np.pi)


class PFDRec(NamedTuple):

    p2: float  # PFD=1e-2
    p3: float  # PFD=1e-3
    p4: float  # PFD=1e-4
    p5: float  # PFD=1e-5
    n_eff: float
    sigma_eff: float


SNR_LIMITS = {
    #
    # P-band
    #
    # Copy from C-band
    ("p", 64, 0.5, 285): PFDRec(5.24, 5.62, 5.97, 6.30, 5.523e04, 1.179),
    #
    # L-band
    #
    # ("l", 64, 1.0, 285): [5.4, 5.9, 6.3, 6.8],
    # 18724 points
    ("l", 64, 1, 285): PFDRec(5.32, 5.70, 6.05, 6.39, 5.877e04, 1.194),
    # ("l", 64, 1.0, 570): [5.5, 5.9, 6.2, 6.6],
    # 4860 points
    ("l", 64, 1, 570): PFDRec(5.49, 5.83, 6.15, 6.45, 8.305e05, 1.140),
    # ("l", 64, 1.0, 870): [5.6, 6.0, 6.3, 6.7],
    # 2678 points
    ("l", 64, 1, 870): PFDRec(5.65, 5.98, 6.30, 6.59, 2.010e06, 1.146),
    # ("l", 64, 1.0, 1170): [5.7, 6.1, 6.4, 6.7],
    # 1203 points
    ("l", 64, 1, 1170): PFDRec(5.74, 6.06, 6.36, 6.65, 6.254e06, 1.130),
    # ("l", 2048, 0.25, 285): [5.7, 6.0, 6.4, 6.7],
    # Copy from C-band
    ("l", 64, 0.5, 870): PFDRec(5.57, 5.90, 6.21, 6.50, 1.742e06, 1.133),
    # 12384 points
    ("l", 2048, 0.25, 285): PFDRec(5.85, 6.21, 6.55, 6.87, 1.004e06, 1.209),
    # ("l", 2048, 0.25, 570): [5.9, 6.2, 6.6, 6.9],
    # 12216 points
    ("l", 2048, 0.25, 570): PFDRec(6.09, 6.42, 6.75, 7.05, 5.887e06, 1.200),
    # ("l", 2048, 0.25, 870): [6.1, 6.5, 6.9, 7.2],
    # 12250 points
    ("l", 2048, 0.25, 870): PFDRec(6.21, 6.54, 6.85, 7.15, 1.883e07, 1.191),
    # ("l", 2048, 0.25, 1170): [6.1, 6.5, 6.8, 7.1],
    # 3112 points
    ("l", 2048, 0.25, 1170): PFDRec(6.28, 6.59, 6.89, 7.18, 6.125e07, 1.172),
    # 161 points
    ("l", 2048, 0.125, 570): PFDRec(6.18, 6.52, 6.85, 7.17, 4.567e06, 1.226),
    # 64 points
    ("l", 2048, 0.125, 870): PFDRec(6.37, 6.72, 7.05, 7.37, 6.604e06, 1.252),
    # 30 points
    ("l", 2048, 0.125, 1170): PFDRec(6.32, 6.61, 6.88, 7.15, 4.911e08, 1.129),
    #
    # C-band
    #
    # ("c", 64, 0.5, 285): [5.2, 5.6, 6.0, 6.3],
    # 31361 points
    ("c", 64, 0.5, 285): PFDRec(5.24, 5.62, 5.97, 6.30, 5.523e04, 1.179),
    # ("c", 64, 0.5, 570): [5.4, 5.8, 6.1, 6.5],
    # 7766 points
    ("c", 64, 0.5, 570): PFDRec(5.39, 5.72, 6.03, 6.32, 1.230e06, 1.108),
    # ("c", 64, 0.5, 870): [5.6, 5.9, 6.2, 6.6],
    # 4488 points
    ("c", 64, 0.5, 870): PFDRec(5.57, 5.90, 6.21, 6.50, 1.742e06, 1.133),
    # Just copy
    ("c", 64, 0.25, 570): PFDRec(5.57, 5.90, 6.21, 6.50, 1.742e06, 1.133),
    # ("c", 64, 0.5, 1170): [5.7, 6.0, 6.3, 6.6],
    # 1626 points
    ("c", 64, 0.5, 1170): PFDRec(5.70, 6.04, 6.35, 6.65, 2.156e06, 1.154),
    # ("c", 256, 0.5, 285.0): [5.4, 5.7, 6.0, 6.3],  # 478 points
    # 1299 points
    ("c", 256, 0.5, 285): PFDRec(5.36, 5.68, 6.00, 6.29, 8.022e05, 1.113),
    # ("c", 256, 0.5, 570): [5.6, 5.9, 6.2, 6.5],  # 229 points
    # 548 points
    ("c", 256, 0.5, 570): PFDRec(5.53, 5.83, 6.11, 6.38, 1.223e07, 1.071),
    # ("c", 2048, 0.125, 285): [5.8, 6.1, 6.5, 6.8],
    # 18883 points
    ("c", 2048, 0.125, 285): PFDRec(5.81, 6.16, 6.49, 6.80, 1.154e06, 1.194),
    # ("c", 2048, 0.125, 570): [6.1, 6.4, 6.7, 7.0],
    # 19157 points
    ("c", 2048, 0.125, 570): PFDRec(6.01, 6.34, 6.65, 6.95, 8.163e06, 1.177),
    # ("c", 2048, 0.125, 870): [6.2, 6.5, 6.8, 7.1],
    # 16650 points
    ("c", 2048, 0.125, 870): PFDRec(6.13, 6.45, 6.75, 7.04, 2.879e07, 1.165),
    # ("c", 2048, 0.125, 1170): [6.3, 6.6, 6.9, 7.2],
    # 4796 points
    ("c", 2048, 0.125, 1170): PFDRec(6.24, 6.55, 6.85, 7.14, 4.982e07, 1.170),
    # 12 points
    ("c", 2048, 0.5, 1170): PFDRec(6.01, 6.28, 6.54, 6.78, 1.097e09, 1.057),
    # Copy from K-band
    ("c", 64, 0.125, 570): PFDRec(5.70, 6.02, 6.34, 6.63, 2.884e06, 1.144),
    ("c", 64, 0.125, 870): PFDRec(5.89, 6.23, 6.54, 6.84, 4.311e06, 1.171),
    ("c", 2048, 0.5, 570): PFDRec(5.80, 6.11, 6.41, 6.70, 1.025e07, 1.128),
    #
    # K-band
    #
    # ("k", 64, 0.125, 285): [5.5, 5.9, 6.2, 6.5],
    # 4465 points
    ("k", 64, 0.125, 285): PFDRec(5.45, 5.79, 6.11, 6.41, 7.125e05, 1.136),
    # ("k", 64, 0.125, 570): [5.8, 6.1, 6.4, 6.7],
    # 1276 points
    ("k", 64, 0.125, 570): PFDRec(5.70, 6.02, 6.34, 6.63, 2.884e06, 1.144),
    # ("k", 64, 0.125, 870): [5.9, 6.3, 6.6, 6.9],
    # 542 points
    ("k", 64, 0.125, 870): PFDRec(5.89, 6.23, 6.54, 6.84, 4.311e06, 1.171),
    # ("k", 64, 0.125, 1170): [6.0, 6.3, 6.7, 7.0],
    # 175 points
    ("k", 64, 0.125, 1170): PFDRec(5.90, 6.20, 6.48, 6.76, 4.536e07, 1.109),
    # ("k", 256, 0.5, 285): [5.5, 5.8, 6.1, 6.5],  # 379 points
    # 18 points
    ("k", 64, 0.5, 570): PFDRec(5.42, 5.71, 5.99, 6.25, 1.135e07, 1.052),
    # 26 points
    ('k', 64, 0.5, 870): PFDRec(5.96, 6.38, 6.77, 7.15, 7.621e+04, 1.328),
    ('k', 64, 0.5, 1170): PFDRec(5.96, 6.38, 6.77, 7.15, 7.621e+04, 1.328),
    # 973 points
    ("k", 256, 0.5, 285): PFDRec(5.42, 5.74, 6.04, 6.33, 2.007e06, 1.099),
    # ("k", 256, 0.5, 570): [5.6, 5.8, 6.2, 6.5],  # 190 points
    # 538 points
    ("k", 256, 0.5, 570): PFDRec(5.68, 5.99, 6.29, 6.58, 5.974e06, 1.120),
    # ("k", 2048, 0.125, 285): [6.0, 6.3, 6.6, 6.9],  # 17893 points
    # 9833 points
    ("k", 2048, 0.125, 285): PFDRec(5.86, 6.19, 6.50, 6.80, 5.014e06, 1.161),
    # ("k", 2048, 0.125, 570): [6.2, 6.5, 6.8, 7.1],  # 10867 points
    # 7217 points
    ("k", 2048, 0.125, 570): PFDRec(6.09, 6.40, 6.70, 6.99, 3.037e07, 1.155),
    # ("k", 2048, 0.125, 870): [6.3, 6.6, 6.9, 7.2],  # 5131 points
    # 5018 points
    ("k", 2048, 0.125, 870): PFDRec(6.22, 6.54, 6.83, 7.12, 5.739e07, 1.164),
    # ("k", 2048, 0.125, 1170): [6.4, 6.7, 7.0, 7.3],  # 595 points
    # 863 points
    ("k", 2048, 0.125, 1170): PFDRec(6.32, 6.63, 6.92, 7.21, 1.004e08, 1.167),
    # ("k", 2048, 0.015625, 570): [6.5, 6.8, 7.1, 7.5],  # 731 points
    # 14 points
    ("k", 2048, 0.015625, 870): PFDRec(6.48, 6.83, 7.15, 7.46, 1.797e07, 1.245),
    # #
    # ("p", 64): [5.5, 5.9, 6.2, 6.6],
    # ("p", 2048): [6.0, 6.3, 6.8, 7.0],
    # ("l", 64): [5.5, 5.9, 6.2, 6.6],
    # ("l", 2048): [6.0, 6.3, 6.8, 7.0],
    # ("c", 64): [5.5, 5.8, 6.1, 6.5],
    # ("c", 2048): [6.1, 6.4, 6.7, 7.0],
    # ("c", 256): [5.6, 5.9, 6.2, 6.5],
    # ("k", 64): [5.7, 6.1, 6.4, 6.7],
    # ("k", 2048): [6.3, 6.6, 6.9, 7.2],
    # ("k", 256): [5.6, 5.9, 6.2, 6.5],
}


# PIMA fri-file obs fields
#  0 -> Obs #
#  1 -> Sca #
#  2 -> Time code %j-%H%M
#  3 -> Source name
#  4 -> Station 1
#  5 -> Station 2
#  7 -> SNR
# 11 -> Scan start time
# 12 -> Scan stop time
# 15 -> Ampl_lsq
# 23 -> Gr_del_lsq value
# 31 -> Ph_rat_lsq value
# 37 -> Ph_acc value
# 43 -> Phs_lsq value
# 79 -> Duration of scan
# 82 -> AP_len
# 84 -> U
# 85 -> V
# 89 -> Elevation of sta1
# 90 -> Elevation of sta2
# 94 -> Reference frequency
class Fri:
    """
    This class represents PIMA fringe file.

    """

    supported_versions = (
        "# PIMA Fringe results  v  1.00  Format version of 2010.04.05",
        "# PIMA Fringe results  v  1.01  Format version of 2014.02.08",
        "# PIMA Fringe results  v  1.02  Format version of 2014.12.24",
        "# PIMA Fringe results  v  1.2   Format version of 2019.04.20",
    )

    def __init__(self, file_name=None):
        """
        Parameters
        ----------
        file_name : str
            Name of the input file. If `None` create empty object.

        """
        self.records = []
        self.aux = {}  # Auxiliary user defined parameters

        if file_name:
            self.parse_file(file_name)

    def parse_file(self, file_name):
        """
        Parse PIMA fri-file.

        """
        started = False
        header = {}
        self.records.clear()

        with open(file_name) as fil:
            line = fil.readline()
            if not line.startswith(self.supported_versions):
                raise ValueError(f"{file_name} is not PIMA fri-file")

            for line in fil:
                if line.startswith("# PIMA_FRINGE started"):
                    started = True
                    header.clear()
                    continue
                if line.startswith("# PIMA_FRINGE ended"):
                    started = False
                    continue
                if not started:
                    continue
                toks = line.split()
                if len(toks) < 3:
                    continue
                if toks[0] == "#":
                    if toks[1] == "Control":
                        header["cnt_file"] = toks[-1]
                    elif toks[1] == "Experiment":
                        header["exper_code"] = toks[-1]
                    elif toks[1] == "Session":
                        header["session_code"] = toks[-1]
                    elif toks[1] == "FRINGE_FITTING_STYLE:":
                        header["FRINGE_FITTING_STYLE"] = toks[-1]
                    elif toks[1] == "FRIB.POLAR:":
                        header["polar"] = toks[2]
                    elif toks[1] == "FRIB.FINE_SEARCH:":
                        header["FRIB.FINE_SEARCH"] = toks[-1]
                    elif toks[1] == "PHASE_ACCELERATION:":
                        header["accel"] = float(toks[2].replace("D", "e"))
                    elif toks[1] == "FRIB.BEG_IFRQ:":
                        header["beg_ifrq"] = int(toks[2])
                    elif toks[1] == "FRIB.END_IFRQ:":
                        header["end_ifrq"] = int(toks[2])
                elif toks[6] != "FAILURE":
                    self.records.append({"status": "u", "pfd": 0.0})
                    self.records[-1].update(header)
                    self.records[-1]["obs"] = int(toks[0])
                    self.records[-1]["scan"] = int(toks[1])
                    self.records[-1]["time_code"] = toks[2]
                    self.records[-1]["source"] = toks[3]
                    self.records[-1]["sta1"] = toks[4]
                    self.records[-1]["sta2"] = toks[5]
                    try:
                        self.records[-1]["SNR"] = float(toks[7])
                    except ValueError:
                        self.records[-1]["SNR"] = 0
                    self.records[-1]["start_time"] = datetime.strptime(
                        toks[11], "%Y.%m.%d-%H:%M:%S.%f,"
                    )
                    self.records[-1]["stop_time"] = datetime.strptime(
                        toks[12], "%Y.%m.%d-%H:%M:%S.%f"
                    )
                    self.records[-1]["ampl_lsq"] = float(toks[15].replace("D", "e"))
                    self.records[-1]["delay"] = float(toks[23].replace("D", "e"))
                    self.records[-1]["rate"] = float(toks[31].replace("D", "e"))
                    self.records[-1]["ph_acc"] = float(toks[37].replace("D", "e"))
                    self.records[-1]["phs_lsq"] = float(toks[43].replace("D", "e"))
                    self.records[-1]["duration"] = float(toks[79].replace("D", "e"))
                    self.records[-1]["ap_len"] = float(toks[82].replace("D", "e"))
                    self.records[-1]["U"] = float(toks[84].replace("D", "e"))
                    self.records[-1]["V"] = float(toks[85].replace("D", "e"))
                    self.records[-1]["ref_freq"] = float(toks[94].replace("D", "e"))
                    self.records[-1]["elevation"] = [float(toks[89]), float(toks[90])]
                    # Calculated parameters
                    # UV-radius in lambda
                    self.records[-1]["uv_rad"] = math.hypot(
                        self.records[-1]["U"], self.records[-1]["V"]
                    )
                    # Position angle
                    self.records[-1]["PA"] = math.degrees(
                        math.atan2(self.records[-1]["U"], self.records[-1]["V"])
                    )
                    wave_len = C / self.records[-1]["ref_freq"]
                    # UV-radius in Earth diameters
                    self.records[-1]["uv_rad_ed"] = (
                        self.records[-1]["uv_rad"] * wave_len / ED
                    )

    def update_status(self, ch_num, snr_det_limits=None):
        """
        Update observations statuses acording to PFD.

        Parameters
        ----------
        ch_num : int
            Number of spectral channels.
        snr_det_limits : PFDRec, optional
            List of SNR values corresponding to PDF 1e-2, 1e-3, 1e-4, 1e-5.

        Notes
        -----
            This function is useful only for RadioAstron AGN survey.

        """
        if not self.records:
            return

        accum_length = self.records[0]["ap_len"]

        # Check if Session code in RA AGN survey format
        if "_" not in self.records[0]["session_code"]:
            return

        band = self.records[0]["session_code"].split("_")[1]

        max_scan_len = self.max_scan_length()
        if max_scan_len <= 300:
            scan_length = 285
        elif max_scan_len <= 600:
            scan_length = 570
        elif max_scan_len <= 900:
            scan_length = 870
        else:
            scan_length = 1170

        if not snr_det_limits:
            snr_det_limits = SNR_LIMITS[band, ch_num, accum_length, scan_length]

        my_dist = pfd_dist(name="my_dist")

        for rec in self.records:
            if ch_num <= 256 and abs(rec["delay"]) < 5e-7 and abs(rec["rate"]) < 1e-12:
                if rec["SNR"] >= snr_det_limits.p3:  # PFD = 1e-3
                    rec["status"] = "y"
            else:
                if rec["SNR"] < snr_det_limits.p2:  # PFD = 1e-2
                    rec["status"] = "n"
                elif rec["SNR"] >= snr_det_limits.p4:  # PFD = 1e-4
                    rec["status"] = "y"

                # Additional check for fringe rate
                if abs(rec["rate"]) > 5e-10 and rec["status"] == "y":
                    rec["status"] = "u"

            rec["pfd"] = 1.0 - my_dist.cdf(
                rec["SNR"], snr_det_limits.n_eff, snr_det_limits.sigma_eff
            )

    def rec_by_obs(self, obs):
        """
        Return record with given observation number.

        Parameters
        ----------
        obs : int
            Observation number.

        """
        result = None

        for rec in self.records:
            if rec["obs"] == obs:
                result = rec
                break

        return result

    def max_snr(self, station=None):
        """
        Return observation record with maximum SNR.

        Parameters
        ----------
        station : str, optional
            If station name is provided select observations with this station
            only.

        Returns
        -------
        result : dict
            Returns fri-record -- dictionary with parameters of the selected
            observation.

        """
        # Select observations with 'station'
        if station:
            records = [
                rec for rec in self.records if station in (rec["sta1"], rec["sta2"])
            ]
        else:
            records = self.records

        # Sort records by SNR
        records = sorted(records, key=lambda rec: rec["SNR"], reverse=True)
        if records:
            result = records[0]
        else:
            result = None

        return result

    def average_scan_length(self):
        """
        Return average length of observations.

        """
        if self.records:
            lengths = [rec["duration"] for rec in self.records]
            aver_len = math.fsum(lengths) / len(lengths)
        else:
            aver_len = 0

        return aver_len

    def max_scan_length(self):
        """
        Return maximum length of observations.

        """
        if self.records:
            lengths = [rec["duration"] for rec in self.records]
            max_len = max(lengths)
        else:
            max_len = 0

        return max_len

    def min_detected_snr(self):
        """
        Return minimum SNR of obsevations with detection status = 'y'.

        """
        if self.records:
            result = min([rec["SNR"] for rec in self.records if rec["status"] == "y"])
        else:
            result = 0

        return result

    def any_detections(self, station=None):
        """
        Return ``True`` if there is at least one observation with status 'y'.

        Parameters
        ----------
        station : str, optional
            If station name is provided select observations with this station
            only.

        """
        # Select observations with `station`
        if station:
            records = [
                rec for rec in self.records if station in (rec["sta1"], rec["sta2"])
            ]
        else:
            records = self.records

        return "y" in [rec["status"] for rec in records]

    def non_detections(self):
        """
        Return list of observation indices with status != 'y'.

        """
        return [rec["obs"] for rec in self.records if rec["status"] != "y"]

    def append(self, rec):
        """
        Append fri-file record to the end of the list.

        """
        self.records.append(rec)

    def remove_obs(self, obs_list):
        """
        Remove observations with indices in `obs_list`.

        Parameters
        ----------
        obs_list : list
            List of observation indices to remove.

        """
        if obs_list:
            self.records = [rec for rec in self.records if rec["obs"] not in obs_list]

    def __str__(self):
        out = "#Obs Timecode   Source      Sta1/Sta2         SNR    Delay    \
Rate      Accel      Base   Base  Length   PFD\n"
        for rec in self.records:
            accel = rec["accel"]
            if rec["FRIB.FINE_SEARCH"] == "ACC":
                accel = rec["ph_acc"]

            line = "{:>3}{:>10} {:>8} {:>8}/{:>8} {:8.2f} {:8.3f} \
{:10.3e} {:9.2e} {:8.2f} {:5.1f} {:6.1f} {:.1e} {:>1}\n".format(
                rec["obs"],
                rec["time_code"],
                rec["source"],
                rec["sta1"],
                rec["sta2"],
                rec["SNR"],
                1e6 * rec["delay"],
                rec["rate"],
                accel,
                1e-6 * rec["uv_rad"],
                rec["uv_rad_ed"],
                rec["duration"],
                rec["pfd"],
                rec["status"],
            )
            out = out + line

        return out

    def __getitem__(self, ind):
        return self.records[ind]

    def __len__(self):
        return len(self.records)
