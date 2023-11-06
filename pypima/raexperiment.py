"""
Created on Sun Dec 29 04:02:35 2013

@author: Petr Voytsik
"""

import logging
import os.path
import shutil
import subprocess
import threading
from datetime import datetime
from io import BytesIO

import numpy as np
import pandas as pd
import pycurl

import pypima.pima

from .fri import Fri
from .pima import ActaFile, Pima, bpas_log_snr_new
from .uvfits import UVFits


class Error(Exception):
    """Raised when RaExperiment error occurs."""

    def __init__(self, exper, band, msg):
        self.exper = exper
        self.band = band
        self.msg = msg
        # self.time = str(datetime.now())

    def __str__(self):
        return f"{self.exper}({self.band}): {self.msg}"


class RaExperiment:
    """Describe experiment in RadioAstron AGN survey."""

    def __init__(
        self,
        experiment_code,
        band,
        data_base,
        data_dir=None,
        uv_fits=None,
        orbit=None,
        gvlbi=False,
    ):
        """
        Parameters
        ----------

        experiment_code : str
            Experiment code.
        band : srt
            One letter frequency band code.
        data_base : pypima.db.DB
            pypima.db.DB instance.
        data_dir : str, optional
            Directory for FITS-IDI. If ``None`` working directory of the
            current experiment is used.
        uv_fits : str or list, optional
            Path(s) to the data file (FITS-IDI). If ``None`` (default) get a
            file name from data base and download file from the FTP archive.
        orbit : str, optional
            Path to a reconstructed orbit file. If ``None`` (default),
            download it from the FTP archive.
        gvlbi : bool, optional
            If ``True``, process ground only of part of the experiment (GVLBI
            FITS file).

        """
        # First, set common variables
        self.exper = experiment_code.lower()
        self.band = band.lower()
        self.db = data_base
        self.gvlbi = gvlbi
        self.sta_ref = "RADIO-AS"
        self.run_id = 0  # Record id in pima_runs database table
        self.lock = threading.Lock()  # Lock for FITS file downloading control
        self.logger = logging.getLogger(f"{self.exper}({self.band})")
        self.antab = None
        self.calibration_loaded = False
        self.split_time_aver = 0
        self.pima = None
        self.uv_fits = uv_fits
        self.orbit = orbit
        self.fri = None  # Result of last fringe fitting
        self.scan_part = 0
        self.bad_obs_set = set()  # Set of bad obs (autospec)
        self.acta_files = {}  # Store autospectra data

        # dict (POLAR, frq_grp) -> BANDPASS_FILE
        self.bpass_files = {
            ("RR", 1): "",
            ("LL", 1): "",
            ("RL", 1): "",
            ("LR", 1): "",
        }

        if self.band not in ("p", "l", "c", "k"):
            self._error(f"unknown band {band}")

        self.pima_dir = os.getenv("PIMA_DIR")
        self.exp_dir = os.getenv("pima_exp_dir")
        if not self.exp_dir:
            self._error("Environment variable $pima_exp_dir is not set")

        self.pima_scr = os.getenv("pima_scr_dir")

        # Work directory path
        self.work_dir = os.path.join(self.exp_dir, self.exper + "_auto")
        if self.gvlbi:
            self.work_dir += "_gvlbi"

        #  Select directory for raw data from a correlator
        if data_dir:
            self.data_dir = os.path.join(data_dir, self.exper)
        else:
            self.data_dir = self.work_dir

        if len(self.data_dir) >= pypima.pima.UVFILE_NAME_LEN - 1:
            self._error(
                "Length of the data_dir path must be less than {} \
bytes".format(
                    pypima.pima.UVFILE_NAME_LEN - 1
                )
            )

        # PIMA control file path
        self.cnt_file_name = os.path.join(
            self.work_dir, f"{self.exper}_{self.band}_pima.cnt"
        )

    def init_workdir(self):
        """Create working directory and PIMA control file."""
        # Create work directory
        os.makedirs(self.work_dir, exist_ok=True)

        os.chdir(self.work_dir)

        # Create PIMA control file
        self._mk_cnt()

        self.pima = Pima(self.exper, self.band, self.work_dir)

        if self.uv_fits:
            self.pima.update_cnt({"UV_FITS:": self.uv_fits})

        # Only one sideband at P-band
        if self.band == "p":
            self.pima.update_cnt({"END_FRQ:": "1"})

        # Do not restrict delay rate window
        self.pima.update_cnt({"FRIB.RATE_WINDOW_WIDTH:": "1.0D-8"})

        # Setup staging dir
        staging_dir = os.getenv("PYPIMA_STAGING_DIR")

        if staging_dir:
            staging_dir = os.path.join(staging_dir, self.exper)
            os.makedirs(staging_dir, exist_ok=True)
            self.pima.update_cnt({"STAGING_DIR:": staging_dir})
        else:
            self.pima.update_cnt({"STAGING_DIR:": "NO"})

        # Select the best fftw wisdom file
        for thread_num in (8, 4, 2, 1):
            for size in ("huge", "big", "small"):
                wis_file = f"pima_{size}_measure_{thread_num:d}thr.wis"
                wis_file = os.path.join(self.pima_dir, "share", "pima", wis_file)

                if os.path.isfile(wis_file):
                    self.pima.update_cnt(
                        {"FFT_CONFIG_FILE:": wis_file, "NUM_THREADS:": thread_num}
                    )
                    os.environ["OMP_NUM_THREADS"] = str(thread_num)

                    return

    def _error(self, msg):
        """Raise pima.Error exception."""
        self.logger.error(msg)
        raise Error(self.exper, self.band, msg)

    def _mk_cnt(self):
        """Make new cnt-file from template."""
        cnt_templ_name = os.path.join(
            self.pima_dir, "share", "pima", "TEMPLATE_pima.cnt"
        )

        sess_code = f"{self.exper}_{self.band}"
        fringe_file = os.path.join(self.work_dir, sess_code + ".fri")
        frires_file = os.path.join(self.work_dir, sess_code + ".frr")

        if self.band in ("l", "p"):
            polar = "RR"
        else:
            polar = "LL"

        with open(cnt_templ_name) as cnt_templ, open(
            self.cnt_file_name, "w"
        ) as cnt_file:
            for line in cnt_templ:
                if "@CDATE@" in line:
                    line = line.replace("@CDATE@", str(datetime.now()))
                elif "@pima_dir@" in line:
                    line = line.replace("@pima_dir@", self.pima_dir)
                elif line.startswith("SESS_CODE:"):
                    line = line.replace("@sess_code@", sess_code)
                elif line.startswith("BAND:"):
                    line = line.replace("@band@", self.band.upper())
                elif line.startswith("EXPER_DIR:"):
                    line = line.replace("@exper_dir@", self.pima_scr)
                elif line.startswith("FRINGE_FILE:"):
                    line = line.replace("@fringe_file@", fringe_file)
                elif line.startswith("FRIRES_FILE:"):
                    line = line.replace("@frires_file@", frires_file)
                elif line.startswith("STA_REF:") and self.sta_ref:
                    line = "{:<20}{}\n".format("STA_REF:", self.sta_ref)
                elif line.startswith("EPHEMERIDES_FILE:") and self.orbit:
                    line = line.replace("@ephemerides_file@", self.orbit)
                elif line.startswith("POLAR:") or line.startswith("SPLT.POLAR:"):
                    line = line.replace("@polar@", polar)
                cnt_file.write(line)

    def _download_fits(self, force_small=False):
        """
        Download FITS-file from the FTP archive.

        """
        fits_url, size = self.db.get_uvfits_url(
            self.exper, self.band, self.gvlbi, force_small
        )

        if not fits_url:
            self._error("Could not find FITS file name in DB")

        # Delete spaces in filename
        uv_fits = os.path.join(
            self.data_dir, os.path.basename(fits_url).replace(" ", "")
        )

        if len(uv_fits) > pypima.pima.UVFILE_NAME_LEN:
            uv_fits = uv_fits[: pypima.pima.UVFILE_NAME_LEN]
            assert uv_fits[-1] == "/"

        if os.path.isfile(uv_fits) and os.path.getsize(uv_fits) == size:
            self.logger.info("File %s already exists", uv_fits)
        else:
            os.makedirs(self.data_dir, exist_ok=True)
            self.logger.info("Start downloading file %s...", fits_url)
            try:
                with open(uv_fits, "wb") as fil:
                    _download_it(fits_url, fil, max_retries=2)
            except pycurl.error as err:
                self._error(f"Could not download file {fits_url}: {err}")

            self.logger.info("FITS-file downloading is complete")

        # We use self.uv_fits as a flag of FITS file existence, so set it at
        # the end of this function
        self.uv_fits = uv_fits

    def _get_orbit(self):
        """Download reconstructed orbit file from FTP server."""
        orbit_url = self.db.get_orbit_url(self.exper)

        if not orbit_url:
            self._error("Could not find reconstructed orbit")

        self.orbit = os.path.join(self.work_dir, os.path.basename(orbit_url))

        self.logger.info("Start downloading orbit file %s ...", orbit_url)

        buffer = BytesIO()
        try:
            _download_it(orbit_url, buffer)
        except pycurl.error as err:
            self._error(f"Could not download file {orbit_url}: {err}")

        orb_data = buffer.getvalue().decode().splitlines()

        with open(self.orbit, "w") as orb_file:
            if not orb_data[0].startswith("CCSDS_OEM_VERS"):
                orb_file.write("CCSDS_OEM_VERS = 2.0\n")

            # Fix meta information
            for line in orb_data:
                if line.startswith("CENTER_NAME"):
                    line = "CENTER_NAME   = Earth Barycenter"
                elif line.startswith("OBJECT_NAME"):
                    line = "OBJECT_NAME   = RADIO-ASTRON"
                elif line.startswith("CREATION"):
                    line = line.replace("CREATION DATE", "CREATION_DATE")
                elif line.startswith("STOP_TIME") and len(line) < 20:
                    for back_line in reversed(orb_data):
                        cols = back_line.split()
                        if cols:
                            break
                    line = line.strip() + " " + cols[0]

                orb_file.write(line + "\n")

        self.logger.info("Orbit downloading is complete")
        self.pima.update_cnt({"EPHEMERIDES_FILE:": self.orbit})

    def _get_antab(self):
        """Download ANTAB-file from the FTP server."""
        # antab_url = self.db.get_antab_url(self.exper, self.band)

        # Make ANTAB url
        date_str1 = self.pima.exper_info.nominal_start.strftime("%Y_%m")
        date_str2 = self.pima.exper_info.nominal_start.strftime("%Y_%m_%d")
        url_base = "ftp://webinet.asc.rssi.ru/radioastron/ampcal"
        antab_url = "{0}/{1}/{2}_{3}/{3}{4}.antab2".format(
            url_base, date_str1, date_str2, self.exper, self.band
        )

        antab_dir = os.path.join(self.work_dir, "antab")
        os.makedirs(antab_dir, exist_ok=True)

        antab_file = os.path.join(antab_dir, os.path.basename(antab_url) + ".orig")
        self.logger.info("Start downloading file %s", antab_url)

        try:
            with open(antab_file, "wb") as fil:
                _download_it(antab_url, fil)

            self.logger.info("ANTAB-file downloading is complete.")
            self.antab = self._fix_antab(antab_file)
        except pycurl.error as err:
            self.antab = None
            self.logger.warning("Could not download file %s: %s", antab_url, err)

    def _fix_antab(self, antab):
        """Fix ANTAB file."""
        if not antab or not os.path.isfile(antab):
            return

        new_antab = antab.replace(".orig", "")

        # ANTAB file already exists and prepared
        if antab == new_antab:
            return new_antab

        freq_setup = self.pima.frequencies
        freq_list = [1e-6 * freq.freq for freq in freq_setup]

        # Should we fix frequency setup?
        fix_freq = False
        if self.band != "p" and freq_setup[0].side_band != -1:
            self.logger.warning("enable sideband fix for ANTAB")
            fix_freq = True

        sta_list = self.pima.station_list(ivs_name=False)

        with open(antab) as inp, open(new_antab, "w") as out:
            magic = inp.readline()
            if not magic.startswith("! Produced by: TSM"):
                self.logger.warning(
                    "antab file %s does NOT have magic in the" "first line", antab
                )
                return

            # Do not forget to write a 'magic' line to the output file
            out.write(magic)

            for line in inp:
                line = line.strip()

                # Skip empty lines
                if not line:
                    continue

                if line.startswith("POLY") and line.endswith("/"):
                    line = line.replace("/", " /")
                elif line.startswith("/") and len(line) > 1:
                    line = line.replace("/", "/ ", 1)
                elif not line.startswith("!") and "!" in line:
                    line = line.split("!")[0].strip()

                # VLA (raes11a and friends)
                if "YY" in sta_list and "Y27" in line:
                    line = line.replace("Y27", "YY")
                elif "KZ" in sta_list and "KL" in line:
                    line = line.replace("KL", "KZ")
                elif "EF" in sta_list and "EB" in line:
                    line = line.replace("EB", "EF")
                elif "WB" in sta_list and "WB1" in line:
                    line = line.replace("WB1", "WB")

                toks = line.split()

                # # Fix EF C-band channels table
                # if len(toks) == 10 and toks[0] == "!" and toks[1].isdigit():
                #     self.logger.warning("fix EF C-band channels table")
                #     toks.insert(2, "6cm")

                if fix_freq and len(toks) > 9 and toks[1].isdigit():
                    if toks[6] == "L":
                        toks[6] = "U"
                        toks[9] = "{:.2f}MHz".format(freq_list[0])

                # Deselect stations
                if toks[0] == "TSYS" and len(toks) > 4:
                    toks[4] = toks[4].upper()
                    if toks[4] not in sta_list:
                        self.logger.warning("deselect %s from ANTAB", toks[4])
                        toks.insert(0, "!")
                elif toks[0] == "GAIN":
                    # EF, L-band GAINs
                    if toks[-1] == "/":
                        for ind in range(len(toks)):
                            if toks[ind].startswith("POLY"):
                                toks[ind] = "\n" + toks[ind]
                                break

                    # Comment out GAIN line for different frequency
                    for tok in toks:
                        if tok.startswith("FREQ=") and "," in tok:
                            fr1, fr2 = tok.replace("FREQ=", "").split(",")
                            fr1 = float(fr1)  # Lower limit
                            fr2 = float(fr2)  # Upper limit
                            if min(freq_list) < fr1 or max(freq_list) > fr2:
                                self.logger.warning(
                                    "deselect GAIN due to " "%s is out of freq range",
                                    tok,
                                )
                                toks.insert(0, "!")
                                break

                elif toks[0] == "/" and len(toks) > 1:
                    toks[1] = "\n" + toks[1]

                out.write(" ".join(toks) + "\n")

        return new_antab

    def load(
        self,
        download_only=False,
        update_db=False,
        scan_length=1200,
        scan_part=1,
        force_small=False,
        beg_frq=None,
        end_frq=None,
    ):
        """
        Download data, run pima load, and do some checks.

        Parameters
        ----------
        download_only : bool, optional
            If ``True``, download FITS-file and return.
        update_db : bool, optional
            If ``True``, update database with experiment information.
        scan_length : float, optional
            Set maximum length of scan. Default is 20 min.
        scan_part : int, optional
            1 is full scan, 2 is half of scan. In general `scan_part` can be
            used as run index.
        force_small : bool, optional
            Force laod 64-channels FITS_IDI files.
        beg_frq : int, optional
            Start IF index.
        end_frq : int, optional
            End IF index.

        """
        self.scan_part = scan_part
        self.bad_obs_set.clear()
        self.acta_files.clear()

        # If self.uv_fits is not None assume FITS file already exists
        with self.lock:
            if not self.uv_fits:
                self._download_fits(force_small)

        if download_only:
            return

        self.pima.update_cnt({"UV_FITS:": self.uv_fits})

        # Check free space in STAGING_DIR
        staging_dir = self.pima.cnt_params["STAGING_DIR:"]

        if os.path.isdir(staging_dir):
            df = shutil.disk_usage(staging_dir)

            if isinstance(self.uv_fits, list):
                fits_file = self.uv_fits[0]

                uv_fits_size = sum([os.path.getsize(file) for file in self.uv_fits])
            else:
                fits_file = self.uv_fits
                uv_fits_size = os.path.getsize(fits_file)

            fits_file_staging = os.path.join(staging_dir, os.path.basename(fits_file))

            # Free space < size of file + 10%
            if not os.path.isfile(fits_file_staging) and df.free < 1.1 * uv_fits_size:
                self.logger.warning("Not enough space in STAGING_DIR")
                self.pima.update_cnt({"STAGING_DIR:": "NO"})
                shutil.rmtree(staging_dir)

        if self.orbit is None:
            self._get_orbit()

        # Set maximum scan length
        self.logger.info("Set maximum scan length to %s s", scan_length)
        self.pima.update_cnt(
            {"MAX_SCAN_LEN:": str(scan_length), "SCAN_LEN_USED:": str(scan_length)}
        )

        if update_db:
            if isinstance(self.uv_fits, list):
                uv_fits = self.uv_fits[0]  # Only first FITS to DB
            else:
                uv_fits = self.uv_fits

            self.run_id = self.db.add_exper_info(
                self.exper, self.band, os.path.basename(uv_fits), scan_part
            )

        self.pima.load()

        if update_db:
            self.db.update_exper_info(self.pima.exper_info, self.run_id)
            if scan_part == 1:
                self.db.model2db(self.run_id, self.pima.clock_model())

        #
        # Various checks and setups
        #
        if self.pima.obs_number == 0:
            self._error("ZERO observations have been loaded")

        sou_dist = self.pima.source_dist()
        for source, distance in sou_dist.items():
            if distance > 1.0:
                self._error(f"Dist = {distance} arcsec for source {source}")

        # Set number of IFs
        if beg_frq:
            if beg_frq < 1 or beg_frq > self.pima.exper_info.if_num:
                self._error(
                    "beg_frq must be in range [1, {}]".format(
                        self.pima.exper_info.if_num
                    )
                )
            else:
                self.pima.update_cnt({"BEG_FRQ:": beg_frq})
        else:
            self.pima.update_cnt({"BEG_FRQ:": 1})

        if end_frq:
            if (
                end_frq < int(self.pima.cnt_params["BEG_FRQ:"])
                or end_frq > self.pima.exper_info.if_num
            ):
                self._error(
                    "end_frq must be in range [{}, {}]".format(
                        self.pima.cnt_params["BEG_FRQ:"], self.pima.exper_info.if_num
                    )
                )
            else:
                self.pima.update_cnt({"END_FRQ:": end_frq})
        else:
            self.pima.update_cnt({"END_FRQ:": self.pima.exper_info.if_num})

        if "RADIO-AS" not in self.pima.station_list():
            self.logger.warning("RADIO-AS is not in station list")
            self.sta_ref = self.pima.station_list()[0]
            self.pima.update_cnt({"STA_REF:": self.sta_ref})

        desel_nam = self.pima.number_of_deselected_points
        if desel_nam > 10:
            self.logger.warning("Total number of deselected points is %s", desel_nam)

        # Save memory by reducing oversampling
        if self.pima.ap_minmax[0] < 0.1:
            self.pima.update_cnt(
                {"FRIB.OVERSAMPLE_MD:": "2", "FRIB.OVERSAMPLE_RT:": "2"}
            )

        # Average all spectral channels in each IF when splitting.
        self.pima.update_cnt({"SPLT.FRQ_MSEG:": str(self.pima.chan_number)})

        if scan_part == 1:
            self.pima.update_cnt(
                {"FRIB.1D_RESFRQ_PLOT:": "TXT", "FRIB.1D_RESTIM_PLOT:": "TXT"}
            )
        else:
            self.pima.update_cnt(
                {"FRIB.1D_RESFRQ_PLOT:": "NO", "FRIB.1D_RESTIM_PLOT:": "NO"}
            )

    def load_antab(self, antab_file=None):
        """Download ANTAB file and load calibration information to PIMA."""
        # Always download antab-file.
        if antab_file:
            self.antab = antab_file
        elif not self.antab:
            self._get_antab()

        # Try to load calibration information from ANTAB
        if self.antab and os.path.isfile(self.antab):
            try:
                self.pima.load_gains(self.antab)
                self.pima.load_tsys(self.antab)
                self.calibration_loaded = True
            except pypima.pima.Error:
                self.logger.warning("Could not load calibration information")
                self.calibration_loaded = False

    def _select_ref_sta(self, fri, ref_sta=None):
        """
        Select reference station for bandpass calibration.

        Parameters
        ----------
        fri : ``Fri`` object
            **PIMA** fringe fitting results as ``Fri`` object.
        ref_sta : str, optional

        """
        if ref_sta and ref_sta not in self.pima.station_list():
            self._error(f"Station {ref_sta} is not in station list")
        else:
            self.sta_ref = ref_sta

        if not self.sta_ref:
            snr_detecton = float(self.pima.cnt_params["FRIB.SNR_DETECTION:"])

            obs = fri.max_snr("RADIO-AS")

            if obs:
                if obs["SNR"] < snr_detecton:
                    self.logger.debug(
                        "SNR is too low on space baseline for \
    bandpass: %s",
                        obs["SNR"],
                    )
                else:
                    if obs["sta1"] == "RADIO-AS":
                        self.sta_ref = obs["sta2"]
                    elif obs["sta2"] == "RADIO-AS":
                        self.sta_ref = obs["sta1"]
            else:
                self.logger.info("No scans with RADIO-AS")

            if not self.sta_ref:
                obs = fri.max_snr()

                if not obs:
                    return False

                if obs["SNR"] < snr_detecton:
                    self.logger.debug("SNR is too low for bandpass: %s", obs["SNR"])
                else:
                    good_stations = ["ARECIBO", "GBT-VLBA", "EFLSBERG", "ATCA-104"]
                    for sta in good_stations:
                        if sta in (obs["sta1"], obs["sta2"]):
                            self.sta_ref = sta
                            break
                    if self.sta_ref is None:
                        self.sta_ref = obs["sta1"]

        if self.sta_ref:
            self.pima.update_cnt(
                {
                    "STA_REF:": self.sta_ref,
                    "BPS.SNR_MIN_ACCUM:": "5.5",
                    "BPS.SNR_MIN_FINE:": "5.5",
                }
            )
            self.logger.info("New reference station is %s", self.sta_ref)

            return True
        else:
            return False

    def _check_bad_autospec_obs(self):
        """Return set of observation numbers with bad autospectrum."""
        self.generate_autospectra()

        if not self.acta_files:
            return

        bad_obs_set = set()

        for obs in self.pima.observations:
            for sta in (obs.sta1, obs.sta2):
                if self.pima.cnt_params["POLAR:"] in ("RR", "LL"):
                    sta_pol = self.pima.cnt_params["POLAR:"]
                elif self.pima.cnt_params["POLAR:"] == "RL":
                    if sta == obs.sta1:
                        sta_pol = "RR"
                    else:
                        sta_pol = "LL"
                elif self.pima.cnt_params["POLAR:"] == "LR":
                    if sta == obs.sta1:
                        sta_pol = "LL"
                    else:
                        sta_pol = "RR"
                else:
                    raise RuntimeError(
                        "unsupported polar {}".format(self.pima.cnt_params["POLAR:"])
                    )

                try:
                    acta_file = self.acta_files[sta_pol, obs.time_code, sta]
                except KeyError as err:
                    self.logger.warning("no ACTA file for %s", err)
                    continue

                if np.median(acta_file.ampl) < 0.5:
                    self.logger.warning("Bad autospec for sta: %s obs: %s", sta, obs.obs)
                    bad_obs_set.add(obs.obs)

        return bad_obs_set

    def _auto_bpas(self, fringe_fit: bool = False) -> None:
        """Iterate over bandpass parameters and select the best case."""
        # log_bps_dict = {}
        snr_dict = {}

        for deg in range(1, 7):
            log_file = self.pima.bpas(
                params={
                    "BPS.DEG_AMP:": str(deg),
                    "BPS.DEG_PHS:": str(deg),
                    "PHASE_ACCEL_MIN:": "0",
                    "PHASE_ACCEL_MAX:": "0",
                    "FRIB.FINE_SEARCH:": "LSQ",
                }
            )
            log_file_deg = f"{log_file}_{deg}"
            os.rename(log_file, log_file_deg)

            bps_file = self.pima.cnt_params["BANDPASS_FILE:"]
            bps_file_deg = f"{bps_file}_{deg}"
            os.rename(bps_file, bps_file_deg)

            # log_bps_dict[log_file_deg] = bps_file_deg

            if fringe_fit:
                fri = Fri(
                    self.pima.fine(
                        params={
                            "PHASE_ACCEL_MIN:": "0",
                            "PHASE_ACCEL_MAX:": "0",
                            "FRIB.FINE_SEARCH:": "LSQ",
                            "BANDPASS_FILE:": bps_file_deg,
                        }
                    )
                )
                if not fri:
                    self._error("fringe fitting fails in _auto_bpas")
                fri.update_status(64)
                # snr_data = {}
                # for rec in fri:
                #     if rec['status'] == 'y':
                #         snr_data[rec['obs']] = rec['SNR']
                if fri.any_detections("RADIO-AS"):
                    snr_data = {
                        rec["obs"]: rec["SNR"]
                        for rec in fri
                        if rec["status"] == "y"
                        if rec["sta1"] == "RADIO-AS"
                    }
                else:
                    snr_data = {
                        rec["obs"]: rec["SNR"] for rec in fri if rec["status"] == "y"
                    }
            else:
                snr_data = bpas_log_snr_new(log_file_deg, mode="INIT")

            if snr_data:
                snr_dict[bps_file_deg] = snr_data

        if not snr_dict:
            self.logger.warning("could not get bpas_accum SNR from logs")
            self.pima.bpas()
        else:
            table = pd.DataFrame.from_records(
                list(snr_dict.values()), index=list(snr_dict.keys())
            )

            self.logger.debug("\n%s", str(table))

            norm_table = table / table.iloc[0]
            scores = norm_table.sum(axis=1)

            # best_bps_file = log_bps_dict[scores.idxmax()]
            best_bps_file = scores.idxmax()
            self.logger.info("Best bps is: %s", best_bps_file)

            self.pima.update_cnt({"BANDPASS_FILE:": best_bps_file})

    def _bandpass(
        self, bandpass_mode=None, ampl_bandpass=True, bandpass_var=0, bandpass_norm="IF"
    ):
        """
        Setup **PIMA** bandpass parameters and run ``bpas`` task.

        """
        bpas_params = {}

        if bandpass_var == 0:
            self.logger.warning("using bandpass parameters from TEMPLATE")
        elif bandpass_var == 1:
            bpas_params = {
                "BPS.MODE:": "FINE",
                "BPS.NOBS_ACCUM:": "8",
                "BPS.MSEG_ACCUM:": "1",
                "BPS.NOBS_FINE:": "12",
                "BPS.MINOBS_FINE:": "8",
                "BPS.MSEG_FINE:": "1",
                "BPS.SNR_MIN_ACCUM:": "200.0",
                "BPS.SNR_MIN_FINE:": "200.0",
                "BPS.AMPL_REJECT:": "0.4",
                "BPS.PHAS_REJECT:": "0.2",
                "BPS.INTRP_METHOD:": "SPLINE",
                "BPS.DEG_AMP:": "17",
                "BPS.DEG_PHS:": "11",
                "BPS.AMP_MIN:": "0.01",
                "BPS.NORML:": "IF",
                "BPS.SEFD_USE:": "NO",
            }
        elif bandpass_var == 2:
            bpas_params = {
                "BPS.MODE:": "FINE",
                "BPS.NOBS_ACCUM:": "8",
                "BPS.MSEG_ACCUM:": "1",
                "BPS.NOBS_FINE:": "12",
                "BPS.MINOBS_FINE:": "8",
                "BPS.MSEG_FINE:": "1",
                "BPS.SNR_MIN_ACCUM:": "50.0",
                "BPS.SNR_MIN_FINE:": "50.0",
                "BPS.AMPL_REJECT:": "0.4",
                "BPS.PHAS_REJECT:": "0.2",
                "BPS.INTRP_METHOD:": "SPLINE",
                "BPS.DEG_AMP:": "5",
                "BPS.DEG_PHS:": "5",
                "BPS.AMP_MIN:": "0.01",
                "BPS.NORML:": "IF",
                "BPS.SEFD_USE:": "NO",
            }
        elif bandpass_var == 3:
            mseg = self.pima.chan_number // 2
            min_snr = 5.0  # Could be tuned

            bpas_params = {
                "BPS.MODE:": "ACCUM",
                "BPS.NOBS_ACCUM:": "6",
                "BPS.MSEG_ACCUM:": mseg,
                "BPS.NOBS_FINE:": "12",
                "BPS.MINOBS_FINE:": "8",
                "BPS.MSEG_FINE:": mseg,
                "BPS.SNR_MIN_ACCUM:": min_snr,
                "BPS.SNR_MIN_FINE:": min_snr,
                "BPS.AMPL_REJECT:": "0.4",
                "BPS.PHAS_REJECT:": "0.2",
                "BPS.INTRP_METHOD:": "LINEAR",
                "BPS.DEG_AMP:": "0",
                "BPS.DEG_PHS:": "1",
                "BPS.AMP_MIN:": "0.01",
                "BPS.NORML:": "IF",
                "BPS.SEFD_USE:": "NO",
                "FRIB.SNR_DETECTION:": min_snr,
            }
        elif bandpass_var in (4, 5):
            mseg = 4
            min_snr = 5.1  # Could be tuned

            bpas_params = {
                "BPS.MODE:": "ACCUM",
                "BPS.NOBS_ACCUM:": "6",
                "BPS.MSEG_ACCUM:": mseg,
                "BPS.NOBS_FINE:": "12",
                "BPS.MINOBS_FINE:": "8",
                "BPS.MSEG_FINE:": mseg,
                "BPS.SNR_MIN_ACCUM:": min_snr,
                "BPS.SNR_MIN_FINE:": min_snr,
                "BPS.AMPL_REJECT:": "0.4",
                "BPS.PHAS_REJECT:": "0.2",
                "BPS.INTRP_METHOD:": "LEGENDRE",
                "BPS.DEG_AMP:": "5",
                "BPS.DEG_PHS:": "5",
                "BPS.AMP_MIN:": "0.01",
                "BPS.NORML:": "IF",
                "BPS.SEFD_USE:": "NO",
                "FRIB.SNR_DETECTION:": min_snr,
            }
        else:
            self._error(f"Unsupported bandpass_var {bandpass_var}")

        if bandpass_mode:
            bpas_params["BPS.MODE:"] = bandpass_mode

        if not ampl_bandpass:
            bpas_params["BPS.DEG_AMP:"] = "0"

        if bandpass_norm:
            bpas_params["BPS.NORML:"] = bandpass_norm
            if bandpass_norm == "NO":
                bpas_params["BPS.MODE:"] = "INIT"

        self.pima.update_cnt(bpas_params)

        try:
            if bandpass_var in (4, 5) and self.pima.cnt_params["BPS.MODE:"] == "ACCUM":
                self.logger.info("starting _auto_bpas")
                if bandpass_var == 4:
                    self._auto_bpas(fringe_fit=False)
                else:
                    self._auto_bpas(fringe_fit=True)
            else:
                bpas_params = {
                    "FRIB.FINE_SEARCH:": "LSQ",
                    "MKDB.FRINGE_ALGORITHM:": "LSQ",
                    "PHASE_ACCEL_MIN:": "0",
                    "PHASE_ACCEL_MAX:": "0",
                }

                self.pima.bpas(bpas_params)
        except pypima.pima.Error:
            self.logger.warning("continue without bandpass")
            self.pima.update_cnt({"BANDPASS_FILE:": "NO"})
            return False

        return True

    def fringe_fitting(
        self,
        bandpass=False,
        accel=False,
        bandpass_mode=None,
        ampl_bandpass=True,
        bandpass_var=0,
        bandpass_use=None,
        bandpass_norm="IF",
        bandpass_renorm=True,
        reference_station=None,
    ):
        """
        Perform a fringe fitting.

        Parameters
        ----------
        bandpass : bool, optional
            If ``True`` try to do a bandpass calibration. Default is ``False``.
        accel : bool, optional
            If ``True`` turn on a phase acceleration fitting.
        bandpass_mode : str, optional
            Set the ``BPS.MODE`` **PIMA** parameter.
        ampl_bandpass : bool, optional
            If ``True``, do the amplitude bandpass calibration. Set polynomial
            degree to zero otherwise.
        bandpass_var : int, optional
            Select predefined bandpass parameters.
        bandpass_use : str, optional
            Set ``BANDPASS_USE`` **PIMA** parameter.
        bandpass_norm : str, optional
            Set ``BPS.NORML`` **PIMA** parameter. This keywords specifies the
            way how the bandpass normalization is made.
        bandpass_renorm: bool, optional
            If ``True``, apply the bandpass renormalization factors for given
            intermediate frequencies using only a part of the bandwidth.
        reference_station : str, optional
            Reference station for bandpass calibration. If ``None`` select
            an optimal station for ground-space baselines.

        Returns
        -------
        fri : ``Fri`` object
            Fringe fitting results as ``Fri`` object.

        """
        if accel:
            self.pima.update_cnt({"FRIB.FINE_SEARCH:": "ACC"})
            if self.band == "l":
                self.pima.update_cnt(
                    {"PHASE_ACCEL_MIN:": "-1.D-13", "PHASE_ACCEL_MAX:": "1.D-13"}
                )
            elif self.band == "k":
                self.pima.update_cnt(
                    {"PHASE_ACCEL_MIN:": "-5.D-15", "PHASE_ACCEL_MAX:": "5.D-15"}
                )
            else:
                self.pima.update_cnt(
                    {"PHASE_ACCEL_MIN:": "-1.D-14", "PHASE_ACCEL_MAX:": "1.D-14"}
                )
        else:
            self.pima.update_cnt(
                {
                    "FRIB.FINE_SEARCH:": "LSQ",
                    "PHASE_ACCEL_MIN:": "0",
                    "PHASE_ACCEL_MAX:": "0",
                }
            )

        if bandpass_use:
            self.pima.update_cnt({"BANDPASS_USE:": bandpass_use})

            if bandpass_use == "NO":
                bandpass = False

        if bandpass and self.pima.chan_number > 512:
            self.logger.warning(
                "Too many spectral channels for bandpass: %s", self.pima.chan_number
            )
            bandpass = False

        polar = self.pima.cnt_params["POLAR:"]
        frq_grp = int(self.pima.cnt_params["FRQ_GRP:"])

        if bandpass:
            if bandpass_renorm:
                self.pima.update_cnt({"SPLT.BPASS_NRML_METHOD:": "WEIGHTED"})
            else:
                self.pima.update_cnt({"SPLT.BPASS_NRML_METHOD:": "NO"})

            # Update list of obs with bad autospectrum
            bad_obs = self._check_bad_autospec_obs()
            if bad_obs:
                self.bad_obs_set = bad_obs
            else:
                self.bad_obs_set.clear()

            # If bps-file already exists -- use it
            if (polar, frq_grp) in self.bpass_files and os.path.isfile(
                self.bpass_files[polar, frq_grp]
            ):
                self.pima.update_cnt(
                    {"BANDPASS_FILE:": self.bpass_files[polar, frq_grp]}
                )
            else:
                self.pima.mk_exclude_obs_file(self.bad_obs_set, "coarse")

                coarse_params = {
                    "FRIB.FINE_SEARCH:": "LSQ",
                    "MKDB.FRINGE_ALGORITHM:": "LSQ",
                    "PHASE_ACCEL_MIN:": "0",
                    "PHASE_ACCEL_MAX:": "0",
                }

                fri_file = self.pima.coarse(coarse_params)
                fri = Fri(fri_file)

                if self.run_id > 0 and fri:
                    self.db.fri2db(fri, self.pima.exper_info, self.run_id, nobps=True)

                # Exclude suspicious observations
                obs_list = []
                for rec in fri:
                    if abs(rec["rate"]) > 1e-11 or abs(rec["delay"]) > 1e-6:
                        obs_list.append(rec["obs"])

                self.pima.mk_exclude_obs_file(obs_list, "bpas")
                fri.remove_obs(obs_list)

                self.pima.update_cnt({"FRIB.SNR_DETECTION:": "5.2"})

                # Now auto select reference station
                if fri and self._select_ref_sta(fri, reference_station):
                    bandpass = self._bandpass(
                        bandpass_mode, ampl_bandpass, bandpass_var, bandpass_norm
                    )
                    self.bpass_files[polar, frq_grp] = self.pima.cnt_params[
                        "BANDPASS_FILE:"
                    ]
                else:
                    self.logger.info(
                        "skip bandpass due to absence of the " "useful scans"
                    )
                    bandpass = False
                    self.bpass_files[polar] = ""
                    self.pima.update_cnt({"BANDPASS_FILE:": "NO"})

        self.pima.mk_exclude_obs_file(self.bad_obs_set, "fine")

        if frq_grp > 1:
            fri_file = f"{self.exper}_{self.band}_{polar}_{frq_grp}.fri"
            frr_file = f"{self.exper}_{self.band}_{polar}_{frq_grp}.frr"
        else:
            fri_file = f"{self.exper}_{self.band}_{polar}.fri"
            frr_file = f"{self.exper}_{self.band}_{polar}.frr"

        fri_file = os.path.join(self.work_dir, fri_file)
        frr_file = os.path.join(self.work_dir, frr_file)
        self.pima.update_cnt({"FRINGE_FILE:": fri_file, "FRIRES_FILE:": frr_file})

        fri_file = self.pima.fine()
        self.fri = Fri(fri_file)
        self.fri.aux["bandpass"] = bandpass

        if not self.fri:
            self.logger.warning("PIMA fri-file is empty after fine")
        else:
            if self.pima.exper_info.sp_chann_num <= 128:
                ch_num = 64
            elif self.pima.exper_info.sp_chann_num == 256:
                ch_num = 256
            else:
                ch_num = 2048

            self.fri.update_status(ch_num)

            if self.run_id > 0:
                self.db.fri2db(self.fri, self.pima.exper_info, self.run_id)

        return self.fri

    def flag_edge_chann(self, number):
        """
        Flag `number` spectral channels at the edges of the bandpass. Must be
        called after ``load``.

        Parameters
        ----------
        number : int
            Number of spectral channels to flag.
        """
        chann_num = self.pima.chan_number
        mask = []

        if number < 0 or number >= chann_num / 2:
            self._error(f"invald number of channels to flag: {number}")
        elif number > 0:
            ind_frq = "1-{}".format(self.pima.exper_info.if_num)
            ind_chn1 = f"1-{number}"
            ind_chn2 = "{}-{}".format(chann_num - number + 1, chann_num)

            mask = [
                ("ALL", "ALL", ind_frq, ind_chn1, "OFF"),
                ("ALL", "ALL", ind_frq, ind_chn2, "OFF"),
            ]

        mask_gen_file = self.pima.mk_bpass_mask_gen(mask)
        mask_file = self.pima.set_mask_file(mask_gen_file)

        if mask_file:
            self.logger.info("Set %s as new mask file", mask_file)

    def split(self, source=None, average=False):
        """
        Split a multi-source uv data set into single-source data files.

        Parameters
        ----------
        source : string, optional
            Do split only for given source. By default split all sources in
            the experiment.
        average : bool, optional
            If ``True`` average data over full scan length.

        """
        # Delete old uv-fits remained from previous run
        exper_dir = self.pima.cnt_params["EXPER_DIR:"]
        sess_code = self.pima.cnt_params["SESS_CODE:"]
        pima_fits_dir = os.path.join(exper_dir, sess_code + "_uvs")

        if os.path.isdir(pima_fits_dir):
            shutil.rmtree(pima_fits_dir)

        if not self.calibration_loaded:
            self.logger.warning(
                "Could not do splitting due to absence of " "calibration information"
            )
            return

        if not self.fri.any_detections():
            self.logger.warning("No useful scans for splitting")
            return

        snr_detection = round(min(7.0, self.fri.min_detected_snr() - 0.05), 2)
        self.logger.info("Set FRIB.SNR_DETECTION to %s", snr_detection)
        split_params = {
            "FRIB.SNR_DETECTION:": f"{snr_detection:.2f}",
            "DEBUG_LEVEL:": "6",
        }

        # Exclude suspicious observations
        obs_list = self.fri.non_detections()
        for rec in self.fri:
            if (
                abs(rec["rate"]) > 1e-10
                or abs(rec["delay"]) > 1e-6
                or rec["duration"] < 30
            ):
                obs_list.append(rec["obs"])

        # Exclude very short observations
        ap_len = self.pima.ap_minmax[0]
        min_scan_len = float(self.pima.cnt_params["MIN_SCAN_LEN:"])
        for obs in self.pima.observations:
            if obs.ap_num * ap_len < min_scan_len:
                obs_list.append(obs.obs)

        self.pima.mk_exclude_obs_file(obs_list, "splt")

        if source:
            split_params["SPLT.SOU_NAME:"] = source
        else:
            split_params["SPLT.SOU_NAME:"] = "ALL"

        if average:
            time_segments = max([obs.ap_num for obs in self.pima.observations])
        else:
            time_segments = 1

        self.split_time_aver = time_segments * ap_len
        self.pima.split(tim_mseg=time_segments, params=split_params)

    def copy_uvfits(self, out_dir):
        """Copy calibrated uv-fits files from pima scratch dir to `out_dir`."""
        exper_dir = self.pima.cnt_params["EXPER_DIR:"]
        sess_code = self.pima.cnt_params["SESS_CODE:"]
        band = self.pima.cnt_params["BAND:"]
        polar = self.pima.cnt_params["POLAR:"]

        pima_fits_dir = os.path.join(exper_dir, sess_code + "_uvs")

        if not os.path.isdir(pima_fits_dir):
            return

        splt_sou_name = self.pima.cnt_params["SPLT.SOU_NAME:"]

        for source_names in self.pima.source_list:
            if splt_sou_name != "ALL" and splt_sou_name not in source_names:
                continue

            pima_fits_name = "{}_{}_uva.fits".format(source_names[1], band)
            pima_fits_path = os.path.join(pima_fits_dir, pima_fits_name)

            if not os.path.isfile(pima_fits_path):
                self.logger.warning('UV-FITS "%s" does not exists.', pima_fits_path)
                continue

            # Use B1950 name for output directory
            b1950_name = source_names[2]

            # Fix source names
            if b1950_name == "OJ287":
                b1950_name = "0851+202"

            out_fits_dir = os.path.join(out_dir, b1950_name)
            os.makedirs(out_fits_dir, exist_ok=True)

            out_fits_name = "{}_{}_{}_{}_{:04d}s_{}_uva.fits".format(
                b1950_name,
                self.exper,
                self.band.upper(),
                polar,
                round(self.split_time_aver),
                self.pima.exper_info.correlator_name,  # correlator name
            )

            if self.scan_part >= 1000:
                scan_part_base = (self.scan_part // 1000) * 1000
                out_fits_name = out_fits_name.replace(
                    "_uva", f"_ALT{scan_part_base}_uva"
                )

            out_fits_path = os.path.join(out_fits_dir, out_fits_name)

            self.logger.info("Copy %s to %s", pima_fits_path, out_fits_path)
            shutil.copy(pima_fits_path, out_fits_path)

            # Run `fits_to_radplot` only for averaged uv-fits
            if self.split_time_aver > 2:
                try:
                    pypima.pima.fits_to_txt(out_fits_path)
                except subprocess.SubprocessError:
                    self._error("fits_to_radplot failed")

                if self.run_id > 0:
                    with UVFits(out_fits_path) as uvfits_file:
                        self.db.uvfits2db(uvfits_file, b1950_name, self.run_id)

    def fringes2db(self):
        """Put fringe fitting information to the database."""
        if self.run_id > 0 and self.fri:
            self.db.fri2db(self.fri, self.pima.exper_info, self.run_id)

    def delete_uvfits(self):
        """Delete UV-FITS file."""
        # Delete FITS file in `data_dir` only
        if (
            isinstance(self.uv_fits, str)
            and os.path.isfile(self.uv_fits)
            and self.uv_fits.startswith(self.data_dir)
        ):
            os.remove(self.uv_fits)

            # Delete data directory if empty
            try:
                os.rmdir(self.data_dir)
            except OSError:
                pass

        # Delete staging directory
        staging_dir = self.pima.cnt_params["STAGING_DIR:"]
        if os.path.isdir(staging_dir):
            shutil.rmtree(staging_dir)

    def generate_autospectra(self, plot=False, out_dir=None, db=False) -> None:
        """
        Generate autocorrelation spectrum.

        Thist function generates autocorrelation spectrum for each station for each
        scan using ``acta`` **PIMA** task and fill `self.acta_files` dict.

        Parameters
        ----------
        plot : bool
            If ``True`` plot autospectra.
        out_dir : str
            Plot output directory.
        db : bool
            If ``True`` store autospectra to the database.

        """
        if self.acta_files:
            self.logger.debug("acta has already been called")
            return

        for polar in ("RR", "LL"):
            # Sometimes PIMA crashes on `acta` task
            try:
                file_list = self.pima.acta(params={"POLAR:": polar})
            except pypima.pima.Error:
                # Remove core dump file.
                if os.path.isfile("core"):
                    os.remove("core")

                return

            for file_name in file_list:
                acta_file = ActaFile(
                    file_name, polar, self.pima.exper_info.utc_minus_tai
                )

                sta = acta_file.header["station"]
                scan_name = acta_file.header["scan_name"]

                assert (polar, scan_name, sta) not in self.acta_files

                self.acta_files[polar, scan_name, sta] = acta_file

                if plot:
                    acta_file.plot(out_dir)

                if db:
                    self.db.autospec2db(acta_file)


def _download_it(url, buffer, max_retries=0, ftp_user=None):
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
