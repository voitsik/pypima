"""
Created on Tue Dec 10 17:21:29 2013

@author: Petr Voytsik
"""

import logging
import os.path
import shutil
import subprocess
from collections import namedtuple
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure


class Error(Exception):
    """Raised when PIMA error occurs."""

    def __init__(self, exper, band, msg):
        self.exper = exper
        self.band = band
        self.msg = msg

    def __str__(self):
        return f"{self.exper}({self.band}): {self.msg}"


class ExperInfo:
    """Some experiment information."""

    def __init__(self, exper_name, band, stt_file=None):
        self.exper = exper_name
        self.band = band

        self._data = {}

        if stt_file:
            self.update(stt_file)

    def update(self, stt_file):
        """
        Fill ExperInfo with info from PIMA stt-file
        """

        with open(stt_file) as fil:
            for line in fil:
                if line.startswith("Correlator_name:"):
                    self._data["correlator_name"] = line.split()[1]
                elif line.startswith("Number of frequencies:"):
                    self._data["if_num"] = int(line.split(":")[1])
                elif line.startswith("Number of frequency groups:"):
                    self._data["frq_grp"] = int(line.split(":")[1])
                elif line.startswith("Number of spectral channels:"):
                    self._data["sp_chann_num"] = int(line.split(":")[1])
                elif line.startswith("Number of time epochs:"):
                    self._data["time_epochs_num"] = int(line.split(":")[1])
                elif line.startswith("Number of scans:"):
                    self._data["scans_num"] = int(line.split(":")[1])
                elif line.startswith("Number of observations:"):
                    self._data["obs_num"] = int(line.split(":")[1])
                elif line.startswith("Total number of UV points:"):
                    self._data["uv_points_num"] = int(line.split(":")[1])
                elif line.startswith("Total number of used UV points:"):
                    self._data["uv_points_used_num"] = int(line.split(":")[1])
                elif line.startswith("Total number of deselected points:"):
                    self._data["deselected_points_num"] = int(line.split(":")[1])
                elif line.startswith(
                    "Number of cross-correl NO_AUTO_1 deselected points:"
                ):
                    no_auto1 = int(line.split(":")[1])
                elif line.startswith(
                    "Number of cross-correl NO_AUTO_2 deselected points:"
                ):
                    no_auto2 = int(line.split(":")[1])
                elif line.startswith("Accummulation period length_min:"):
                    acc_min = float(line.split(":")[1])
                elif line.startswith("Accummulation period length_max:"):
                    acc_max = float(line.split(":")[1])
                elif line.startswith("UTC_minus_TAI:"):
                    self._data["utc_minus_tai"] = timedelta(
                        seconds=float(line.split(":")[1])
                    )
                elif line.startswith("Experiment nominal start:"):
                    self._data["nominal_start"] = datetime.strptime(
                        line.split(":", 1)[1].strip()[:23], "%Y.%m.%d-%H:%M:%S.%f"
                    )
                elif line.startswith("Experiment nominal end:"):
                    self._data["nominal_end"] = datetime.strptime(
                        line.split(":", 1)[1].strip()[:23], "%Y.%m.%d-%H:%M:%S.%f"
                    )
                elif line.startswith("# Generated    at"):
                    cols = line.split()
                    if len(cols) == 4:
                        self._data["hostname"] = cols[3]
                    else:
                        self._data["hostname"] = ""
                elif line.startswith("#                 by PIMA"):
                    if "version" in line:
                        self._data["pima_version"] = line.split()[5]
                    else:
                        self._data["pima_version"] = line.split()[4]

        self._data["accum_length"] = (acc_min + acc_max) / 2.0
        self._data["no_auto_points_num"] = no_auto1 + no_auto2

    def __getitem__(self, key):
        return self._data[key]


# Some global variables
UVFILE_NAME_LEN = 128


class Pima:
    """
    Pima class is analog of pima_fringe.csh script.

    """

    def __init__(self, experiment_code, band, work_dir=None):
        # First, set common variables
        self.exper = experiment_code.lower()
        self.band = band.lower()

        self.logger = logging.getLogger(f"{self.exper}({self.band})")

        if not work_dir:
            self.work_dir = os.getcwd()
        else:
            self.work_dir = work_dir

        self.pima_dir = os.getenv("PIMA_DIR")

        if not self.pima_dir or not os.path.isdir(self.pima_dir):
            self._error(
                "Could not find PIMA directory. Please, set $PIMA_DIR "
                "environment variable."
            )

        self.pima_exec = os.path.join(self.pima_dir, "bin", "pima")

        if not os.path.isfile(self.pima_exec):
            self._error("Could not find pima executable. Check your PIMA installation!")

        # PIMA control file path
        self.cnt_file_name = f"{self.exper}_{self.band}_pima.cnt"
        self.cnt_file_name = os.path.join(self.work_dir, self.cnt_file_name)

        # Dictionary with all parameters from cnt-file
        self.cnt_params = {}
        self._update_cnt_params()
        self.exper_info = ExperInfo(self.exper, self.band)
        # Update exper_info if experiment already loaded.
        stt_file = os.path.join(
            self.cnt_params["EXPER_DIR:"], self.cnt_params["SESS_CODE:"] + ".stt"
        )
        if os.path.isfile(stt_file):
            self.exper_info.update(stt_file)

    def _update_cnt_params(self):
        """Read cnt-file and fill cnt_params dictionary."""
        self.cnt_params.clear()
        self.cnt_params["UV_FITS:"] = []

        with open(self.cnt_file_name) as cnt_file:
            for line in cnt_file:
                line = line.split("#")[0].strip()
                if len(line) < 8:
                    continue
                key, val = line.split(None, 1)

                # Special case
                if key == "UV_FITS:":
                    self.cnt_params[key].append(val)
                else:
                    self.cnt_params[key] = val

    def update_cnt(self, opts):
        """Update pima control file according to `opts` dictionary."""
        if not opts:
            return

        uv_fits_flag = False

        lines = []

        with open(self.cnt_file_name) as file:
            for line in file:
                key = line.split()[0].strip()

                if key in opts.keys():
                    val = opts[key]

                    # Skip empty values
                    if not val:
                        continue

                    # Modify UV_FITS only once
                    if key == "UV_FITS:":
                        if uv_fits_flag:
                            continue
                        else:
                            uv_fits_flag = True

                    # In case of many FITS-files
                    if key == "UV_FITS:" and isinstance(val, list):
                        line = "".join([f"{key:<10} {item}\n" for item in val])
                    else:
                        line = f"{key:<20} {val}\n"
                elif line.startswith("# Last update on"):
                    line = "# Last update on  {}\n".format(str(datetime.now()))

                lines.append(line)

        with open(self.cnt_file_name, "w") as file:
            file.write("".join(lines))

        self._update_cnt_params()

    def _exec(self, operation, options=None, log_name=None) -> int:
        """
        Execute PIMA binary.

        Parameters
        ----------
        operation : str
            ``PIMA`` task to execute.
        options : dict, optional
            Additional ``PIMA`` parameters.
        log_name : str, optional
            Name of the log file.

        Returns
        -------
        ret : int
            Return ``pima`` process exit code.

        """
        cmd_line = [self.pima_exec, self.cnt_file_name, operation]
        if options:
            opt_list = [
                item for key in options for item in (str(key), str(options[key]))
            ]
            cmd_line.extend(opt_list)

        if not log_name:
            log_name = os.path.join(
                self.work_dir, f"{self.exper}_{self.band}_{operation}.log"
            )

        with open(log_name, "w") as log:
            print(datetime.now(), end="\n\n", file=log)
            print("cmd: ", " ".join(cmd_line), end="\n\n", file=log, flush=True)
            self.logger.debug("execute: %s", " ".join(cmd_line))

            # Run `pima` with minimum priority
            with subprocess.Popen(
                cmd_line, stdout=log, universal_newlines=True
            ) as proc:
                os.setpriority(os.PRIO_PROCESS, proc.pid, 19)
                os.sched_setscheduler(proc.pid, os.SCHED_BATCH, os.sched_param(0))
                ret = proc.wait()

            print("", file=log)
            print(datetime.now(), file=log)

        return ret

    def _error(self, msg):
        """Raise pima.Error exception."""
        self.logger.error(msg)
        raise Error(self.exper, self.band, msg)

    def load(self):
        """
        Run pima load.

        This function removes all existing index files of the current session and
        the runs ``load`` task.

        """
        exper_dir_path = Path(self.cnt_params["EXPER_DIR:"])

        # Delete existing PIMA auxiliary files
        auxiliary_files = exper_dir_path.glob(self.cnt_params["SESS_CODE:"] + "*")

        for aux_file in auxiliary_files:
            if aux_file.is_file():
                aux_file.unlink()
            elif aux_file.is_dir():
                shutil.rmtree(aux_file)

        opts = {"BANDPASS_FILE:": "NO", "POLARCAL_FILE:": "NO"}
        ret = self._exec("load", options=opts)
        if ret:
            self._error(f"load failed with code {ret}")

        stt_file = exper_dir_path / (self.cnt_params["SESS_CODE:"] + ".stt")
        self.exper_info.update(stt_file)

        self.logger.info("load ok")

    def coarse(self, params=None):
        """
        Do coarse fringe fitting.

        This function runs 'pima frib' with spetial parameters:

        1. Disable bandpass

        2. Use fast algorithm for fringe search

        3. Set fri-file name

        Parameters
        ----------
        params : dict, optional
            Dictionary with optional pima parameters.

        Returns
        -------
        fri_file : str
            Name of the fri-file

        """
        if params and "POLAR:" in params:
            polar = params["POLAR:"]
        else:
            polar = self.cnt_params["POLAR:"]

        frq_grp = int(self.cnt_params["FRQ_GRP:"])

        if frq_grp > 1:
            name_base = f"{self.exper}_{self.band}_{polar}_{frq_grp}"
        else:
            name_base = f"{self.exper}_{self.band}_{polar}"

        log_name = f"{name_base}_coarse.log"
        log_name = os.path.join(self.work_dir, log_name)
        fri_file = f"{name_base}_nobps.fri"
        fri_file = os.path.join(self.work_dir, fri_file)
        frr_file = f"{name_base}_nobps.frr"
        frr_file = os.path.join(self.work_dir, frr_file)
        exc_obs_file = f"{name_base}_coarse_obs.exc"
        exc_obs_file = os.path.join(self.work_dir, exc_obs_file)

        if os.path.isfile(fri_file):
            os.remove(fri_file)

        if os.path.isfile(frr_file):
            os.remove(frr_file)

        opts = {
            "FRINGE_FILE:": fri_file,
            "FRIRES_FILE:": frr_file,
            "BANDPASS_USE:": "NO",
            "BANDPASS_FILE:": "NO",
            "POLARCAL_FILE:": "NO",
        }

        if os.path.isfile(exc_obs_file):
            opts["EXCLUDE_OBS_FILE:"] = exc_obs_file

        if params:
            opts.update(params)

        ret = self._exec("frib", opts, log_name)

        if ret:
            self._error(f"coarse failed with code {ret}")

        self.logger.info("coarse ok")

        return fri_file

    def fine(self, params=None):
        """
        Do file fringe fitting.

        This function runs 'pima frib' with parameters from cnt-file plus user
        defined params.

        Parameters
        ----------
        params : dict, optional
            Dictionary with optional pima parameters.

        Returns
        -------
        fri_file : str
            Name of the fri-file.

        """
        if params and "POLAR:" in params:
            polar = params["POLAR:"]
        else:
            polar = self.cnt_params["POLAR:"]

        log_name = f"{self.exper}_{self.band}_{polar}_fine.log"
        log_name = os.path.join(self.work_dir, log_name)
        exc_obs_file = f"{self.exper}_{self.band}_{polar}_fine_obs.exc"
        exc_obs_file = os.path.join(self.work_dir, exc_obs_file)

        if params and "FRINGE_FILE:" in params:
            fri_file = params["FRINGE_FILE:"]
        else:
            fri_file = self.cnt_params["FRINGE_FILE:"]

        if params and "FRIRES_FILE:" in params:
            frr_file = params["FRIRES_FILE:"]
        else:
            frr_file = self.cnt_params["FRIRES_FILE:"]

        if os.path.isfile(fri_file):
            os.remove(fri_file)

        if os.path.isfile(frr_file):
            os.remove(frr_file)

        opts = {}

        if os.path.isfile(exc_obs_file):
            opts["EXCLUDE_OBS_FILE:"] = exc_obs_file

        if params:
            opts.update(params)

        ret = self._exec("frib", opts, log_name=log_name)

        if ret:
            self._error(f"fine failed with code {ret}")

        self.logger.info("fine ok")

        return fri_file

    def bpas(self, params=None) -> str:
        """
        Do bandpass calibration.

        This function runs 'pima bpas'.

        Parameters
        ----------
        params : dict, optional
            Dictionary with optional pima parameters.

        Returns
        -------
        log_file : str
            Path to the bandpass log file.

        """
        if params and "POLAR:" in params:
            polar = params["POLAR:"]
        else:
            polar = self.cnt_params["POLAR:"]

        if polar in ["I", "RPL"]:
            polar = "RR"

        frq_grp = int(self.cnt_params["FRQ_GRP:"])

        if frq_grp > 1:
            name_base = f"{self.exper}_{self.band}_{polar}_{frq_grp}"
        else:
            name_base = f"{self.exper}_{self.band}_{polar}"

        fri_file = f"{name_base}_nobps.fri"
        fri_file = os.path.join(self.work_dir, fri_file)
        log_file = f"{name_base}_bps.log"
        log_file = os.path.join(self.work_dir, log_file)
        exc_obs_file = f"{name_base}_bpas_obs.exc"
        exc_obs_file = os.path.join(self.work_dir, exc_obs_file)

        bps_file = f"{name_base}.bps"
        bps_file = os.path.join(self.work_dir, bps_file)
        self.update_cnt({"BANDPASS_FILE:": bps_file})

        opts = {
            "FRINGE_FILE:": fri_file,
            "DEBUG_LEVEL:": "3",
            # "PHASE_ACCEL_MIN:": "0",
            # "PHASE_ACCEL_MAX:": "0",
            # "FRIB.FINE_SEARCH:": "LSQ",
            "POLAR:": polar,
        }

        if os.path.isfile(exc_obs_file):
            opts["EXCLUDE_OBS_FILE:"] = exc_obs_file

        if params:
            opts.update(params)

        ret = self._exec("bpas", opts, log_file)
        if ret:
            self._error(f"bpas failed with code {ret}")

        self.logger.info("bpas ok")

        return log_file

    def split(self, tim_mseg=1, params=None):
        """
        Do SPLIT.

        This function runs 'pima splt' with given integration time.

        Parameters
        ----------
        tim_mseg : int, optional
            Number of time segments to integrate.

        params : dict, optional
            Dictionary with optional pima parameters.

        """
        if params and "POLAR:" in params:
            polar = params["POLAR:"]
        else:
            polar = self.cnt_params["POLAR:"]

        exc_obs_file = f"{self.exper}_{self.band}_{polar}_splt_obs.exc"
        exc_obs_file = os.path.join(self.work_dir, exc_obs_file)
        log_file = f"{self.exper}_{self.band}_{polar}_splt.log"
        log_file = os.path.join(self.work_dir, log_file)

        opts = {"SPLT.TIM_MSEG:": str(tim_mseg)}

        if os.path.isfile(exc_obs_file):
            opts["EXCLUDE_OBS_FILE:"] = exc_obs_file

        if params:
            opts.update(params)

        ret = self._exec("splt", opts, log_file)

        if ret:
            self._error(f"splt failed with code {ret}")

        self.logger.info("split ok")

    def load_gains(self, gain_file, params=None) -> None:
        """
        Load antenna gains from the `gain_file`.

        Parameters
        ----------
        gain_file : str
            Name of the AIPS ANTAB file with gain information.

        params : dict, optional
            Dictionary with optional pima parameters.

        """
        opts = {"evn_gain": gain_file, "DEBUG_LEVEL:": "6"}

        if params:
            opts.update(params)

        log_file = f"{self.exper}_{self.band}_gain.log"
        log_file = os.path.join(self.work_dir, log_file)

        ret = self._exec("gean", opts, log_file)

        if ret:
            self._error(f"evn_gain failed with code {ret}")

        self.logger.info("evn_gain ok")

    def load_tsys(self, tsys_file, params=None) -> None:
        """
        Load Tsys from the `tsys_file`.

        Parameters
        ----------
        tsys_file : str
            Name of the AIPS ANTAB file with Tsys data.

        params : dict, optional
            Dictionary with optional pima parameters.

        """
        opts = {"vlba_log_file": tsys_file, "DEBUG_LEVEL:": "6"}

        if params:
            opts.update(params)

        log_file = f"{self.exper}_{self.band}_tsys.log"
        log_file = os.path.join(self.work_dir, log_file)

        ret = self._exec("gean", opts, log_file)

        if ret:
            self._error(f"vlba_log_file failed with code {ret}")

        self.logger.info("vlba_log_file ok")

    def acta(self, params=None):
        """
        Run `acta` pima task for autocorrelation spectrum generation.

        Parameters
        ----------
        params : dict, optional
            Dictionary with optional pima parameters.

        Returns
        -------
        file_list : list
            List of the names of the generated files.

        """
        if params and "POLAR:" in params:
            polar = params["POLAR:"]
        else:
            polar = self.cnt_params["POLAR:"]

        opts = {"DEBUG_LEVEL:": "2"}

        if params:
            opts.update(params)

        log_file = f"{self.exper}_{self.band}_{polar}_acta.log"
        log_file = os.path.join(self.work_dir, log_file)

        ret = self._exec("acta", opts, log_file)

        if ret:
            self._error(f"acta failed with code {ret}")
        else:
            self.logger.info("acta ok")

        file_list = []
        with open(log_file) as fil:
            for line in fil:
                if line.startswith("PIMA_ACTA created file:"):
                    file_name = line.split(":")[1].strip()
                    file_list.append(file_name)

        return file_list

    # Additional useful utilites
    def set_polar(self, polar):
        """
        Set polarization in control file.

        Parameters
        ----------
        polar : str
            Polarization code.

        """
        polar = polar.upper()

        if polar not in ("RR", "RL", "LR", "LL", "I"):
            self._error(f"Wrong polarization: {polar}")

        self.logger.info("Set polarization to %s", polar)
        self.update_cnt({"POLAR:": polar, "SPLT.POLAR:": polar})

    def set_frq_grp(self, frq_grp):
        """
        Set frequency group in control file.

        Parameters
        ----------
        frq_grp : int
            Frequency group number.

        """
        if frq_grp < 1 or frq_grp > self.exper_info["frq_grp"]:
            self._error(f"Invalid frequency group {frq_grp}")

        self.logger.info("Set frequency group to %s", frq_grp)
        self.update_cnt({"FRQ_GRP:": frq_grp})

    @property
    def ap_minmax(self):
        """
        Return minimum and maximum accummulation periods in experiment.

        """
        ap_min = ap_max = 0
        stt_file = os.path.join(
            self.cnt_params["EXPER_DIR:"], self.cnt_params["SESS_CODE:"] + ".stt"
        )

        if os.path.isfile(stt_file):
            with open(stt_file) as fil:
                for line in fil:
                    if line.startswith("Accummulation period length_min"):
                        ap_min = float(line.split()[3])
                    elif line.startswith("Accummulation period length_max"):
                        ap_max = float(line.split()[3])

        return ap_min, ap_max

    @property
    def number_of_deselected_points(self):
        """
        Return total number of deselected points

        """
        return self.exper_info["deselected_points_num"]

    def station_list(self, ivs_name=True):
        """
        Return a list of the station names.

        This function returns list of station names which participated in the
        experiment.

        Parameters
        ----------
        ivs_name : bool, optional
            If `ivs_name` is False this function returns list of 2-letter
            station codes instead of IVS names.

        Returns
        -------
        list
            List of the station names.

        """
        sta_l = []
        sta_file = os.path.join(
            self.cnt_params["EXPER_DIR:"], self.cnt_params["SESS_CODE:"] + ".sta"
        )

        if os.path.isfile(sta_file):
            with open(sta_file) as fil:
                for line in fil:
                    toks = line.split()
                    if ivs_name:
                        sta_l.append(toks[3])
                    else:
                        sta_l.append(toks[5])

        return sta_l

    @property
    def source_list(self):
        """
        Return a list of the source names.

        Returns
        -------
        out : list
            Each item in the out is a tuple of 3 names: IVS, J2000, B1950.

        """
        sou_list = []
        sou_file = os.path.join(
            self.cnt_params["EXPER_DIR:"], self.cnt_params["SESS_CODE:"] + ".sou"
        )

        if os.path.isfile(sou_file):
            with open(sou_file) as fil:
                for line in fil:
                    toks = line.split()
                    sou_list.append((toks[2], toks[3], toks[4]))

        return sou_list

    def source_dist(self):
        """
        Return distance between correlator phase center and source position
        from catalog.

        """
        dist = {}
        sou_file = os.path.join(
            self.cnt_params["EXPER_DIR:"], self.cnt_params["SESS_CODE:"] + ".sou"
        )

        if os.path.isfile(sou_file):
            with open(sou_file) as fil:
                for line in fil:
                    toks = line.split()
                    try:
                        dist[toks[2]] = float(toks[12])
                    except ValueError:
                        dist[toks[2]] = 9999999.9

        return dist

    @property
    def obs_number(self):
        """
        Return number of the observations in the experiment.

        """
        return self.exper_info["obs_num"]

    @property
    def chan_number(self):
        """
        Return number of the spectral channels in uv-data.

        """
        return self.exper_info["sp_chann_num"]

    @property
    def frequencies(self):
        """
        Return list of frequencies used in the experiment.

        Returns
        -------
        freqs : list
            The function returns a list of named tuples.

        """
        Freq = namedtuple("Freq", "freq band_width chan_width side_band")
        freqs = []

        frq_file = os.path.join(
            self.cnt_params["EXPER_DIR:"], self.cnt_params["SESS_CODE:"] + ".frq"
        )

        if os.path.isfile(frq_file):
            with open(frq_file) as fil:
                for line in fil:
                    if len(line) < 16:
                        break
                    if line.startswith("Ind_grp:"):
                        toks = line.split()
                        freq = Freq(
                            freq=float(toks[6]),
                            band_width=float(toks[8]),
                            chan_width=float(toks[10]),
                            side_band=int(toks[12]),
                        )
                        freqs.append(freq)

        return freqs

    @property
    def observations(self):
        """
        Return list of observations with some information.

        """
        Obs = namedtuple(
            "Obs", "obs scan time_code source sta1 sta2 start_time stop_time ap_num"
        )
        obs_list = []

        obs_file = os.path.join(
            self.cnt_params["EXPER_DIR:"], self.cnt_params["SESS_CODE:"] + ".obs"
        )

        if os.path.isfile(obs_file):
            with open(obs_file) as file:
                for line in file:
                    toks = line.split()

                    if len(toks) < 13:
                        continue

                    try:
                        ap_num = int(toks[12])
                    except ValueError as err:
                        self.logger.warning(
                            "could not get ap num for obs %s: %s", toks[0], err
                        )
                        ap_num = 0

                    obs = Obs(
                        obs=int(toks[0]),
                        scan=int(toks[1]),
                        time_code=toks[2],
                        source=toks[8],
                        sta1=toks[9],
                        sta2=toks[10],
                        start_time=datetime.strptime(toks[4], "%Y.%m.%d-%H:%M:%S.%f,"),
                        stop_time=datetime.strptime(toks[5], "%Y.%m.%d-%H:%M:%S.%f"),
                        ap_num=ap_num,
                    )

                    obs_list.append(obs)

        return obs_list

    def clock_model(self):
        """
        Read clock model components from the PIMA mdc file and retrun
        they as list of records.

        Returns
        -------
        clock_model : list


        """
        clock_model = []
        mdc_file = os.path.join(
            self.cnt_params["EXPER_DIR:"], self.cnt_params["SESS_CODE:"] + ".mdc"
        )

        if os.path.isfile(mdc_file):
            with open(mdc_file) as fil:
                for line in fil:
                    if not line.startswith("CLOCK_MODEL"):
                        continue

                    cols = line.strip().split()
                    if len(cols) != 19:
                        continue

                    sta = cols[2]
                    time = datetime.strptime(cols[6], "%Y.%m.%d-%H:%M:%S.%f")
                    try:
                        clock_offset = float(cols[8].replace("D", "e"))
                    except ValueError:
                        clock_offset = 0.0

                    try:
                        clock_rate = float(cols[10].replace("D", "e"))
                    except ValueError:
                        clock_rate = 0.0

                    try:
                        group_delay = float(cols[16].replace("D", "e"))
                    except ValueError:
                        group_delay = 0.0

                    try:
                        delay_rate = float(cols[18].replace("D", "e"))
                    except ValueError:
                        delay_rate = 0.0

                    clock_model.append(
                        (sta, time, clock_offset, clock_rate, group_delay, delay_rate)
                    )

        return clock_model

    def mk_exclude_obs_file(self, obs_list, suffix, polar=None):
        """
        Create ``EXCLUDE_OBS_FILE`` file using list of the observation indices.
        If `obs_list` is empty delete the file.

        Parameters
        ----------
        obs_list : list
            List of observation indices.
        suffix : str
            ``bpas`` or ``splt``.
        polar : str, optional
            Create file for given polarization.

        Returns
        -------
        exc_obs_file : srt
            Return name of the generated file or ``NO`` if `obs_list` is empty.

        """
        if not polar:
            polar = self.cnt_params["POLAR:"]

        frq_grp = int(self.cnt_params["FRQ_GRP:"])

        if frq_grp > 1:
            name_base = f"{self.exper}_{self.band}_{polar}_{frq_grp}"
        else:
            name_base = f"{self.exper}_{self.band}_{polar}"

        exc_obs_file = f"{name_base}_{suffix}_obs.exc"
        exc_obs_file = os.path.join(self.work_dir, exc_obs_file)

        if obs_list:
            with open(exc_obs_file, "w") as file:
                for obs in sorted(set(obs_list)):
                    print(obs, file=file)
        else:
            if os.path.exists(exc_obs_file):
                os.remove(exc_obs_file)

            exc_obs_file = "NO"

        return exc_obs_file

    def set_mask_file(self, mask_gen_file):
        """
        Create mask file based on `mask_gen_file` and update control file.

        Prameters
        ---------
        mask_gen_file : str
            **PIMA** mask definition file name. If empty or ``None`` set
            ``BANDPASS_MASK_FILE`` to ``NO`` and return ``None``.

        Returns
        -------
        mask_file : str
            Name of the generated mask file.

        """
        mask_file = None

        if mask_gen_file and os.path.isfile(mask_gen_file):
            mask_file = os.path.join(self.work_dir, f"{self.exper}_{self.band}.mask")
            if os.path.exists(mask_file):
                os.remove(mask_file)

            self.update_cnt({"BANDPASS_MASK_FILE:": mask_file})

            options = {"mask_gen": mask_gen_file}
            self._exec("bmge", options)
        else:
            self.update_cnt({"BANDPASS_MASK_FILE:": "NO"})

        return mask_file

    def mk_bpass_mask_gen(self, params):
        """
        Create PIMA ``BPASS_MASK_GEN`` file.

        Parameters
        ----------
        params : list
            List of tuples with parameters. Each tuple represents one row in
            mask gen file. If `params` empty or ``None`` remove
            ``BPASS_MASK_GEN`` file if exists.

        """
        mask_gen_file = os.path.join(
            self.work_dir, f"{self.exper}_{self.band}_mask.gen"
        )

        if params:
            with open(mask_gen_file, "w") as file:
                print("# PIMA BPASS_MASK_GEN  v 0.90 2009.02.05", file=file)
                print("#", file=file)
                print(
                    "#  Control for bandpass mask generation for \
experiment {}".format(
                        self.exper
                    ),
                    file=file,
                )
                print("#", file=file)
                print(f"#  Created on {datetime.now()}", file=file)
                print("#", file=file)
                print("BOTH   ALL:   ON", file=file)
                print("#", file=file)

                for row in params:
                    print(
                        "{} STA: {} IND_FRQ: {} IND_CHN: {} {}".format(*row), file=file
                    )

                print("#", file=file)
        else:
            if os.path.exists(mask_gen_file):
                os.remove(mask_gen_file)
            mask_gen_file = None

        return mask_gen_file


def fits_to_txt(fits_file: str) -> str:
    """
    Dump visibilites from the uvf file using the fits_to_radplot utility.

    ``fits_to_radplot`` reads visibility data from the input file and write
    amplitude, phase, and weight for each baseline to the text file.

    Parameters
    ----------
    fits_file : str
        Name of the input FITS file. FITS file must be in uvf format.

    Returns
    -------
    out : str
        Returns the stdout of the ``fits_to_radplot`` utility.

    """
    fits_file_path = Path(fits_file)

    if not fits_file_path.is_file():
        raise FileNotFoundError(f"file {fits_file_path} does not exist")

    suffix = fits_file_path.suffix

    if suffix:
        txt_file = fits_file.replace(suffix, ".txt")
    else:
        txt_file = fits_file + ".txt"

    cmd_line = ["fits_to_radplot", "-o", txt_file, fits_file]
    out = subprocess.check_output(cmd_line, universal_newlines=True)

    return out


class ActaFile:
    """
    This class represents PIMA ``ACTA`` file.

    """

    def __init__(self, input_file_name, polar=None, utc_tai=None):
        """
        Parameters
        ----------
        input_file_name : str
            Name of the text file generated by PIMA ``acta`` task.
        polar : str
            Polarization
        utc_tai : timedelta obj
            UTC-TAI as timedelta object.

        """
        self._freq = []
        self._if = []
        self._channel = []
        self._ampl = []
        self._header = {"polar": polar}

        with open(input_file_name) as file:
            magic = file.readline().strip()
            if magic != "# ACTA Output.  Format version of 2014.04.19":
                logging.error("File %s has invalid format")
                return

            for line in file:
                cols = line.split()

                if not cols:
                    continue

                if len(cols) == 3 and cols[0] == "#":
                    if cols[1] == "Experiment:":
                        self._header["experiment"] = cols[2]
                    elif cols[1] == "Station:":
                        self._header["station"] = cols[2]
                    elif cols[1] == "Scan_name:":
                        self.header["scan_name"] = cols[2]
                        logging.debug("scan_name = %s", cols[2])
                    elif cols[1] == "Scan_index:":
                        self.header["scan"] = int(cols[2])
                    elif cols[1] == "Observation_index:":
                        self.header["obs"] = int(cols[2])
                    elif cols[1] == "Start_date:":
                        self._header["start_date"] = datetime.strptime(
                            cols[2][:23], "%Y.%m.%d-%H:%M:%S.%f"
                        )
                        if utc_tai:
                            self._header["start_date"] += utc_tai
                    elif cols[1] == "Stop_date:":
                        self._header["stop_date"] = datetime.strptime(
                            cols[2][:23], "%Y.%m.%d-%H:%M:%S.%f"
                        )
                        if utc_tai:
                            self._header["stop_date"] += utc_tai
                    elif cols[1] == "Number_of_points:":
                        self.header["num_of_points"] = int(cols[2])
                elif cols[0] == "ACRL":
                    self._if.append(int(cols[2]))
                    self._channel.append(int(cols[4]))
                    self._freq.append(1e-6 * float(cols[6]))
                    self._ampl.append(float(cols[9]))

    @property
    def header(self):
        """
        Return dictionary with header information.

        """
        return self._header

    @property
    def if_num(self):
        """
        Return list of the IF numbers.

        """
        return self._if

    @property
    def channel(self):
        """
        Return list of the channel numbers,

        """
        return self._channel

    @property
    def freq(self):
        """
        Return list of the frequencies of the autospectrum.

        """
        return self._freq

    @property
    def ampl(self):
        """
        Return list of the amplitudes of the autospectrum.
        """
        return self._ampl

    def plot(self, out_dir_base: str, out_format: str = "pdf") -> None:
        """
        Plot autocorrelation spectrum to PDF file.

        Parameters
        ----------
        out_dir_base : str
            Base output directory.
        out_format : str, optional
            Format of the output file. The default is 'pdf'.

        """
        fig = Figure()
        FigureCanvas(fig)
        ax = fig.add_subplot(111)
        ax.set_ymargin(0.05)
        ax.xaxis.set_ticks(np.arange(-16, 16 + 1, 4))

        date = self.header["start_date"]
        date_str = date.strftime("%Y-%m-%d %H:%M:%S")

        ax.set_title(
            "{} - {} - {} - {}".format(
                self.header["station"],
                self.header["experiment"],
                self.header["polar"],
                date_str,
            )
        )
        ax.set_xlabel("Frequency (MHz)")
        ax.set_ylabel("Amplitude")
        ax.grid(True)
        central_freq = np.mean(self.freq)
        logging.debug("plot_autospectra: central_freq = %s", central_freq)
        freq = np.asarray(self.freq) - central_freq
        ax.plot(freq, self.ampl, marker="o", ms=2)
        ax.set_xlim(-16, 16)

        date_str = date.strftime("%Y%m%dT%H%M")
        out_file = "AUTOSPEC_{}_{}_{}_{}.{}".format(
            date_str,
            self.header["experiment"],
            self.header["polar"],
            self.header["station"],
            out_format,
        )

        out_dir = os.path.join(out_dir_base, self.header["experiment"])
        os.makedirs(out_dir, exist_ok=True)
        out_file = os.path.join(out_dir, out_file)
        fig.savefig(out_file, format=out_format)


class Text1D:
    """
    This class represents PIMA ``1D text table`` file.

    """

    def __init__(self, file_name=None):
        """
        Parameters
        ----------
        file_name : str
            Input file name.

        """
        self.header = {}
        self._ax1 = []
        self._ax2 = []

        if file_name:
            self.load_file(file_name)

    def load_file(self, file_name):
        """
        Read and parse file in 1D text table format.

        Parameters
        ----------
        file_name : str
            Input file name.

        """
        self.header.clear()
        self._ax1.clear()
        self._ax2.clear()

        with open(file_name) as file:
            magic = file.readline().strip()

            if magic != "# 1D text table.  Format version of 2012.12.30":
                logging.error("File %s has invalid format", file_name)
                return

            for line in file:
                if line[0] == "#":
                    continue

                cols = line.split(":", 1)
                if len(cols) != 2:
                    continue

                key = cols[0].strip()
                val = cols[1].strip()

                # print(key, val)

                if key == "POINT":
                    cols = val.split()
                    self._ax1.append(float(cols[1]))
                    self._ax2.append(float(cols[2]))
                else:
                    self.header[key] = val

    @property
    def axis1_data(self):
        """
        Return axis 1 array.

        """
        return self._ax1

    @property
    def axis2_data(self):
        """
        Return axis 2 array.

        """
        return self._ax2


def acta_plot(input_file: str, output_file: str) -> str:
    """Run ``acta_plot``."""
    if not os.path.isfile(input_file):
        raise FileNotFoundError(f"file {input_file} does not exist")

    cmd_line = ["acta_plot", input_file, output_file]
    out = subprocess.check_output(cmd_line, universal_newlines=True)

    return out


def bpas_log_snr_new(file_name: str, mode: str = "INIT"):
    """
    Retrieve ``new`` SNR values from bandpass log file.

    Parameters
    ----------
    file_name : str
        Log file name.
    mode : str, optional
        Bandpass stage. The default is 'INIT' that is all stages.

    Returns
    -------
    res : dict
        Return dictionary {obs: snr}.

    """
    res = {}
    in_bpass = False

    bpas_mode_line = f"PIMA_BPASS_{mode}"

    with open(file_name) as file:
        for line in file:
            if bpas_mode_line in line:
                in_bpass = True
                continue

            if in_bpass and line.startswith("Obs:"):
                cols = line.split()

                obs = int(cols[1])
                snr_new_ind = cols.index("new:") + 1

                try:
                    res[obs] = float(cols[snr_new_ind])
                except ValueError:
                    res[obs] = 0.0

    return res
