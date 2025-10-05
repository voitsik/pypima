"""
Module to run PIMA tasks and handle PIMA-related files.

Created on Tue Dec 10 17:21:29 2013

@author: Petr Voytsik
"""

import logging
import os
import os.path
import shutil
import subprocess
from collections import namedtuple
from dataclasses import InitVar, dataclass, field
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any, NamedTuple, NoReturn

import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

logger = logging.getLogger(__name__)

# Some global variables
MAX_UVFILE_NAME_LEN = 128


class Error(Exception):
    """Raised when PIMA error occurs."""

    def __init__(self, exper, band, msg):
        self.exper = exper
        self.band = band
        self.msg = msg

    def __str__(self):
        return f"{self.exper}({self.band}): {self.msg}"


@dataclass
class ExperInfo:
    """Some experiment information."""

    exper: str
    band: str
    stt_file: InitVar[str] = ""

    correlator_name: str = field(init=False)
    if_num: int = field(init=False)
    frq_grp: int = field(init=False)
    sp_chann_num: int = field(init=False)
    time_epochs_num: int = field(init=False)
    scans_num: int = field(init=False)
    obs_num: int = field(init=False)
    uv_points_num: int = field(init=False)
    uv_points_used_num: int = field(init=False)
    deselected_points_num: int = field(init=False)
    utc_minus_tai: timedelta = field(init=False)
    nominal_start: datetime = field(init=False)
    nominal_end: datetime = field(init=False)
    hostname: str = field(init=False)
    pima_version: str = field(init=False)
    accum_length: float = field(init=False)
    no_auto_points_num: int = field(init=False)

    def __post_init__(self, stt_file: str) -> None:
        if stt_file:
            self.update(stt_file)

    def update(self, stt_file: str) -> None:
        """Fill ExperInfo with info from PIMA stt-file."""
        acc_min = 0.0
        acc_max = 0.0
        no_auto1 = 0
        no_auto2 = 0

        with open(stt_file) as fil:
            for line in fil:
                if line.startswith("Correlator_name:"):
                    self.correlator_name = line.split()[1]
                elif line.startswith("Number of frequencies:"):
                    self.if_num = int(line.split(":")[1])
                elif line.startswith("Number of frequency groups:"):
                    self.frq_grp = int(line.split(":")[1])
                elif line.startswith("Number of spectral channels:"):
                    self.sp_chann_num = int(line.split(":")[1])
                elif line.startswith("Number of time epochs:"):
                    self.time_epochs_num = int(line.split(":")[1])
                elif line.startswith("Number of scans:"):
                    self.scans_num = int(line.split(":")[1])
                elif line.startswith("Number of observations:"):
                    self.obs_num = int(line.split(":")[1])
                elif line.startswith("Total number of UV points:"):
                    self.uv_points_num = int(line.split(":")[1])
                elif line.startswith("Total number of used UV points:"):
                    self.uv_points_used_num = int(line.split(":")[1])
                elif line.startswith("Total number of deselected points:"):
                    self.deselected_points_num = int(line.split(":")[1])
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
                    self.utc_minus_tai = timedelta(seconds=float(line.split(":")[1]))
                elif line.startswith("Experiment nominal start:"):
                    self.nominal_start = datetime.strptime(
                        line.split(":", 1)[1].strip()[:23], "%Y.%m.%d-%H:%M:%S.%f"
                    )
                elif line.startswith("Experiment nominal end:"):
                    self.nominal_end = datetime.strptime(
                        line.split(":", 1)[1].strip()[:23], "%Y.%m.%d-%H:%M:%S.%f"
                    )
                elif line.startswith("# Generated    at"):
                    cols = line.split()
                    if len(cols) == 4:
                        self.hostname = cols[3]
                    else:
                        self.hostname = ""
                elif line.startswith("#                 by PIMA"):
                    if "version" in line:
                        self.pima_version = line.split()[5]
                    else:
                        self.pima_version = line.split()[4]

        self.accum_length = (acc_min + acc_max) / 2.0
        self.no_auto_points_num = no_auto1 + no_auto2


@dataclass
class ActaFileHeader:
    """Header of the PIMA ACTA file."""

    experiment: str = ""
    station: str = ""
    source: str = ""
    scan_name: str = ""
    scan: int = 0
    obs: int = 0
    start_date: datetime = datetime.min
    stop_date: datetime = datetime.min
    duration: float = 0.0
    num_of_points: int = 0
    polar: str | None = None  # not in ACTA file, set externally

    def update(self, key: str, val: str) -> None:
        """Update header field according to key."""
        match key:
            case "Experiment":
                self.experiment = val
            case "Station":
                self.station = val
            case "Source":
                self.source = val
            case "Scan_name":
                self.scan_name = val
            case "Scan_index":
                self.scan = int(val)
            case "Observation_index":
                self.obs = int(val)
            case "Start_date":
                self.start_date = datetime.strptime(val[:23], "%Y.%m.%d-%H:%M:%S.%f")
            case "Stop_date":
                self.stop_date = datetime.strptime(val[:23], "%Y.%m.%d-%H:%M:%S.%f")
            case "Number_of_points":
                self.num_of_points = int(val)
            case "Effective duration":
                self.duration = float(val)
            case _:
                logger.debug("unused key in ACTA header: %s", key)


class ActaFile:
    """Represents PIMA ``ACTA`` file format and its contents."""

    supported_versions = (
        "# ACTA Output.  Format version of 2014.04.19",
        "# ACTA Output.  Format version of 2021.07.15",
    )

    def __init__(
        self,
        input_file_name: str,
        polar: str | None = None,
        utc_tai: timedelta | None = None,
    ):
        """
        Read and parse ACTA-file.

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
        self._header = ActaFileHeader(polar=polar)

        with open(input_file_name) as file:
            magic = file.readline().strip()
            if not magic.startswith(self.supported_versions):
                raise ValueError(f"{input_file_name} is not PIMA ACTA-file")

            for line in file:
                if line.startswith("#"):
                    if ":" in line:
                        key, val = line[1:].split(":", 1)
                        key = key.strip()
                        val = val.strip()
                        self._header.update(key, val)
                elif line.startswith("ACRL"):
                    cols = line.split()
                    self._if.append(int(cols[2]))
                    self._channel.append(int(cols[4]))
                    self._freq.append(1e-6 * float(cols[6]))
                    self._ampl.append(float(cols[9]))

        # Convert TAI to UTC if needed
        if utc_tai:
            if self._header.start_date:
                self._header.start_date += utc_tai
            if self._header.stop_date:
                self._header.stop_date += utc_tai

    @property
    def header(self) -> ActaFileHeader:
        """Return dictionary with header information."""
        return self._header

    @property
    def if_num(self) -> list[int]:
        """Return list of the IF numbers."""
        return self._if

    @property
    def channel(self) -> list[int]:
        """Return list of the channel numbers."""
        return self._channel

    @property
    def freq(self) -> list[float]:
        """Return list of the frequencies of the autospectrum."""
        return self._freq

    @property
    def ampl(self) -> list[float]:
        """Return list of the amplitudes of the autospectrum."""
        return self._ampl

    def plot(self, out_dir_base: str | os.PathLike, out_format: str = "pdf") -> None:
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

        date = self.header.start_date
        date_str = date.strftime("%Y-%m-%d %H:%M:%S")
        exper, band = self.header.experiment.split("_")

        ax.set_title(
            "{} - {}({}) - {} - {}".format(
                self.header.station,
                exper,
                band.upper(),
                self.header.polar,
                date_str,
            )
        )
        ax.set_xlabel("Frequency (MHz)")
        ax.set_ylabel("Amplitude")
        ax.grid(True)
        central_freq = np.mean(self.freq)
        logger.debug("plot_autospectra: central_freq = %s", central_freq)
        freq = np.asarray(self.freq) - central_freq
        ax.plot(freq, self.ampl, marker="o", ms=2)
        ax.set_xlim(-16, 16)

        date_str = date.strftime("%Y%m%dT%H%M")
        out_file = "AUTOSPEC_{}_{}_{}_{}.{}".format(
            date_str,
            self.header.experiment,
            self.header.polar,
            self.header.station,
            out_format,
        )

        out_dir = os.path.join(out_dir_base, self.header.experiment)
        os.makedirs(out_dir, exist_ok=True)
        out_file = os.path.join(out_dir, out_file)
        fig.savefig(out_file, format=out_format)


@dataclass
class TextTable1D:
    """Read PIMA 1D text table."""

    file_path: InitVar[str | os.PathLike]

    plot_title: str = ""
    subtitle: str = ""
    num_points: int = 0
    axis1_title: str = ""
    axis1_units: str = ""
    axis1_min: float = 0.0
    axis1_max: float = 0.0
    axis2_title: str = ""
    axis2_units: str = ""
    axis2_min: float = 0.0
    axis2_max: float = 0.0
    axis1_data: list[float] = field(default_factory=list)
    axis2_data: list[float] = field(default_factory=list)

    def __post_init__(self, file_path: str | os.PathLike):
        FORMAT_STRING = "# 1D text table.  Format version of 2012.12.30"

        with open(file_path) as file:
            magic = file.readline().strip()
            if magic != FORMAT_STRING:
                raise ValueError(f'Bad format string in file "{file_path}"')
            for line in file:
                if line.startswith("#"):
                    continue

                key, val = line.split(":", 1)
                key = key.strip()
                val = val.strip()
                match key:
                    case "PLOT_TITLE":
                        self.plot_title = val
                    case "SUBTITLE":
                        self.subtitle = val
                    case "NUM_POINTS":
                        self.num_points = int(val)
                    case "AXIS1_TITLE":
                        self.axis1_title = val
                    case "AXIS1_UNITS":
                        self.axis1_units = val
                    case "AXIS1_MIN":
                        self.axis1_min = float(val.replace("D", "e"))
                    case "AXIS1_MAX":
                        self.axis1_max = float(val.replace("D", "e"))
                    case "AXIS2_NAME":
                        self.axis2_title = val
                    case "AXIS2_UNITS":
                        self.axis2_units = val
                    case "AXIS2_MIN":
                        self.axis2_min = float(val.replace("D", "e"))
                    case "AXIS2_MAX":
                        self.axis2_max = float(val.replace("D", "e"))
                    case "POINT":
                        data = val.split()
                        self.axis1_data.append(float(data[1].replace("D", "e")))
                        self.axis2_data.append(float(data[2].replace("D", "e")))
                    case _:
                        logger.debug("unknown key in file %s: %s", file_path, key)


class Obs(NamedTuple):
    """PIMA observation information."""

    obs: int
    scan: int
    time_code: str
    source: str
    sta1: str
    sta2: str
    start_time: datetime
    stop_time: datetime
    ap_num: int


class ObsFriPlots(NamedTuple):
    """Observation fringe fitting plots."""

    obs: Obs
    frq_amp_file: str | None
    frq_amp_data: TextTable1D | None
    frq_phs_file: str | None
    frq_phs_data: TextTable1D | None
    tim_amp_file: str | None
    tim_amp_data: TextTable1D | None
    tim_phs_file: str | None
    tim_phs_data: TextTable1D | None


class ClockModelRec(NamedTuple):
    """Clock model record."""

    sta: str
    time: datetime
    clock_offset: float
    clock_rate: float
    group_delay: float
    delay_rate: float


class Pima:
    """Pima class is analog of pima_fringe.csh script."""

    def __init__(
        self, experiment_code: str, band: str, work_dir: str | os.PathLike | None = None
    ):
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

        # Dictionary with all parameters from PIMA control file
        self.cnt_params = {}
        self._update_cnt_params()
        self.exper_info = ExperInfo(self.exper, self.band)
        # Update exper_info if experiment already loaded.
        stt_file = os.path.join(
            self.cnt_params["EXPER_DIR:"], self.cnt_params["SESS_CODE:"] + ".stt"
        )
        if os.path.isfile(stt_file):
            self.exper_info.update(stt_file)

    def _update_cnt_params(self) -> None:
        """Read *PIMA* control file and fill cnt_params dictionary."""
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

    def update_cnt(self, opts: dict[str, Any]) -> None:
        """Update *PIMA* control file according to `opts` dictionary."""
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
                    line = f"# Last update on  {str(datetime.now())}\n"

                lines.append(line)

        with open(self.cnt_file_name, "w") as file:
            file.write("".join(lines))

        self._update_cnt_params()

    def _exec(
        self,
        operation: str,
        options: dict[str, Any] | None = None,
        log_name: str | None = None,
    ) -> int:
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

    def _error(self, msg: str) -> NoReturn:
        """Raise pima.Error exception."""
        self.logger.error(msg)
        raise Error(self.exper, self.band, msg)

    def load(self) -> None:
        """
        Run pima load.

        This function removes all existing index files of the current session and
        runs ``load`` task.

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

    def coarse(self, params: dict[str, Any] | None = None) -> str:
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

    def fine(self, params: dict[str, Any] | None = None) -> str:
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

    def bpas(self, params: dict[str, Any] | None = None) -> str:
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

    def split(self, tim_mseg: int = 1, params: dict[str, Any] | None = None):
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

    def load_gains(self, gain_file: str, params: dict[str, Any] | None = None) -> None:
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

    def load_tsys(self, tsys_file: str, params: dict[str, Any] | None = None) -> None:
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

    def acta(self, params: dict[str, Any] | None = None) -> list[str]:
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
    def set_polar(self, polar: str) -> None:
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

    def set_frq_grp(self, frq_grp: int) -> None:
        """
        Set frequency group in control file.

        Parameters
        ----------
        frq_grp : int
            Frequency group number.

        """
        if frq_grp < 1 or frq_grp > self.exper_info.frq_grp:
            self._error(f"Invalid frequency group {frq_grp}")

        self.logger.info("Set frequency group to %s", frq_grp)
        self.update_cnt({"FRQ_GRP:": frq_grp})

    @property
    def ap_minmax(self) -> tuple[float, float]:
        """Return minimum and maximum accummulation periods in experiment."""
        ap_min = ap_max = 0.0
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
    def number_of_deselected_points(self) -> int:
        """Return total number of deselected points."""
        return self.exper_info.deselected_points_num

    def station_list(self, ivs_name=True) -> list[str]:
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
    def source_list(self) -> list[tuple[str, str, str]]:
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

    def source_dist(self) -> dict[str, float]:
        """Return distance between correlator phase center and source position."""
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
    def obs_number(self) -> int:
        """Return number of observations in the experiment."""
        return self.exper_info.obs_num

    @property
    def chan_number(self) -> int:
        """Return number of spectral channels in uv-data."""
        return self.exper_info.sp_chann_num

    @property
    def frequencies(self) -> list:
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
    def observations(self) -> list[Obs]:
        """Return list of observations with some information."""
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

                    obs_list.append(
                        Obs(
                            obs=int(toks[0]),
                            scan=int(toks[1]),
                            time_code=toks[2],
                            source=toks[8],
                            sta1=toks[9],
                            sta2=toks[10],
                            start_time=datetime.strptime(
                                toks[4], "%Y.%m.%d-%H:%M:%S.%f,"
                            ),
                            stop_time=datetime.strptime(
                                toks[5], "%Y.%m.%d-%H:%M:%S.%f"
                            ),
                            ap_num=ap_num,
                        )
                    )

        return obs_list

    def clock_model(self) -> list[ClockModelRec]:
        """
        Read clock model components from the **PIMA** ``mdc`` file.

        Returns
        -------
        clock_model : list
            List of tuples with clock model components.

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
                        ClockModelRec(
                            sta, time, clock_offset, clock_rate, group_delay, delay_rate
                        )
                    )

        return clock_model

    def mk_exclude_obs_file(
        self, obs_list: list[int], suffix: str, polar: str | None = None
    ) -> str:
        """
        Create ``EXCLUDE_OBS_FILE`` file using list of the observation indices.

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

        Notes
        -----
            If `obs_list` is empty or ``None`` function removes existing file.

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

    def set_mask_file(self, mask_gen_file: str | None) -> str | None:
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

    def mk_bpass_mask_gen(self, params: list[tuple]) -> str | None:
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
                    "#  Control for bandpass mask generation for experiment {}".format(
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

    def fringe_plots(self) -> list[ObsFriPlots]:
        """Return list of plot file generated by ``frib`` task."""
        extentions = {
            "XW": "",
            "SAV": "sav",
            "GIF": "gif",
            "PS": "ps",
            "TXT": "txt",
            "NO": "",
        }
        frq_extention = extentions[self.cnt_params["FRIB.1D_RESFRQ_PLOT:"]]
        tim_extention = extentions[self.cnt_params["FRIB.1D_RESTIM_PLOT:"]]
        fpl_dir = os.path.join(
            self.cnt_params["EXPER_DIR:"], f"{self.cnt_params['SESS_CODE:']}_fpl"
        )
        if not os.path.isdir(fpl_dir):
            logger.debug("fringe_plots: %s does not exist", fpl_dir)
            return []

        obs_list = []
        for obs in self.observations:
            sta1 = f"{obs.sta1.lower():_<8}"
            sta2 = f"{obs.sta2.lower():_<8}"
            time_code = f"{obs.time_code:_<9}"

            frq_amp_file = None
            frq_phs_file = None
            tim_amp_file = None
            tim_phs_file = None

            frq_amp_data = None
            frq_phs_data = None
            tim_amp_data = None
            tim_phs_data = None

            if frq_extention:
                frq_amp_file = os.path.join(
                    fpl_dir,
                    "fr1d_frq_{}_{}_{}_{}_amp.{}".format(
                        time_code, self.band, sta1, sta2, frq_extention
                    ),
                )
                self.logger.debug("fringe_plots: frq_amp_file: %s", frq_amp_file)
                if not os.path.isfile(frq_amp_file):
                    frq_amp_file = None
                elif frq_extention == "txt":
                    frq_amp_data = TextTable1D(frq_amp_file)

                frq_phs_file = os.path.join(
                    fpl_dir,
                    "fr1d_frq_{}_{}_{}_{}_phs.{}".format(
                        time_code, self.band, sta1, sta2, frq_extention
                    ),
                )
                self.logger.debug("fringe_plots: frq_phs_file: %s", frq_phs_file)
                if not os.path.isfile(frq_phs_file):
                    frq_phs_file = None
                elif frq_extention == "txt":
                    frq_phs_data = TextTable1D(frq_phs_file)

            if tim_extention:
                tim_amp_file = os.path.join(
                    fpl_dir,
                    "fr1d_tim_{}_{}_{}_{}_amp.{}".format(
                        time_code, self.band, sta1, sta2, tim_extention
                    ),
                )
                if not os.path.isfile(tim_amp_file):
                    tim_amp_file = None
                elif tim_extention == "txt":
                    tim_amp_data = TextTable1D(tim_amp_file)

                tim_phs_file = os.path.join(
                    fpl_dir,
                    "fr1d_tim_{}_{}_{}_{}_phs.{}".format(
                        time_code, self.band, sta1, sta2, tim_extention
                    ),
                )
                if not os.path.isfile(tim_phs_file):
                    tim_phs_file = None
                elif tim_extention == "txt":
                    tim_phs_data = TextTable1D(tim_phs_file)

            obs_list.append(
                ObsFriPlots(
                    obs,
                    frq_amp_file,
                    frq_amp_data,
                    frq_phs_file,
                    frq_phs_data,
                    tim_amp_file,
                    tim_amp_data,
                    tim_phs_file,
                    tim_phs_data,
                )
            )

        return obs_list


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

    txt_file = str(fits_file_path.with_suffix(".txt"))

    cmd_line = ["fits_to_radplot", "-o", txt_file, fits_file]

    return subprocess.check_output(cmd_line, text=True)


def acta_plot(input_file: str, output_file: str) -> str:
    """Run ``acta_plot``."""
    if not os.path.isfile(input_file):
        raise FileNotFoundError(f"file {input_file} does not exist")

    cmd_line = ["acta_plot", input_file, output_file]
    out = subprocess.check_output(cmd_line, text=True)

    return out


def bpas_log_snr_new(file_name: str, mode: str = "INIT") -> dict[int, float]:
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
