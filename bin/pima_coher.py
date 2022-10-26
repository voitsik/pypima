#!/usr/bin/env python3
"""
Created on Fri Apr  4 18:22:19 2014

@author: Petr Voytsik
"""

import argparse
import logging
import math
import os.path
import sys
from datetime import datetime
from tempfile import NamedTemporaryFile

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from pypima import Fri, Pima, PimaError


class CoherAnalyzer:
    def __init__(self, exper, band):
        self.sta1 = ""
        self.sta2 = ""
        self.obs = 0
        self.start_time = None
        # self.dur_arr = None
        # self.ampl_arr = None
        # self.snr_arr = None
        self.data = {}

        self.pim = Pima(exper, band)

    def _mk_txt_file_name(self):
        """Return file name."""
        return "{}_{}_{:02d}_coher.{}".format(
            self.pim.exper, self.pim.band, self.obs, "txt"
        )

    def analyze_obs(self, obs, max_dur, accel=True, force=False):
        """
        Analyze coherent loses in observation `obs`.

        Parameters
        ----------
        obs : int
            Observation id.
        max_dur : float
            Maximum scal length.
        force : bool
            Force rerun ``pima`` even if the output txt file exists.

        """
        self.obs = obs
        # self.dur_arr = None
        # self.ampl_arr = None
        # self.snr_arr = None
        self.data.clear()

        if self.obs <= 0 or self.obs > self.pim.obs_number:
            logging.error(
                "Incorrect observation number %s, must be in range [%s %s]",
                self.obs,
                1,
                self.pim.obs_number,
            )
            raise ValueError

        # Read info from obs-file
        obs_info = self.pim.observations[self.obs - 1]
        self.start_time = obs_info.start_time + self.pim.exper_info["utc_minus_tai"]
        self.sta1 = obs_info.sta1
        self.sta2 = obs_info.sta2

        txt_file_name = self._mk_txt_file_name()
        if os.path.isfile(txt_file_name):
            logging.info("file %s already exists", txt_file_name)
            if force:
                os.remove(txt_file_name)
                logging.info("remove it!")
            else:
                values = np.loadtxt(txt_file_name, unpack=True)
                if values.shape[0] == 4:
                    (
                        self.data["dur"],
                        self.data["skip"],
                        self.data["ampl"],
                        self.data["snr"],
                    ) = values
                elif values.shape[0] == 3:
                    (
                        self.data["dur"],
                        self.data["ampl"],
                        self.data["snr"],
                    ) = values
                else:
                    logging.error("Invalid input file %s", txt_file_name)

                return

        # Prepare temporary fri and frr files
        with NamedTemporaryFile(suffix=".fri", delete=False) as tmp:
            tmp_fri = tmp.name
        with NamedTemporaryFile(suffix=".frr", delete=False) as tmp:
            tmp_frr = tmp.name

        params = {
            "FRIB.OBS:": str(self.obs),
            "SCAN_LEN_SKIP:": "0.0",
            "SCAN_LEN_USED:": str(max_dur),
            "FRINGE_FILE:": tmp_fri,
            "FRIRES_FILE:": tmp_frr,
        }

        if not accel:
            params["PHASE_ACCELERATION:"] = "0"
            params["PHASE_ACCEL_MIN:"] = "0.0"
            params["PHASE_ACCEL_MAX:"] = "0.0"
            params["FRIB.FINE_SEARCH:"] = "LSQ"

        fri = Fri(self.pim.fine(params=params))

        if not fri:
            logging.error("PIMA fri-file is empty")
            return

        if fri[0]["SNR"] < 7.0:
            logging.warning("SNR of obs #%s is too low, skip it.", obs)
            return

        # Update from fri-file
        self.start_time = fri[0]["start_time"] + self.pim.exper_info["utc_minus_tai"]
        self.sta1 = fri[0]["sta1"]
        self.sta2 = fri[0]["sta2"]

        dur_arr = []
        skip_arr = []
        ampl_arr = []
        snr_arr = []
        full_duration = fri[0]["duration"]

        for divisor in [20, 15, 12, 10, 8, 6, 5, 4, 3, 2, 1.7, 1.5, 1.3, 1.1, 1]:
            dur = round(full_duration / divisor)

            #        if dur < 60:
            #            continue

            logging.debug("Set SCAN_LEN_USED %s", dur)

            for j in range(math.ceil(divisor)):
                if divisor < 2 and j == 1:
                    skip = full_duration - dur
                else:
                    skip = j * dur

                logging.debug("Set SCAN_LEN_SKIP to %s", skip)
                params["SCAN_LEN_SKIP:"] = str(skip)
                params["SCAN_LEN_USED:"] = str(dur)
                fri = Fri(self.pim.fine(params=params))

                if not fri:
                    logging.warning("Skip empty fri-file")
                    continue

                logging.debug("SNR = %s", fri[0]["SNR"])

                if fri[0]["SNR"] < 5.7:
                    continue

                dur_arr.append(fri[0]["duration"])
                skip_arr.append(skip)
                ampl_arr.append(fri[0]["ampl_lsq"])
                snr_arr.append(fri[0]["SNR"])

                print(f"{dur_arr[-1]} {ampl_arr[-1]} {snr_arr[-1]}")

        self.data["dur"] = np.asarray(dur_arr)
        self.data["skip"] = np.asanyarray(skip_arr)
        self.data["ampl"] = np.asarray(ampl_arr)
        self.data["snr"] = np.asarray(snr_arr)

        # Delete temporary files
        if os.path.isfile(tmp_fri):
            os.remove(tmp_fri)
        if os.path.isfile(tmp_frr):
            os.remove(tmp_frr)

    def save_txt(self):
        """Save data as a text table."""
        # assert isinstance(self.dur_arr, np.ndarray)
        # assert isinstance(self.ampl_arr, np.ndarray)
        # assert isinstance(self.snr_arr, np.ndarray)

        file_name = self._mk_txt_file_name()
        if os.path.isfile(file_name):
            logging.warning("file %s already exists", file_name)
            return

        # data = np.vstack((self.dur_arr, self.ampl_arr, self.snr_arr)).T
        data = np.vstack(
            (self.data["dur"], self.data["skip"], self.data["ampl"], self.data["snr"])
        ).T

        header = "{:^4} {:^6} {:^10} {:^8}".format("solint", "skip", "ampl", "SNR")

        np.savetxt(
            file_name, data, fmt=["%6.1f", "%6d", "%10.3e", "%8.1f"], header=header
        )

    def _plot_theo_curves(self, ax_ampl, ax_snr, max_coher_dur=200):
        """Plot theoretical curves for amplitude and SNR."""
        # Mask to select values with int time less then coher time
        coher_mask = self.data["dur"] < max_coher_dur

        # At least 2 points to fit
        if np.count_nonzero(coher_mask) <= 1:
            return

        def sqrt_func(time, ampl):
            """Return A * sqrt(T)."""
            return ampl * np.sqrt(time)

        fit_params, _ = curve_fit(
            sqrt_func, self.data["dur"][coher_mask], self.data["snr"][coher_mask]
        )

        dur_theo_arr = np.linspace(0, self.data["dur"].max() * 1.01, num=300)

        if ax_ampl is not None:
            ampl_theo_arr = np.full_like(
                dur_theo_arr, self.data["ampl"][coher_mask].mean()
            )
            ax_ampl.plot(dur_theo_arr, ampl_theo_arr, "r-")

        if ax_snr is not None:
            snr_theo_arr = sqrt_func(dur_theo_arr, *fit_params)
            ax_snr.plot(dur_theo_arr, snr_theo_arr, "r-")

    def plot(self, out_format="pdf"):
        """Plot curves."""
        # assert isinstance(self.dur_arr, np.ndarray)
        # assert isinstance(self.ampl_arr, np.ndarray)
        # assert isinstance(self.snr_arr, np.ndarray)

        if self.data["dur"].size == 0:
            logging.warning("Nothing to plot...")
            return

        # Plotting
        fig, (ax1, ax2) = plt.subplots(2, sharex=True)

        ax1.yaxis.get_major_formatter().set_powerlimits((-4, 6))

        ax1.plot(self.data["dur"], self.data["ampl"], ".")
        ax1.plot(0, 0, alpha=0)  # Dirty hack
        ax2.plot(self.data["dur"], self.data["snr"], ".")

        self._plot_theo_curves(ax1, ax2)

        ax1.set_ylabel("Amplitude")
        ax1.grid(True)

        ax2.set_xlabel("Integration time (s)")
        ax2.set_ylabel("SNR")
        ax2.grid(True)

        ax2.text(
            0.95,
            0.01,
            f"Generated at {datetime.now():%Y-%m-%d %H:%M}",
            color="gray",
            size="xx-small",
            ha="right",
            transform=ax2.transAxes,
        )

        title = "{:8}({:1}) {:8} / {:8} [{:%Y-%m-%d %H:%M} UTC]".format(
            self.pim.exper.lower(),
            self.pim.band.upper(),
            self.sta1,
            self.sta2,
            self.start_time,
        )
        fig.suptitle(title)

        plot_file_name = "{}_{}_{:02d}_coher.{}".format(
            self.pim.exper, self.pim.band, self.obs, out_format
        )
        if out_format == "png":
            dpi = 300
        else:
            dpi = None

        fig.savefig(plot_file_name, bbox_inches="tight", pad_inches=0.1, dpi=dpi)

    def plot_snr(self, out_format="pdf"):
        """Plot SNR vs integration time."""
        band2freq = {"l": 1.668, "c": 4.836, "k": 22.236}
        ivs2sta = {"RADIO-AS": "SRT", "GBT-VLBA": "GBT", "EFLSBERG": "Effelsberg"}

        if self.data["dur"].size == 0:
            logging.warning("Nothing to plot...")
            return

        fig, ax = plt.subplots(figsize=(5, 3))

        ax.plot(self.data["dur"], self.data["snr"], ".")

        self._plot_theo_curves(None, ax)

        ax.set_xlabel("Integration time (s)")
        ax.set_ylabel("SNR")
        ax.grid(True, linestyle=":")

        try:
            sta1_txt = ivs2sta[self.sta1]
        except KeyError:
            sta1_txt = self.sta1.capitalize()

        try:
            sta2_txt = ivs2sta[self.sta2]
        except KeyError:
            sta2_txt = self.sta2.capitalize()

        ax.text(
            0.02,
            0.79,
            "{:%Y-%m-%d %H:%M} UT\n{}-{}\n{:.1f} GHz".format(
                self.start_time, sta1_txt, sta2_txt, band2freq[self.pim.band]
            ),
            bbox=dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9),
            transform=ax.transAxes,
        )

        plot_file_name = "{}_{}_{:02d}_coher_snr.{}".format(
            self.pim.exper, self.pim.band, self.obs, out_format
        )
        if out_format == "png":
            dpi = 300
        else:
            dpi = None

        fig.savefig(plot_file_name, bbox_inches="tight", pad_inches=0.1, dpi=dpi)


def valid_obs_list(string):
    """Convert a string with comma-separated values to list of integers."""
    try:
        return [int(obs) for obs in string.split(",")]
    except ValueError as exc:
        msg = f"Not a valid comma-separated list of integers: '{string}'."
        raise argparse.ArgumentTypeError(msg) from exc


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("exper", help="experiment code")
    parser.add_argument("band", help="frequency band")
    parser.add_argument(
        "obs_list",
        type=valid_obs_list,
        help="comma-separated list of observation numbers",
    )

    parser.add_argument(
        "--scan-length", type=float, default=1500.0, help="full scan length"
    )
    parser.add_argument(
        "--no-accel", action="store_true", help="disable acceleration fit"
    )
    parser.add_argument(
        "--plot-format", default="pdf", help="output plot file format (PDF by default)"
    )
    parser.add_argument(
        "--for-paper", action="store_true", help="plot in format suitable for paper"
    )
    parser.add_argument("--force", action="store_true", help="delete existing txt file")
    parser.add_argument("-v", "--verbose", action="store_true", help="be more verbose")

    return parser.parse_args()


def main():
    """Main."""
    args = parse_args()

    logging.basicConfig(
        format="%(asctime)s %(levelname)s: %(message)s", level=logging.INFO
    )

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    if args.plot_format not in plt.gcf().canvas.get_supported_filetypes().keys():
        logging.error("Format %s does not supported by matplotlib", args.plot_format)

        return 1

    try:
        analyzer = CoherAnalyzer(args.exper, args.band)
    except PimaError as err:
        logging.error("PIMA Error: %s", err)
        return 1
    except OSError as err:
        logging.error("OSError: %s", err)
        return 1

    for obs_id in args.obs_list:
        try:
            analyzer.analyze_obs(
                obs_id, args.scan_length, accel=not args.no_accel, force=args.force
            )
        except PimaError as err:
            logging.error("PIMA Error: %s", err)
            return 1

        analyzer.save_txt()

        if args.for_paper:
            analyzer.plot_snr(args.plot_format)
        else:
            analyzer.plot(args.plot_format)

    return 0


if __name__ == "__main__":
    sys.exit(main())
