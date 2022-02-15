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
        self.dur_arr = None
        self.ampl_arr = None
        self.snr_arr = None

        self.pim = Pima(exper, band)

    def analyze_obs(self, obs, max_dur, accel=True):
        """
        Analyze coherent loses in observation `obs`.

        Parameters
        ----------
        obs : int
            Observation id.
        max_dur : float
            Maximum scal length.

        """
        self.obs = obs
        self.dur_arr = None
        self.ampl_arr = None
        self.snr_arr = None

        if self.obs <= 0 or self.obs > self.pim.obs_number:
            logging.error(
                "Incorrect observation number %s, must be in range [%s %s]",
                self.obs,
                1,
                self.pim.obs_number,
            )
            raise ValueError

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

        full_duration = fri[0]["duration"]
        self.start_time = fri[0]["start_time"] + self.pim.exper_info["utc_minus_tai"]
        self.sta1 = fri[0]["sta1"]
        self.sta2 = fri[0]["sta2"]

        dur_arr = []
        ampl_arr = []
        snr_arr = []

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
                ampl_arr.append(fri[0]["ampl_lsq"])
                snr_arr.append(fri[0]["SNR"])

                print("{} {} {}".format(dur_arr[-1], ampl_arr[-1], snr_arr[-1]))

        self.dur_arr = np.asarray(dur_arr)
        self.ampl_arr = np.asarray(ampl_arr)
        self.snr_arr = np.asarray(snr_arr)

        # Delete temporary files
        if os.path.isfile(tmp_fri):
            os.remove(tmp_fri)
        if os.path.isfile(tmp_frr):
            os.remove(tmp_frr)

    def save_txt(self):
        """Save data as a text table."""
        assert isinstance(self.dur_arr, np.ndarray)
        assert isinstance(self.ampl_arr, np.ndarray)
        assert isinstance(self.snr_arr, np.ndarray)

        data = np.vstack((self.dur_arr, self.ampl_arr, self.snr_arr)).T

        file_name = "{}_{}_{:02d}_coher.{}".format(
            self.pim.exper, self.pim.band, self.obs, "txt"
        )

        header = "{:^4} {:^10} {:^8}".format("solint", "ampl", "SNR")

        np.savetxt(file_name, data, fmt=["%6.1f", "%10.3e", "%8.1f"], header=header)

    def _plot_theo_curves(self, ax_ampl, ax_snr, max_coher_dur=400):
        """Plot theoretical curves for amplitude and SNR."""
        # Mask to select values with int time less then coher time
        coher_mask = self.dur_arr < max_coher_dur

        # At least 2 points to fit
        if np.count_nonzero(coher_mask) <= 1:
            return

        def sqrt_func(time, ampl):
            """Return A * sqrt(T)."""
            return ampl * np.sqrt(time)

        fit_params, _ = curve_fit(
            sqrt_func, self.dur_arr[coher_mask], self.snr_arr[coher_mask]
        )

        dur_theo_arr = np.linspace(0, self.dur_arr.max(), num=200)
        snr_theo_arr = sqrt_func(dur_theo_arr, *fit_params)
        ampl_theo_arr = np.full_like(dur_theo_arr, self.ampl_arr[coher_mask].mean())

        ax_ampl.plot(dur_theo_arr, ampl_theo_arr, "r-")
        ax_snr.plot(dur_theo_arr, snr_theo_arr, "r-")

    def plot(self, out_format="pdf"):
        """Plot curves."""
        assert isinstance(self.dur_arr, np.ndarray)
        assert isinstance(self.ampl_arr, np.ndarray)
        assert isinstance(self.snr_arr, np.ndarray)

        if self.dur_arr.size == 0 or self.ampl_arr.size == 0 or self.snr_arr.size == 0:
            logging.warning("Nothing to plot...")
            return

        # Plotting
        fig, (ax1, ax2) = plt.subplots(2, sharex=True)

        ax1.yaxis.get_major_formatter().set_powerlimits((-4, 6))

        ax1.plot(self.dur_arr, self.ampl_arr, ".")
        ax1.plot(0, 0, alpha=0)  # Dirty hack
        ax2.plot(self.dur_arr, self.snr_arr, ".")

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

        plt.savefig(plot_file_name, bbox_inches="tight", pad_inches=0.1, dpi=dpi)


def valid_obs_list(string):
    """Convert a string with comma-separated values to list of integers."""
    try:
        return [int(obs) for obs in string.split(",")]
    except ValueError:
        msg = f"Not a valid comma-separated list of integers: '{string}'."
        raise argparse.ArgumentTypeError(msg)


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
            analyzer.analyze_obs(obs_id, args.scan_length, accel=not args.no_accel)
        except PimaError as err:
            logging.error("PIMA Error: %s", err)
            return 1

        analyzer.save_txt()
        analyzer.plot(args.plot_format)

    return 0


if __name__ == "__main__":
    sys.exit(main())
