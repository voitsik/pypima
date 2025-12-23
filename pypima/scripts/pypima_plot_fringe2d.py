#!/usr/bin/env python3
"""Plot 2D image of the fringe on delay/rate plane."""

import argparse
import logging
import sys
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from pypima.pima import Error as PimaError
from pypima.pima import Pima, TextTable2D

logger = logging.getLogger(Path(__file__).stem)

FORMAT = "# 2D text table.  Format version of 2012.12.30"


def plot(data: TextTable2D, title=None, save=False, plot_format="pdf"):
    """Plot PIMA 2D file."""
    mean = np.mean(data.axis3_data)
    std = np.std(data.axis3_data)
    mean_corr = np.mean(data.axis3_data[data.axis3_data < mean + 2 * std])

    print(f"Mean = {mean:.3e}")
    print(f"Mean_corr = {mean_corr:.3e}")
    print(f"Peak = {data.axis3_data.max():.3e}")
    print(f"std = {std:.3e}")
    print(f"Peak/mean = {data.axis3_data.max() / mean_corr:.1f}")
    print(f"Peak/std = {data.axis3_data.max() / std:.1f}")
    # print ("Peak/mean1 = " + str(Z.max() / mean1))

    fig, ax = plt.subplots(layout="constrained", figsize=(8, 6), dpi=150)

    if title is not None:
        ax.set_title(title, size="large")
    else:
        fig.suptitle(data.plot_title)

    ax.set_ylabel(r"Group delay ($\mu s$)", size="large")
    ax.set_xlabel(r"Fringe rate ($\mu s/s$)", size="large")

    extent = (
        1e6 * data.axis2_min,
        1e6 * data.axis2_max,
        1e6 * data.axis1_min,
        1e6 * data.axis1_max,
    )
    dx = (extent[1] - extent[0]) / float(data.axis2_num_points)
    dy = (extent[3] - extent[2]) / float(data.axis1_num_points)
    aspect = dx / dy
    logger.debug("dx = %.3e us", dx)
    logger.debug("dy = %.3e us/s", dy)
    logger.debug("aspect = %.3e", aspect)

    im = ax.imshow(
        data.axis3_data / mean_corr,
        extent=extent,
        origin="lower",
        cmap="viridis",
        interpolation="none",
        aspect=aspect,  # "equal",
        # vmin=0.0,
        # vmax=data.axis3_data.max(),
    )

    cbar = fig.colorbar(im)
    cbar.ax.set_ylabel("Signal-to-Noise Ratio", size="large")

    if save:
        file_name = f"Fringe_plot2d.{plot_format}"

        if plot_format == "png":
            dpi = 150
        else:
            dpi = None

        plt.savefig(file_name, bbox_inches="tight", dpi=dpi, pad_inches=0.1)
    else:
        plt.show()


def run_frib(
    pim: Pima,
    obs: int,
    delay_window: float | None = None,
    rate_window: float | None = None,
) -> None:
    """Run pima frib task for given observation."""
    if delay_window is not None:
        delay_window = delay_window * 1e-6
    else:
        delay_window = 1e-6  # 1 microsecond

    if rate_window is not None:
        rate_window = rate_window
    elif pim.band == "k":
        rate_window = 2e-12
    elif pim.band == "c":
        rate_window = 4e-12
    else:
        rate_window = 8e-12

    params = {
        "FRIB.OBS:": str(obs),
        "FRIB.2D_FRINGE_PLOT:": "TXT",
        "FRIB.PLOT_DELAY_WINDOW_WIDTH:": f"{delay_window:E}",
        "FRIB.PLOT_RATE_WINDOW_WIDTH:": f"{rate_window:E}",
    }

    with tempfile.NamedTemporaryFile(suffix=".fri") as tmp_fri:
        params.update({"FRINGE_FILE:": tmp_fri.name})
        pim.fine(params)


def get_data2d(
    pim: Pima,
    obs_num: int,
    delay_window: float | None = None,
    rate_window: float | None = None,
) -> TextTable2D:
    """Return 2D fringe data."""
    obs = pim.observations[obs_num - 1]

    path = (
        Path(pim.cnt_params["EXPER_DIR:"])
        / f"{pim.cnt_params['SESS_CODE:']}_fpl"
        / f"fr2d_{obs.time_code}__{pim.band}_{obs.sta1.lower():_<8}_{obs.sta2.lower()}.txt"
    )

    logger.debug("2D fringe file path: %s", path)

    if not path.is_file():
        logger.info(
            "2D fringe file not found (%s), running FRIB to generate it...", path
        )
        run_frib(pim, obs_num, delay_window, rate_window)

    return TextTable2D(path)


def setup_logging(args):
    """Set up logging."""
    root = logging.getLogger()
    root.setLevel(logging.WARNING)  # Use WARNING by default

    if args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    handler = logging.StreamHandler()
    formatter = logging.Formatter("%(levelname)s: %(message)s")
    handler.setFormatter(formatter)
    root.addHandler(handler)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("exper", type=str.lower, help="experiment code")
    parser.add_argument("band", type=str.lower, help="frequency band code")
    parser.add_argument("obs", type=int, help="observation number")

    parser.add_argument(
        "--delay-window", type=float, help="delay window in microseconds"
    )
    parser.add_argument(
        "--rate-window", type=float, help="rate window in seconds per second"
    )

    parser.add_argument("-t", "--title", help="custom plot title")
    parser.add_argument("-s", "--save", action="store_true", help="save plot to file")
    parser.add_argument(
        "-f",
        "--plot-format",
        default="pdf",
        help="format of the output file, default is PDF",
    )

    parser.add_argument(
        "-d", "--debug", action="store_true", help="enable debug logging"
    )

    return parser.parse_args()


def main() -> int:
    """Run main."""
    args = parse_args()
    setup_logging(args)

    pim = Pima(args.exper, args.band)

    if pim.obs_number <= 0:
        logger.error(
            "No observations found for experiment %s(%s)", args.exper, args.band
        )
        logger.error("Please run 'pima load' first to load data")
        return 1

    if args.obs <= 0 or args.obs > pim.obs_number:
        logger.error(
            "Incorrect observation number %s must be in range [ 1 %s ]",
            args.obs,
            pim.obs_number,
        )
        return 1

    try:
        data = get_data2d(pim, args.obs, args.delay_window, args.rate_window)
    except PimaError as err:
        logger.error("PIMA Error: %s", err)
        return 1
    except OSError as err:
        logger.error("OSError: %s", err)
        return 1

    plot(data, args.title, args.save, args.plot_format)

    return 0


if __name__ == "__main__":
    sys.exit(main())
