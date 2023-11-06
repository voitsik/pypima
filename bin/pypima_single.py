#!/usr/bin/env python3
"""
Created on Fri Dec 13 17:50:20 2013

@author: Petr Voytsik
"""

import argparse
import logging
import os.path
import sys
import tempfile

import psycopg2

from pypima import DB, PimaError, RaExperiment, RaExperimentError


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("exper", help="experiment code")
    parser.add_argument("band", help="frequency band")

    # Optional arguments
    parser.add_argument("--polar", help="polarization")
    parser.add_argument(
        "--gvlbi", "-g", action="store_true", help="process GVLBI FITS-file"
    )
    parser.add_argument(
        "--no-accel", action="store_true", help="disable acceleration term fitting"
    )
    parser.add_argument(
        "--force-small",
        action="store_true",
        help="force to use 64-channel FITS file (if any)",
    )
    parser.add_argument(
        "--autospec-only",
        action="store_true",
        help="generate autocorrelation spectra only",
    )
    parser.add_argument(
        "--individual-ifs",
        action="store_true",
        help="do fringe fittig for individual IFs",
    )
    parser.add_argument("--split", action="store_true", help="do SPLIT")
    parser.add_argument("--fits", nargs="+", help="external FITS-IDI file(s)")
    parser.add_argument("--orbit", help="reconstructed orbit file")
    parser.add_argument("--antab", help="ANTAB file")
    parser.add_argument("--scan-length", type=float, help="set scan length in seconds")
    parser.add_argument(
        "--beg-frq", type=int, help="start intermediate frequency (IF) index"
    )
    parser.add_argument(
        "--end-frq", type=int, help="end intermediate frequency (IF) index"
    )
    parser.add_argument(
        "--frequency-group", type=int, default=1, help="frequency group"
    )
    parser.add_argument(
        "--debug", "-d", action="store_true", help="enable debug output"
    )

    bpas_group = parser.add_argument_group("bandpass settings")
    bpas_group.add_argument("--ref-sta", metavar="STA", help="reference station")
    bpas_group.add_argument(
        "--no-bandpass", action="store_true", help="disable bandpass calibration"
    )
    bpas_group.add_argument(
        "--bpas-mode",
        metavar="MODE",
        choices=["INIT", "ACCUM", "FINE"],
        help="set bandpass calibration mode",
    )
    bpas_group.add_argument(
        "--bpas-use",
        metavar="BANDPASS_USE",
        default="PHS",
        choices=["AMP", "PHS", "AMP_PHS", "NO"],
        help="set BANDPASS_USE PIMA parameter",
    )
    bpas_group.add_argument(
        "--no-ampl-bpas",
        action="store_true",
        help="disable amplitude bandpass calibration",
    )
    bpas_group.add_argument(
        "--bpas-var",
        type=int,
        choices=[0, 1, 2, 3, 4, 5],
        default=3,
        help="predefined bandpass parameters",
    )
    bpas_group.add_argument(
        "--flag-chann",
        type=int,
        default=2,
        metavar="N",
        help="flag N edge spectral channels of the bandpass",
    )
    bpas_group.add_argument(
        "--bpas-norm",
        choices=["NO", "IF", "BAND"],
        default="IF",
        help="the way how the bandpass normalization is made",
    )
    bpas_group.add_argument(
        "--no-bpas-renorm", action="store_true", help="disable bandpass renormalization"
    )

    ff_group = parser.add_argument_group("fringe fitting settings")
    ff_group.add_argument(
        "--delay-window-width",
        metavar="DELAY_WIDTH",
        type=float,
        default=64.0,
        help="delay window width in microseconds",
    )
    ff_group.add_argument(
        "--rate-window-width",
        metavar="RATE_WIDTH",
        type=float,
        default=1e-8,
        help="rate window width in s/s",
    )

    return parser.parse_args()


def main():
    """Main"""
    args = parse_args()

    log_format = "%(asctime)s %(levelname)s: %(name)s: %(message)s"
    logging.basicConfig(format=log_format, level=logging.INFO)

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)

    exper = args.exper.lower()
    band = args.band.lower()
    polar = args.polar

    data_dir = os.getenv(
        "PYPIMA_DATA_DIR", default=os.path.join(os.getenv("HOME"), "data", "pima_data")
    )

    #
    # Check command line arguments
    #
    if args.scan_length:
        if args.scan_length <= 0:
            logging.error("scan_length must be positive")
            return 1

        scan_length = args.scan_length
    else:
        scan_length = 1500

    if args.beg_frq is not None and args.beg_frq <= 0:
        logging.error("beg_frq must be positive")
        return 1

    if args.end_frq is not None and args.end_frq <= 0:
        logging.error("beg_frq must be positive")
        return 1

    if args.orbit:
        if os.path.isfile(args.orbit):
            orbit = os.path.abspath(args.orbit)
        else:
            logging.error("orbit file %s does not exist", args.orbit)
            return 1
    else:
        orbit = None

    if args.antab:
        if os.path.isfile(args.antab):
            antab_file = os.path.abspath(args.antab)
        else:
            logging.error("ANTAB file %s does not exist", args.antab)
            return 1
    else:
        antab_file = None

    database = DB()

    try:
        ra_exp = RaExperiment(
            exper,
            band,
            database,
            gvlbi=args.gvlbi,
            data_dir=data_dir,
            uv_fits=args.fits,
            orbit=orbit,
        )
        ra_exp.init_workdir()

        ra_exp.load(
            update_db=False,
            force_small=args.force_small,
            scan_length=scan_length,
            beg_frq=args.beg_frq,
            end_frq=args.end_frq,
        )
        ra_exp.flag_edge_chann(args.flag_chann)

        if not polar:
            if band == "l":
                polar = "RR"
            else:
                polar = "LL"

        ra_exp.pima.set_polar(polar)
        ra_exp.pima.set_frq_grp(args.frequency_group)

        if args.autospec_only:
            # Define and create directory for auto spectrum plot files
            spec_out_dir = os.getenv(
                "PYPIMA_AUTOSPEC_DIR",
                default=os.path.join(os.getenv("HOME"), "pima_autospec"),
            )
            if not os.path.exists(spec_out_dir):
                os.mkdir(spec_out_dir)

            ra_exp.generate_autospectra(plot=True, out_dir=spec_out_dir, db=True)
        else:
            ra_exp.load_antab(antab_file)

            # Set fringe fitting window
            ra_exp.pima.update_cnt(
                {
                    "FRIB.DELAY_WINDOW_WIDTH:": f"{args.delay_window_width * 1e-6:e}",
                    "FRIB.RATE_WINDOW_WIDTH:": f"{args.rate_window_width:e}",
                }
            )

            if args.individual_ifs:
                if_num = ra_exp.pima.exper_info.if_num
                for ind in range(if_num):
                    ra_exp.pima.update_cnt(
                        {"BEG_FRQ:": str(ind + 1), "END_FRQ:": str(ind + 1)}
                    )
                    fri = ra_exp.fringe_fitting(
                        bandpass=not args.no_bandpass,
                        accel=not args.no_accel,
                        bandpass_mode=args.bpas_mode,
                        ampl_bandpass=not args.no_ampl_bpas,
                        bandpass_use=args.bpas_use,
                        reference_station=args.ref_sta,
                    )
                    print("IF #{}".format(ind + 1))
                    print(fri)

                # Restore
                ra_exp.pima.update_cnt({"BEG_FRQ:": str(1), "END_FRQ:": str(if_num)})

            fri = ra_exp.fringe_fitting(
                bandpass=not args.no_bandpass,
                accel=not args.no_accel,
                bandpass_mode=args.bpas_mode,
                ampl_bandpass=not args.no_ampl_bpas,
                bandpass_var=args.bpas_var,
                bandpass_use=args.bpas_use,
                bandpass_norm=args.bpas_norm,
                bandpass_renorm=not args.no_bpas_renorm,
                reference_station=args.ref_sta,
            )
            print(fri)

            max_scan_len = fri.max_scan_length()
            logging.debug("DEBUG: max_scan_len = %s", max_scan_len)
            if args.split:
                for aver in (False, True):
                    ra_exp.split(average=aver)

                    # Copy final UV-FITS files to the system tmp directory
                    ra_exp.copy_uvfits(tempfile.gettempdir())

    except PimaError as err:
        return 1
    except RaExperimentError as err:
        return 1
    except psycopg2.Error as err:
        logging.error("DBError: %s", err)
        return 1
    except OSError as err:
        logging.error("OSError: %s", err)
        return 1
    except KeyboardInterrupt:
        logging.warning("KeyboardInterrupt")
        return 1
    except Exception:
        logging.error("Unexpected error: %s", sys.exc_info()[0])
        raise
    finally:
        database.close()

    return 0


if __name__ == "__main__":
    sys.exit(main())
