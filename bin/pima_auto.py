#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 18.02.2014

@author: Petr Voytsik
"""

import argparse
import logging
import os.path
import psycopg2
import shutil
import sys
import threading

PATH = os.path.normpath(os.path.join(os.path.dirname(sys.argv[0]), '..'))
sys.path.insert(0, PATH)
import pypima
from pypima.raexperiment import RaExperiment
from pypima.db import DB


def download_it(ra_exps, force_small):
    """
    Download FITS files for the list of the experiments.

    """
    logger = logging.getLogger('download_thread')

    logger.info('started')

    for exp in ra_exps:
        # Hack: make data directory for shutil.disk_usage
        os.makedirs(exp.data_dir, exist_ok=True)
        df = shutil.disk_usage(exp.data_dir)
        if df.free / df.total < 0.1:
            logger.warning('less then 10% of free space left on disk')
            return

        try:
            exp.load(download_only=True, force_small=force_small)
        except pypima.pima.Error:
            continue
        except pypima.raexperiment.Error:
            continue
        except KeyboardInterrupt:
            logger.warn('KeyboardInterrupt')
            return
        except:
            logger.error('Unexpected error: %s', sys.exc_info()[0])
            raise


def generate_autospec(ra_exp, spec_out_dir, force_small=False):
    """
    Run ``PIMA`` only for average autocorrelation spectrum computing by
    ``acta`` task.

    Parameters
    ----------
    ra_exp : RaExperiment object
        Experiment setup.
    spec_out_dir : str
        Directory name for autocorrelation plots.
    force_small : bool
        Use FITS-IDI with 64 spectarl channels only.

    """
    ra_exp.load(update_db=False, scan_part=1, force_small=force_small)

    for polar in ('RR', 'RL', 'LR', 'LL'):
        ra_exp.pima.set_polar(polar)
        ra_exp.generate_autospectra(plot=True, out_dir=spec_out_dir, db=True)

    ra_exp.delete_uvfits()


def process_gvlbi(ra_exp, accel=False, force_small=False):
    """
    Process ground-only part of the experiment.

    """
    ra_exp.load(update_db=True, scan_part=1, force_small=force_small)

    for polar in ('RR', 'RL', 'LR', 'LL'):
        ra_exp.pima.set_polar(polar)
        print(ra_exp.fringe_fitting(True, accel))
        ra_exp.fringes2db()

    ra_exp.delete_uvfits()


def process_ind_ifs(ra_exp, accel=False, force_small=False):
    """
    Do fringe fitting for each IF separatly.

    """
    ra_exp.load(update_db=True, scan_part=-1, force_small=force_small)

    for polar in ('RR', 'RL', 'LR', 'LL'):
        ra_exp.pima.set_polar(polar)
        if_num = ra_exp.pima.exper_info['if_num']

        for ind in range(if_num):
            ra_exp.pima.update_cnt({'BEG_FRQ:': str(ind+1),
                                    'END_FRQ:': str(ind+1)})
            fri = ra_exp.fringe_fitting(True, accel)
            print('IF #{}'.format(ind+1))
            print(fri)
            ra_exp.fringes2db()

        # Restore IFs
        ra_exp.pima.update_cnt({'BEG_FRQ:': str(1),
                                'END_FRQ:': str(if_num)})

    ra_exp.delete_uvfits()


def process_radioastron(ra_exp, uv_fits_out_dir, spec_out_dir, accel=True,
                        force_small=False, alt=False):
    """
    Process space-ground part of the experiment.

    Parameters
    ----------
    ra_exp : RaExperiment object
        Experiment setup.
    uv_fits_out_dir : str
        Directory name for final UV-FITS files.
    spec_out_dir : str
        Directory name for autocorrelation plots.
    accel : bool
        If True, do the parabolic term fitting.
    force_small : bool
        Use FITS-IDI with 64 spectral channels only.

    """
    if alt:
        scan_part_base = 1000
    else:
        scan_part_base = 0

    # First run on full scan
    scan_part = scan_part_base + 1
    ra_exp.load(update_db=True, scan_part=scan_part, force_small=force_small)
    ra_exp.load_antab()

    scan_len_list = []
    for polar in ('RR', 'RL', 'LR', 'LL'):
        ra_exp.pima.set_polar(polar)
        ra_exp.generate_autospectra(plot=True, out_dir=spec_out_dir, db=True)
        fri = ra_exp.fringe_fitting(True, accel)
        print(fri)
        scan_len_list.append(fri.max_scan_length())
        ra_exp.fringes2db()

        if ra_exp.pima.chan_number() <= 128 and ra_exp.calibration_loaded and \
                polar in ('RR', 'LL'):
            for aver in (0, round(scan_len_list[-1])):
                ra_exp.split(average=aver)
                ra_exp.copy_uvfits(uv_fits_out_dir)

    # Second run on a scan half
    max_scan_len = max(scan_len_list)
    scan_len = round(max_scan_len / 2)
    scan_part = scan_part_base + 2
    ra_exp.load(update_db=True, scan_length=scan_len, scan_part=scan_part,
                force_small=force_small)
    ra_exp.load_antab()

    detections = False
    for polar in ('RR', 'RL', 'LR', 'LL'):
        ra_exp.pima.set_polar(polar)
        fri = ra_exp.fringe_fitting(True, accel)
        print(fri)
        if fri.any_detections():
            detections = True
        ra_exp.fringes2db()

        if ra_exp.pima.chan_number() <= 128 and ra_exp.calibration_loaded and \
                polar in ('RR', 'LL'):
            ra_exp.split(average=scan_len)
            ra_exp.copy_uvfits(uv_fits_out_dir)

    # For good experiments more runs
    if ra_exp.pima.chan_number() <= 128 and ra_exp.calibration_loaded:
        scan_part = scan_part_base + 3
        scan_len = round(max_scan_len / (scan_part-scan_part_base))

        while scan_len > 100 and detections:
            ra_exp.load(update_db=True, scan_length=scan_len,
                        scan_part=scan_part, force_small=force_small)
            ra_exp.load_antab()
            detections = False

            for polar in ('RR', 'LL'):
                ra_exp.pima.set_polar(polar)
                fri = ra_exp.fringe_fitting(True, accel)
                print(fri)
                if fri.any_detections():
                    detections = True
                ra_exp.fringes2db()
                ra_exp.split(average=scan_len)
                ra_exp.copy_uvfits(uv_fits_out_dir)

            scan_part += 1
            scan_len = round(max_scan_len / (scan_part-scan_part_base))

    # Special run for ground-ground baselines with 60 s scan length
    if ra_exp.pima.chan_number() <= 128 and ra_exp.calibration_loaded and \
            detections:
        scan_part = scan_part_base + 100
        scan_len = 60
        ra_exp.load(update_db=True, scan_length=scan_len,
                    scan_part=scan_part, force_small=force_small)
        ra_exp.load_antab()

        for polar in ('RR', 'LL'):
            ra_exp.pima.set_polar(polar)
            fri = ra_exp.fringe_fitting(True, accel)
            print(fri)
            ra_exp.fringes2db()
            ra_exp.split(average=scan_len)
            ra_exp.copy_uvfits(uv_fits_out_dir)

    ra_exp.delete_uvfits()


def main(args):
    """
    Main function.

    Parameters:
    -----------
    args.file_name : str
        Name of the file with list of the experiments and bands. Each line in
        this file must have two words: experiment code and band code.

    """
    log_format = '%(asctime)s %(levelname)s: %(name)s: %(message)s'
    logging.basicConfig(format=log_format,
                        level=logging.INFO, filename=args.log_file)
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)

    exp_list = []

    try:
        database = DB()
    except psycopg2.Error as err:
        logging.error('DBError: %s', err)
        return 1

    data_dir = os.getenv('PYPIMA_DATA_DIR',
                         default=os.path.join(os.getenv('HOME'),
                                              'data', 'pima_data'))

    # Parse input file.
    # The input file is a table with three columns: exper_name, band, and
    # FITS-IDI file name (optional)
    try:
        with open(args.file_name, 'r') as in_file:
            for line in in_file:
                line = line.split('#')[0].strip()

                if not line:
                    continue

                exp_band = line.split()

                if len(exp_band) in (2, 3):
                    if len(exp_band) == 3:
                        fits = exp_band[2]
                    else:
                        fits = None

                    exp_list.append(RaExperiment(exp_band[0], exp_band[1],
                                                 database, uv_fits=fits,
                                                 data_dir=data_dir,
                                                 gvlbi=args.gvlbi))
    except OSError as err:
        logging.error('OSError: %s', err)
        return 1

    out_dir = os.getenv('PYPIMA_SPLIT_DIR',
                        default=os.path.join(os.getenv('HOME'),
                                             'pima_auto_split'))

    # Define and create directory for auto spectrum plot files
    spec_out_dir = os.getenv('PYPIMA_AUTOSPEC_DIR',
                             default=os.path.join(os.getenv('HOME'),
                                                  'pima_autospec'))

    load_thread = threading.Thread(target=download_it, args=(exp_list,
                                                             args.force_small))
    load_thread.daemon = True
    load_thread.start()

    for ra_exp in exp_list:
        try:
            ra_exp.init_workdir()
            if args.autospec_only:
                generate_autospec(ra_exp, spec_out_dir, args.force_small)
            elif args.individual_ifs:
                process_ind_ifs(ra_exp, not args.no_accel, args.force_small)
            else:
                if ra_exp.gvlbi:
                    process_gvlbi(ra_exp, not args.no_accel, args.force_small)
                else:
                    process_radioastron(ra_exp, out_dir, spec_out_dir,
                                        accel=not args.no_accel,
                                        force_small=args.force_small,
                                        alt=args.alt)
        except pypima.pima.Error as err:
            database.set_error_msg(ra_exp.run_id, str(err))
            ra_exp.delete_uvfits()
            continue
        except pypima.raexperiment.Error as err:
            continue
        except psycopg2.Error as err:
            logging.error('DBError: %s', err)
            return 1
        except OSError as err:
            logging.error('OSError: %s', err)
            return 1
        except KeyboardInterrupt:
            logging.warning('KeyboardInterrupt')
            return 1
        except Exception:
            logging.error("Unexpected error: %s", sys.exc_info()[0])
            raise

    load_thread.join()

    logging.info("Quitting normally")

    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file_name', metavar='FILE',
                        help='File with list of experiments and bands')

    # Optional arguments
    parser.add_argument('-l', '--log-file', metavar='LOG',
                        help='log file')
    parser.add_argument('--gvlbi', '-g', action='store_true',
                        help='process ground-only part of the experiments')
    parser.add_argument('--no-accel', action='store_true',
                        help='disable parabolic term fitting')
    parser.add_argument('--force-small', action='store_true',
                        help='force to use 64-channel FITS file (if any)')
    parser.add_argument('--alt', action='store_true',
                        help='use alternative scan_part_base')

    group = parser.add_mutually_exclusive_group()
    group.add_argument('--autospec-only', action='store_true',
                       help='generate autocorrelation spectra only')
    group.add_argument('--individual-ifs', action='store_true',
                       help='do fringe fittig for individual IFs')
    parser.add_argument('--debug', '-d', action='store_true',
                        help='enable debug output')

    sys.exit(main(parser.parse_args()))
