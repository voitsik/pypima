#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 18.02.2014

@author: Petr Voytsik
"""

import argparse
import logging
import os.path
import shutil
import sys
import threading

import psycopg2

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
        disk_usage = shutil.disk_usage(exp.data_dir)

        if disk_usage.free / disk_usage.total < 0.1:
            logger.warning('less then 10% of free space left on disk')
            return

        try:
            exp.load(download_only=True, force_small=force_small)
        except pypima.pima.Error:
            continue
        except pypima.raexperiment.Error:
            continue
        except KeyboardInterrupt:
            logger.warning('KeyboardInterrupt')
            return
        except Exception:
            logger.error('Unexpected error: %s', sys.exc_info()[0])
            raise


def generate_autospec(ra_exp, spec_out_dir, force_small=False):
    """
    Run **PIMA** only for average autocorrelation spectrum computing by
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


def process_gvlbi(ra_exp, **kwargs):
    """
    Process ground-only part of the experiment.

    """
    accel = kwargs.pop('accel', False)
    force_small = kwargs.pop('force_small', False)
    # scan_part_base = kwargs.pop('scan_part_base', 0)
    reference_station = kwargs.pop('ref_sta', None)
    bandpass_mode = kwargs.pop('bandpass_mode', None)
    ampl_bandpass = kwargs.pop('ampl_bandpass', True)
    bandpass_var = kwargs.pop('bandpass_var', 0)
    bandpass_use = kwargs.pop('bandpass_use', None)
    flag_chann = kwargs.pop('flag_chann', 0)
    scan_len = kwargs.pop('max_scan_length', 1500)

    ff_params = {
        'bandpass': True,
        'accel': accel,
        'reference_station': reference_station,
        'bandpass_mode': bandpass_mode,
        'ampl_bandpass': ampl_bandpass,
        'bandpass_var': bandpass_var,
        'bandpass_use': bandpass_use,
        }

    ra_exp.load(update_db=True, scan_length=scan_len, scan_part=1,
                force_small=force_small)
    ra_exp.flag_edge_chann(flag_chann)

    for polar in ('RR', 'RL', 'LR', 'LL'):
        ra_exp.pima.set_polar(polar)
        print(ra_exp.fringe_fitting(**ff_params))
        ra_exp.fringes2db()

    ra_exp.delete_uvfits()


def process_ind_ifs(ra_exp, **kwargs):
    """
    Do fringe fitting for each IF separatly.

    """
    accel = kwargs.pop('accel', True)
    force_small = kwargs.pop('force_small', False)
    # scan_part_base = kwargs.pop('scan_part_base', 0)
    reference_station = kwargs.pop('ref_sta', None)
    bandpass_mode = kwargs.pop('bandpass_mode', None)
    ampl_bandpass = kwargs.pop('ampl_bandpass', True)
    bandpass_var = kwargs.pop('bandpass_var', 0)
    bandpass_use = kwargs.pop('bandpass_use', None)
    flag_chann = kwargs.pop('flag_chann', 0)
    scan_len = kwargs.pop('max_scan_length', 1500)

    ff_params = {
        'bandpass': True,
        'accel': accel,
        'reference_station': reference_station,
        'bandpass_mode': bandpass_mode,
        'ampl_bandpass': ampl_bandpass,
        'bandpass_var': bandpass_var,
        'bandpass_use': bandpass_use,
        }

    ra_exp.load(update_db=True, scan_length=scan_len, scan_part=-1,
                force_small=force_small)
    ra_exp.flag_edge_chann(flag_chann)

    for polar in ('RR', 'RL', 'LR', 'LL'):
        ra_exp.pima.set_polar(polar)
        if_num = ra_exp.pima.exper_info['if_num']

        for ind in range(if_num):
            ra_exp.pima.update_cnt({'BEG_FRQ:': str(ind+1),
                                    'END_FRQ:': str(ind+1)})
            fri = ra_exp.fringe_fitting(**ff_params)
            print('IF #{}'.format(ind+1))
            print(fri)
            ra_exp.fringes2db()

        # Restore IFs
        ra_exp.pima.update_cnt({'BEG_FRQ:': str(1),
                                'END_FRQ:': str(if_num)})

    ra_exp.delete_uvfits()


def process_radioastron(ra_exp, uv_fits_out_dir, spec_out_dir, **kwargs):
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
    accel : bool, optional
        If True, do the parabolic term fitting.
    bandpass_mode : str, optional
        Set the ``BPS.MODE`` **PIMA** parameter. If ``None``, a value from the
        control file is used.
    force_small : bool, optional
        Use FITS-IDI with 64 spectral channels only.
    scan_part_base : int, optional
        Default value is 0.
    bandpass_var : int, optional
        Default value is 0.
    flag_chann : int, optional
        Flag edge spectral channels.
    reference_station : str, optional
        Reference station for bandpass calibration. Default value is ``None``.

    """
    accel = kwargs.pop('accel', True)
    force_small = kwargs.pop('force_small', False)
    scan_part_base = kwargs.pop('scan_part_base', 0)
    reference_station = kwargs.pop('ref_sta', None)
    bandpass_mode = kwargs.pop('bandpass_mode', None)
    ampl_bandpass = kwargs.pop('ampl_bandpass', True)
    bandpass_var = kwargs.pop('bandpass_var', 0)
    bandpass_use = kwargs.pop('bandpass_use', None)
    flag_chann = kwargs.pop('flag_chann', 0)
    scan_len = kwargs.pop('max_scan_length', 1500)

    # First run on full scan
    scan_part = scan_part_base + 1
    ra_exp.load(update_db=True, scan_length=scan_len, scan_part=scan_part,
                force_small=force_small)
    ra_exp.flag_edge_chann(flag_chann)
    ra_exp.load_antab()

    ff_params = {
        'bandpass': True,
        'accel': accel,
        'reference_station': reference_station,
        'bandpass_mode': bandpass_mode,
        'ampl_bandpass': ampl_bandpass,
        'bandpass_var': bandpass_var,
        'bandpass_use': bandpass_use,
        }

    scan_len_list = []
    for polar in ('RR', 'RL', 'LR', 'LL'):
        ra_exp.pima.set_polar(polar)
        if scan_part == 1:
            ra_exp.generate_autospectra(plot=True, out_dir=spec_out_dir,
                                        db=True)
        fri = ra_exp.fringe_fitting(**ff_params)
        print(fri)
        scan_len_list.append(fri.max_scan_length())
        ra_exp.fringes2db()

        if ra_exp.pima.chan_number() <= 128 and ra_exp.calibration_loaded and \
                polar in ('RR', 'LL'):
            for aver in (0, round(scan_len_list[-1])):
                ra_exp.split(average=aver)
                ra_exp.copy_uvfits(uv_fits_out_dir)

    #
    # Second run on a scan half
    #
    max_scan_len = max(scan_len_list)
    scan_len = round(max_scan_len / 2)
    scan_part = scan_part_base + 2
    ra_exp.load(update_db=True, scan_length=scan_len, scan_part=scan_part,
                force_small=force_small)
    ra_exp.load_antab()

    detections = False
    for polar in ('RR', 'RL', 'LR', 'LL'):
        ra_exp.pima.set_polar(polar)
        fri = ra_exp.fringe_fitting(**ff_params)
        print(fri)
        if fri.any_detections():
            detections = True
        ra_exp.fringes2db()

        if ra_exp.pima.chan_number() <= 128 and ra_exp.calibration_loaded and \
                polar in ('RR', 'LL'):
            ra_exp.split(average=scan_len)
            ra_exp.copy_uvfits(uv_fits_out_dir)

    #
    # For good experiments more runs
    #
    if ra_exp.pima.chan_number() <= 128:
        scan_part = scan_part_base + 3
        scan_len = round(max_scan_len / (scan_part-scan_part_base))

        while scan_len > 100 and detections:
            ra_exp.load(update_db=True, scan_length=scan_len,
                        scan_part=scan_part, force_small=force_small)
            ra_exp.load_antab()
            detections = False

            for polar in ('RR', 'LL'):
                ra_exp.pima.set_polar(polar)
                fri = ra_exp.fringe_fitting(**ff_params)
                print(fri)
                if fri.any_detections():
                    detections = True
                ra_exp.fringes2db()

                if ra_exp.calibration_loaded:
                    ra_exp.split(average=scan_len)
                    ra_exp.copy_uvfits(uv_fits_out_dir)

            scan_part += 1
            scan_len = round(max_scan_len / (scan_part-scan_part_base))

    # Special run for ground-ground baselines with 60 s scan length
    if ra_exp.pima.chan_number() <= 128 and detections:
        scan_part = scan_part_base + 100
        scan_len = 60
        ra_exp.load(update_db=True, scan_length=scan_len,
                    scan_part=scan_part, force_small=force_small)
        ra_exp.load_antab()

        for polar in ('RR', 'LL'):
            ra_exp.pima.set_polar(polar)
            fri = ra_exp.fringe_fitting(**ff_params)
            print(fri)
            ra_exp.fringes2db()

            if ra_exp.calibration_loaded:
                ra_exp.split(average=scan_len)
                ra_exp.copy_uvfits(uv_fits_out_dir)

    ra_exp.delete_uvfits()


def parser_input_file(file_name, database, data_dir, gvlbi):
    """
    Parse input file.

    The input file is a table with three columns: exper_name, band, and
    FITS-IDI file name (optional).

    Returns
    -------
    exp_list : list
        List of RaExperiment instances.

    """
    exp_list = []

    with open(file_name, 'r') as file:
        for line in file:
            line = line.split('#')[0].strip()  # Strip comments out

            columns = line.split()

            if len(columns) in (2, 3):
                exper_name = columns[0]
                band = columns[1]

                if len(columns) == 3:
                    fits = columns[2]
                else:
                    fits = None

                exp_list.append(RaExperiment(exper_name, band,
                                             database, uv_fits=fits,
                                             data_dir=data_dir,
                                             gvlbi=gvlbi))

    return exp_list


def parse_args():
    """
    Parse command line arguments.

    """
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
                        help='force to use 64-channel FITS file')
    parser.add_argument('--scan-part-base', type=int, default=0, metavar='ALT',
                        choices=[1000 * x for x in range(20)],
                        help='use alternative scan_part_base')
    parser.add_argument('--ref-sta', metavar='STA',
                        help='reference station for bandpass calibration')
    parser.add_argument('--bpas-mode', metavar='MODE',
                        choices=['INIT', 'ACCUM', 'FINE'],
                        help='bandpass calibration mode')
    parser.add_argument('--bpas-use', metavar='BANDPASS_USE', default='PHS',
                        choices=['AMP', 'PHS', 'AMP_PHS', 'NO'],
                        help='BANDPASS_USE PIMA parameter (default is PHS)')
    parser.add_argument('--no-ampl-bpas', action='store_true',
                        help='disable amplitude bandpass calibration')
    parser.add_argument('--bpas-var', type=int, choices=[0, 1, 2, 3],
                        default=3,
                        help='predefined bandpass parameters (default is 3)')
    parser.add_argument('--flag-chann', type=int, default=2, metavar='N',
                        help='flag N edge spectral channels of the bandpass '
                        '(default is 2)')
    parser.add_argument('--max-scan-length', type=float, default=1500,
                        help='maximum scan length in seconds')

    group = parser.add_mutually_exclusive_group()
    group.add_argument('--autospec-only', action='store_true',
                       help='generate autocorrelation spectra only')
    group.add_argument('--individual-ifs', action='store_true',
                       help='do fringe fittig for individual IFs')

    parser.add_argument('--debug', '-d', action='store_true',
                        help='enable debug output')

    return parser.parse_args()


def main():
    """Main"""
    args = parse_args()

    if not args.log_file:
        log_file = args.file_name.rsplit('.', 1)[0] + '.log'
    else:
        log_file = args.log_file

    log_format = '%(asctime)s %(levelname)s: %(name)s: %(message)s'
    logging.basicConfig(format=log_format,
                        level=logging.INFO, filename=log_file)
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)

    # Check input parameters
    if args.max_scan_length <= 0:
        logging.error('Max scan length must be positive')
        return 1

    if args.flag_chann < 0:
        logging.error('Number of flagged channels must not be negative')

    # Connect to database
    try:
        database = DB()
    except psycopg2.Error as err:
        logging.error('DBError: %s', err)
        return 1

    data_dir = os.getenv('PYPIMA_DATA_DIR',
                         default=os.path.join(os.getenv('HOME'),
                                              'data', 'pima_data'))

    try:
        exp_list = parser_input_file(args.file_name, database, data_dir,
                                     args.gvlbi)
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

    params = {
        'scan_part_base': args.scan_part_base,
        'force_small': args.force_small,
        'accel': not args.no_accel,
        'ref_sta': args.ref_sta,
        'bandpass_mode': args.bpas_mode,
        'ampl_bandpass': not args.no_ampl_bpas,
        'bandpass_var': args.bpas_var,
        'bandpass_use': args.bpas_use,
        'flag_chann': args.flag_chann,
        'max_scan_length': args.max_scan_length,
        }

    for ra_exp in exp_list:
        try:
            ra_exp.init_workdir()

            if args.autospec_only:
                generate_autospec(ra_exp, spec_out_dir, args.force_small)
            elif args.individual_ifs:
                process_ind_ifs(ra_exp, **params)
            else:
                if ra_exp.gvlbi:
                    process_gvlbi(ra_exp, **params)
                else:
                    process_radioastron(ra_exp, out_dir, spec_out_dir,
                                        **params)
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
    sys.exit(main())
