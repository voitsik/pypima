#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 18.02.2014

@author: Petr Voytsik
"""

import argparse
import os.path
import psycopg2
import sys
import threading

PATH = os.path.normpath(os.path.join(os.path.dirname(sys.argv[0]), '..'))
sys.path.insert(0, PATH)
import pypima
from pypima.raexperiment import RaExperiment
from pypima.db import DB
from pypima.fri import Fri


def download_it(ra_exps):
    """
    Download FITS files for the list of experiments.

    """
    my_name = 'Download thread'

    print(my_name, ' started')

    for exp in ra_exps:
        try:
            exp.load(download_only=True)
        except pypima.pima.Error as err:
            print(my_name, 'PIMA Error: ', err)
            continue
        except pypima.raexperiment.Error as err:
            print(my_name, 'RaExperiment Error: ', err)
            continue
        except KeyboardInterrupt:
            print(my_name, 'KeyboardInterrupt')
            return
        except:
            print(my_name, "Unexpected error: ", sys.exc_info()[0])
            raise


def main(args):
    """
    Main function.

    Parameters:
    -----------
    args.file_name : str
        Name of the file with list of the experiments and bands. Each line in
        this file must have two words: experiment code and band code.

    """
    exp_list = []

    try:
        db = DB()
    except psycopg2.Error as err:
        print('DBError: ', err, file=sys.stderr)
        return 1

    data_dir = os.getenv('PYPIMA_DATA_DIR',
                         default=os.path.join(os.getenv('HOME'),
                                              'data', 'pima_data'))

    try:
        with open(args.file_name, 'r') as in_file:
            for line in in_file:
                line = line.split('#')[0].strip()
                if not line:
                    continue
                exp_band = line.split()
                if len(exp_band) == 2:
                    exp_list.append(RaExperiment(exp_band[0], exp_band[1], db,
                                                 data_dir=data_dir,
                                                 gvlbi=args.gvlbi))
    except OSError as err:
        print('OSError: ', err, file=sys.stderr)
        return 1

    load_thread = None

    if len(exp_list) > 1:
        load_thread = threading.Thread(target=download_it,
                                       args=(exp_list[1:],))
        load_thread.daemon = True
        load_thread.start()

    out_dir = os.getenv('PYPIMA_SPLIT_DIR',
                        default=os.path.join(os.getenv('HOME'),
                                             'pima_auto_split'))
    # out_dir = os.path.join(os.getenv('HOME'), 'pima_auto_split_new')
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    # Define and create directory for auto spectrum plot files
    spec_out_dir = os.getenv('PYPIMA_AUTOSPEC_DIR',
                             default=os.path.join(os.getenv('HOME'),
                                                  'pima_autospec'))
    if not os.path.exists(spec_out_dir):
        os.mkdir(spec_out_dir)

    for ra_exp in exp_list:
        try:
            # First run on full scan
            ra_exp.load(update_db=True, scan_part=1)

            for polar in ('RR', 'RL', 'LR', 'LL'):
                ra_exp.pima.set_polar(polar)
                if not args.gvlbi:
                    ra_exp.generate_autospectra(spec_out_dir)
                fri_file = ra_exp.fringe_fitting(True, True)
                fri = Fri(fri_file)
                print(fri)
                max_scan_len = fri.max_scan_length()
                ra_exp.fringes2db()

                if ra_exp.pima.chan_number() < 512 and \
                        ra_exp.calibration_loaded and not args.gvlbi:
                    for aver in (0, round(max_scan_len)):
                        ra_exp.split(average=aver)
                        ra_exp.copy_uvfits(out_dir)

            if args.gvlbi:
                continue

            # Second run on half scan
            scan_len = round(max_scan_len/2)
            ra_exp.load(update_db=True, scan_length=scan_len, scan_part=2)

            for polar in ('RR', 'RL', 'LR', 'LL'):
                ra_exp.pima.set_polar(polar)
                fri_file = ra_exp.fringe_fitting(True, True)
                fri = Fri(fri_file)
                print(fri)
                ra_exp.fringes2db()

                if ra_exp.pima.chan_number() < 512 and \
                        ra_exp.calibration_loaded:
                    ra_exp.split(average=scan_len)
                    ra_exp.copy_uvfits(out_dir)

            # For good experiments more runs
            if ra_exp.pima.chan_number() < 512 and ra_exp.calibration_loaded:
                for scan_part in (3, 4, 5):
                    scan_len = round(max_scan_len / scan_part)
                    ra_exp.load(update_db=False, scan_length=scan_len,
                                scan_part=scan_part)
                    for polar in ('RR', 'LL'):
                        ra_exp.pima.set_polar(polar)
                        fri_file = ra_exp.fringe_fitting(True, True)
                        fri = Fri(fri_file)
                        print(fri)
                        ra_exp.split(average=scan_len)
                        ra_exp.copy_uvfits(out_dir)

            ra_exp.delete_uvfits()
        except pypima.pima.Error as err:
            print('PIMA Error: ', err)
            db.set_error_msg(ra_exp.run_id, str(err))
            continue
        except pypima.raexperiment.Error as err:
            print('RaExperiment Error: ', err)
            continue
        except KeyboardInterrupt:
            print('KeyboardInterrupt', file=sys.stderr)
            return 1
        except:
            print("Unexpected error: ", sys.exc_info()[0])
            raise

    if load_thread:
        load_thread.join()

    print("Quitting normally")

    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file_name', metavar='FILE',
                        help='File with list of experiments and bands')
    parser.add_argument('--gvlbi', '-g', action='store_true',
                        help='process ground-only part of the experiments')
    args = parser.parse_args()
    sys.exit(main(args))
