#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 18.02.2014

@author: Petr Voytsik
"""

from __future__ import print_function
import threading
import os.path
import sys
PATH = os.path.normpath(os.path.join(os.path.dirname(sys.argv[0]), '..'))
sys.path.insert(0, PATH)
import pypima
from pypima.raexperiment import RaExperiment
from pypima.db import DB


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


def main(in_file_name):
    """
    Main function.

    Parameters:
    -----------
    in_file_name : str
        Name of the file with list of the experiments and bands. Each line in
        this file must have two words: experiment code and band.

    """
    exp_list = []

    db = DB()

    with open(in_file_name, 'r') as in_file:
        for line in in_file:
            line = line.strip()
            if line.startswith('#'):
                continue
            exp_band = line.split()
            if len(exp_band) == 2:
                exp_list.append(RaExperiment(exp_band[0], exp_band[1], db))

    load_thread = None

    if len(exp_list) > 1:
        load_thread = threading.Thread(target=download_it,
                                       args=(exp_list[1:],))
        load_thread.daemon = True
        load_thread.start()

    out_dir = os.path.join(os.getenv('HOME'), 'pima_auto_split')
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    spec_out_dir = os.path.join(os.getenv('HOME'), 'public_ftp', 'ra_data',
                                'pima_autospec')

    for ra_exp in exp_list:
        try:
            ra_exp.load(update_db=True)

            for polar in ['RR', 'RL', 'LR', 'LL']:
                ra_exp.pima.set_polar(polar)
                ra_exp.generate_autospectra(spec_out_dir)
                ra_exp.fringe_fitting(True, True)
                ra_exp.fringes2db()
                if ra_exp.pima.chan_number() < 512:
                    for aver in [0, 120, 300, 600, 1200]:
                        ra_exp.split(average=aver)
                        ra_exp.copy_uvfits(out_dir)

            ra_exp.delete_uvfits()
        except pypima.pima.Error as err:
            print('PIMA Error: ', err)
            db.set_error_msg(ra_exp.exper, ra_exp.band, str(err))
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
    if len(sys.argv) != 2:
        print('Usage: {} <FILE>'.format(os.path.basename(sys.argv[0])),
              file=sys.stderr)
        print('Read list of experiments and bands from FILE', file=sys.stderr)
        sys.exit(2)

    sys.exit(main(sys.argv[1]))
