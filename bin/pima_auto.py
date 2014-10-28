#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 18.02.2014

@author: Petr Voytsik
"""

from __future__ import print_function
import multiprocessing
import os.path
import sys
PATH = os.path.normpath(os.path.join(os.path.dirname(sys.argv[0]), '..'))
sys.path.insert(0, PATH)
import pypima
from pypima.raexperiment import RaExperiment
#import signal
import time


#def init_pool():
#    signal.signal(signal.SIGINT, signal.SIG_IGN)


def download_it(ra_exp):
    """
    Download all necessary files for given experiment.

    """
    my_name = multiprocessing.current_process().name
    print(my_name, 'started ({} {})'.format(ra_exp.exper, ra_exp.band))
    try:
        ra_exp.load(download_only=True)
    except pypima.pima.Error as err:
        print(my_name, 'PIMA Error: ', err)
        return
    except pypima.raexperiment.Error as err:
        print(my_name, 'RaExperiment Error: ', err)
        return
    except KeyboardInterrupt:
        print(my_name, 'KeyboardInterrupt')
        return
    except:
        print(my_name, "Unexpected error: ", sys.exc_info()[0])
        raise
    finally:
        ra_exp.db.close()


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

    with open(in_file_name, 'r') as in_file:
        for line in in_file:
            line = line.strip()
            if line.startswith('#'):
                continue
            exp_band = line.split()
            if len(exp_band) == 2:
                exp_list.append(RaExperiment(exp_band[0], exp_band[1]))

    pool = multiprocessing.Pool(processes=2)
    pool.map_async(download_it, exp_list)

    pool.close()
    time.sleep(1)

    out_dir = os.path.join(os.getenv('HOME'), 'pima_auto_split')
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    for ra_exp in exp_list:
        try:
            ra_exp.load()

            for polar in ['RR', 'RL', 'LR', 'LL']:
                ra_exp.pima.set_polar(polar)
                ra_exp.fringe_fitting(True, True)
                ra_exp.fringes2db()
                ra_exp.split()
                ra_exp.copy_uvfits(out_dir)

            ra_exp.delete_uvfits()
            ra_exp.db.close()
        except pypima.pima.Error as err:
            print('PIMA Error: ', err)
            ra_exp.db.set_error_msg(str(err))
            continue
        except pypima.raexperiment.Error as err:
            print('RaExperiment Error: ', err)
            continue
        except KeyboardInterrupt:
            print('KeyboardInterrupt', file=sys.stderr)
            pool.terminate()
            pool.join()
            return 1
        except:
            print("Unexpected error: ", sys.exc_info()[0])
            raise

    pool.join()
    print("Quitting normally")

    return 0


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: {} <FILE>'.format(os.path.basename(sys.argv[0])),
              file=sys.stderr)
        print('Read list of experiments and bands from FILE', file=sys.stderr)
        sys.exit(2)

    sys.exit(main(sys.argv[1]))
