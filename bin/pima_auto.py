#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 18.02.2014

@author: Petr Voytsik
"""

from __future__ import print_function
from multiprocessing import Pool
import sys
import os.path
path = os.path.normpath(os.path.join(os.path.dirname(sys.argv[0]), '..'))
sys.path.insert(0, path)
import pypima
from pypima.raexperiment import RaExperiment
import time


def download_it(exper):
    """
    Download all necessary files for given experiment.

    """
    try:
        ra_exp = RaExperiment(exper[0], exper[1])
        ra_exp.load(True)
    except pypima.pima.Error as err:
        print('PIMA Error: ', err)
        return
    except pypima.raexperiment.Error as err:
        print('RaExperiment Error: ', err)
        return
#    except KeyboardInterrupt:
#        print('KeyboardInterrupt', file=sys.stderr)
#        exit(1)
    except:
        print("Unexpected error: ", sys.exc_info()[0])
        raise


def main(in_file):
    """
    Main function.

    Parameters:
    -----------
    in_file : str
        Name of the file with list of the experiments and bands. Each line in
        this file must have two words: experiment code and band.

    """
    exp_list = []

    with open(in_file, 'r') as inf:
        for line in inf:
            line = line.strip()
            if line.startswith('#'):
                continue
            exp_band = line.split()
            if len(exp_band) == 2:
                exp_list.append(exp_band)

    pool = Pool(processes=2)
    pool.map_async(download_it, exp_list)

    pool.close()
    time.sleep(1)

    for item in exp_list:
        try:
            p = RaExperiment(item[0], item[1])
            p.load()

            for polar in ['RR', 'RL', 'LR', 'LL']:
                p.pima.set_polar(polar)
                p.fringe_fitting(True, True)
                p.fringes2db()

            p.delete_uvfits()
        except pypima.pima.Error as err:
            print('PIMA Error: ', err)
            p.db.set_error_msg(str(err))
            continue
        except pypima.raexperiment.Error as err:
            print('RaExperiment Error: ', err)
            continue
        except KeyboardInterrupt:
            print('KeyboardInterrupt', file=sys.stderr)
            pool.terminate()
            pool.join()
            exit(1)
        except:
            print("Unexpected error: ", sys.exc_info()[0])
            raise

    pool.join()

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: {} <FILE>'.format(os.path.basename(sys.argv[0])),
              file=sys.stderr)
        print('Read list of experiments and bands from FILE', file=sys.stderr)
        sys.exit(2)

    main(sys.argv[1])
