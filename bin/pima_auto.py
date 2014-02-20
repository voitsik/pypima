#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 18.02.2014

@author: Petr Voytsik
"""

from __future__ import print_function
import sys
import os.path
path = os.path.normpath(os.path.join(os.path.dirname(sys.argv[0]), '..'))
sys.path.insert(0, path)
import pypima
from pypima.raexperiment import RaExperiment


def download_it_all(exp_list):
    """
    Download all necessary files for the list of experiments.

    """
    for item in exp_list:
        try:
            ra_exp = RaExperiment(item[0], item[1])
            ra_exp.load(True)
        except pypima.pima.Error as err:
            print('PIMA Error: ', err)
            sys.exit(1)
        except pypima.raexperiment.Error as err:
            print('RaExperiment Error: ', err)
            sys.exit(1)
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
            exp_list.append(line.split())

    download_it_all(exp_list)

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
            sys.exit(1)
        except pypima.raexperiment.Error as err:
            print('RaExperiment Error: ', err)
            sys.exit(1)
        except:
            print("Unexpected error: ", sys.exc_info()[0])
            raise

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: {} <FILE>'.format(os.path.basename(sys.argv[0])),
              file=sys.stderr)
        print('Read list of experiments and bands from FILE', file=sys.stderr)
        sys.exit(2)

    main(sys.argv[1])
