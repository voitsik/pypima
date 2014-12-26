#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 17:50:20 2013

@author: Petr Voytsik
"""

from __future__ import print_function
import sys
import os.path
PATH = os.path.normpath(os.path.join(os.path.dirname(sys.argv[0]), '..'))
sys.path.insert(0, PATH)
import pypima
from pypima.raexperiment import RaExperiment
from pypima.db import DB


def main():
    """Main"""
    if len(sys.argv) != 3:
        print('Usage: {} <exper> <band>'.format(sys.argv[0]), file=sys.stderr)
        sys.exit(2)

    exper = sys.argv[1].lower()
    band = sys.argv[2].lower()

    try:
        p = RaExperiment(exper, band, DB())
        p.load()

        polar = 'LL'
        if band == 'l':
            polar = 'RR'

        p.pima.set_polar(polar)
        p.fringe_fitting(True, True)

        spec_out_dir = os.path.join(os.getenv('HOME'), 'public_ftp', 'ra_data',
                                    'pima_autospec')
        if not os.path.isdir(spec_out_dir):
            os.mkdir(spec_out_dir)

        p.generate_autospectra(spec_out_dir)

        p.split()

    #    p.delete_uvfits()
    except pypima.pima.Error as err:
        print('PIMA Error: ', err)
        sys.exit(1)
    except pypima.raexperiment.Error as err:
        print('RaExperiment Error: ', err)
        sys.exit(1)
    except KeyboardInterrupt:
        print('KeyboardInterrupt', file=sys.stderr)
        exit(1)
    except:
        print("Unexpected error: ", sys.exc_info()[0])
        raise


if __name__ == "__main__":
    main()
