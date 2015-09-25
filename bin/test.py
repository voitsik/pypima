#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 17:50:20 2013

@author: Petr Voytsik
"""

import argparse
import os.path
import psycopg2
import sys
import tempfile

PATH = os.path.normpath(os.path.join(os.path.dirname(sys.argv[0]), '..'))
sys.path.insert(0, PATH)
import pypima
from pypima.raexperiment import RaExperiment
from pypima.db import DB
from pypima.fri import Fri


def main(args):
    """Main"""
    exper = args.exper.lower()
    band = args.band.lower()
    polar = args.polar

    data_dir = os.getenv('PYPIMA_DATA_DIR',
                         default=os.path.join(os.getenv('HOME'),
                                              'data', 'pima_data'))

    try:
        ra_exp = RaExperiment(exper, band, DB(), gvlbi=args.gvlbi,
                              data_dir=data_dir)
        ra_exp.load(update_db=False)

        if not polar:
            if band == 'l':
                polar = 'RR'
            else:
                polar = 'LL'

        ra_exp.pima.set_polar(polar)
        fri_file = ra_exp.fringe_fitting(True, not args.no_accel)
        fri = Fri(fri_file)
        print(fri)
#        p.fringes2db()
        max_scan_len = fri.max_scan_length()
        print('DEBUG: max_scan_len = ', max_scan_len, file=sys.stderr)
        ra_exp.split(average=0)

        # Copy final UV-FITS files to the system tmp directory
        ra_exp.copy_uvfits(tempfile.gettempdir())
    except pypima.pima.Error as err:
        print('PIMA Error: ', err)
        return 1
    except pypima.raexperiment.Error as err:
        print('RaExperiment Error: ', err)
        return 1
    except psycopg2.Error as err:
        print('DBError: ', err)
        return 1
    except OSError as err:
        print('OSError: ', err)
        return 1
    except KeyboardInterrupt:
        print('KeyboardInterrupt', file=sys.stderr)
        return 1
    except:
        print("Unexpected error: ", sys.exc_info()[0])
        raise

    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('exper', help='experiment code')
    parser.add_argument('band', help='frequency band')
    parser.add_argument('polar', help='polarization', nargs='?')

    # Optional arguments
    parser.add_argument('--gvlbi', '-g', action='store_true',
                        help='process GVLBI FITS-file')
    parser.add_argument('--no-accel', action='store_true',
                        help='disable acceleration term fitting')

    sys.exit(main(parser.parse_args()))
