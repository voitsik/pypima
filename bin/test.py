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

    try:
        p = RaExperiment(exper, band, DB())
        p.load(update_db=False, gvlbi=args.g)

        if not polar:
            if band == 'l':
                polar = 'RR'
            else:
                polar = 'LL'

        p.pima.set_polar(polar)
        fri_file = p.fringe_fitting(True, True)
        fri = Fri(fri_file)
        print(fri)
#        p.fringes2db()
        max_scan_len = fri.max_scan_length()
        print('DEBUG: max_scan_len = ', max_scan_len, file=sys.stderr)
        p.split(average=0)
        p.copy_uvfits('/home/voitsik/tmp')
#        p.split(average=max_scan_len)
#        p.copy_uvfits('/home/voitsik/tmp')

#        if p.pima.chan_number() < 512:
#            for part in (2, 3):
#                scan_len = max_scan_len / part
#                p.load(update_db=False, scan_length=scan_len, scan_part=part)
#                fri_file = p.fringe_fitting(True, True)
#                print(Fri(fri_file))
##                p.fringes2db()
#                p.split(average=scan_len)
#                p.copy_uvfits('/home/voitsik/tmp')

    #    p.delete_uvfits()
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
    parser.add_argument('-g', help='Process GVLBI FITS-file',
                        action="store_true")
    args = parser.parse_args()
    sys.exit(main(args))
