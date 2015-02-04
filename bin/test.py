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
from pypima.fri import Fri


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
        fri_file = p.fringe_fitting(True, True)
        fri = Fri(fri_file)
        print(fri)
        max_scan_len = fri.max_scan_length()
        print('DEBUG: max_scan_len = ', max_scan_len, file=sys.stderr)
        p.split(average=max_scan_len)
        p.copy_uvfits('/home/voitsik/tmp')

        if p.pima.chan_number() < 512:
            for scan_len in (round(max_scan_len/2),
                             round(max_scan_len/3), round(max_scan_len/4)):
                p.load(update_db=False, scan_length=scan_len)
                fri_file = p.fringe_fitting(True, True)
                print(Fri(fri_file))
                p.split(average=scan_len)
                p.copy_uvfits('/home/voitsik/tmp')

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
