#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 17:50:20 2013

@author: Petr Voytsik
"""

import argparse
import logging
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
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(name)s: %(message)s',
                        level=logging.INFO)

    exper = args.exper.lower()
    band = args.band.lower()
    polar = args.polar

    data_dir = os.getenv('PYPIMA_DATA_DIR',
                         default=os.path.join(os.getenv('HOME'),
                                              'data', 'pima_data'))

    try:
        ra_exp = RaExperiment(exper, band, DB(), gvlbi=args.gvlbi,
                              data_dir=data_dir)
        ra_exp.init_workdir()
        ra_exp.load(update_db=False, force_small=args.force_small)

        if not polar:
            if band == 'l':
                polar = 'RR'
            else:
                polar = 'LL'

        ra_exp.pima.set_polar(polar)

        if args.autospec_only:
            # Define and create directory for auto spectrum plot files
            spec_out_dir = os.getenv('PYPIMA_AUTOSPEC_DIR',
                                     default=os.path.join(os.getenv('HOME'),
                                                          'pima_autospec'))
            if not os.path.exists(spec_out_dir):
                os.mkdir(spec_out_dir)

            ra_exp.generate_autospectra(spec_out_dir)
        else:
            ra_exp.load_antab()

            if args.individual_ifs:
                if_num = ra_exp.pima.exper_info['if_num']
                for ind in range(if_num):
                    ra_exp.pima.update_cnt({'BEG_FRQ:': str(ind+1),
                                            'END_FRQ:': str(ind+1)})
                    fri_file = ra_exp.fringe_fitting(True, not args.no_accel)
                    print('IF #{}'.format(ind+1))
                    print(Fri(fri_file))

                # Restore
                ra_exp.pima.update_cnt({'BEG_FRQ:': str(1),
                                        'END_FRQ:': str(if_num)})

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
        return 1
    except pypima.raexperiment.Error as err:
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
    parser.add_argument('--force-small', action='store_true',
                        help='force to use 64-channel FITS file (if any)')
    parser.add_argument('--autospec-only', action='store_true',
                        help='generate autocorrelation spectra only')
    parser.add_argument('--individual-ifs', action='store_true',
                        help='do fringe fittig for individual IFs')

    sys.exit(main(parser.parse_args()))
