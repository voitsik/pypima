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


def main(args):
    """Main"""
    log_format = '%(asctime)s %(levelname)s: %(name)s: %(message)s'
    logging.basicConfig(format=log_format, level=logging.INFO)

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)

    exper = args.exper.lower()
    band = args.band.lower()
    polar = args.polar

    data_dir = os.getenv('PYPIMA_DATA_DIR',
                         default=os.path.join(os.getenv('HOME'),
                                              'data', 'pima_data'))

    try:
        ra_exp = RaExperiment(exper, band, DB(), gvlbi=args.gvlbi,
                              data_dir=data_dir, uv_fits=args.fits)
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

            ra_exp.generate_autospectra(plot=True, out_dir=spec_out_dir,
                                        db=True)
        else:
            ra_exp.load_antab()

            if args.individual_ifs:
                if_num = ra_exp.pima.exper_info['if_num']
                for ind in range(if_num):
                    ra_exp.pima.update_cnt({'BEG_FRQ:': str(ind+1),
                                            'END_FRQ:': str(ind+1)})
                    fri = ra_exp.fringe_fitting(bandpass=True,
                                                accel=not args.no_accel,
                                                bandpass_mode=args.bpas_mode,
                                                ampl_bandpass=not args.no_ampl_bpas)
                    print('IF #{}'.format(ind+1))
                    print(fri)

                # Restore
                ra_exp.pima.update_cnt({'BEG_FRQ:': str(1),
                                        'END_FRQ:': str(if_num)})

            fri = ra_exp.fringe_fitting(bandpass=True, accel=not args.no_accel,
                                        bandpass_mode=args.bpas_mode,
                                        ampl_bandpass=not args.no_ampl_bpas)
            print(fri)

            max_scan_len = fri.max_scan_length()
            logging.debug('DEBUG: max_scan_len = %s', max_scan_len)
            if args.split:
                ra_exp.split(average=0)

                # Copy final UV-FITS files to the system tmp directory
                ra_exp.copy_uvfits(tempfile.gettempdir())
    except pypima.pima.Error as err:
        return 1
    except pypima.raexperiment.Error as err:
        return 1
    except psycopg2.Error as err:
        logging.error('DBError: %s', err)
        return 1
    except OSError as err:
        logging.error('OSError: %s', err)
        return 1
    except KeyboardInterrupt:
        logging.warning('KeyboardInterrupt')
        return 1
    except Exception:
        logging.error('Unexpected error: %s', sys.exc_info()[0])
        raise

    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('exper', help='experiment code')
    parser.add_argument('band', help='frequency band')

    # Optional arguments
    parser.add_argument('--polar', help='polarization')
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
    parser.add_argument('--split', action='store_true',
                        help='do SPLIT')
    parser.add_argument('--fits',
                        help='external FITS-IDI file')
    parser.add_argument('--bpas-mode', metavar='MODE',
                        choices=['INIT', 'ACCUM', 'FINE'],
                        help='set bandpass calibration mode')
    parser.add_argument('--no-ampl-bpas', action='store_true',
                        help='disable amplitude bandpass calibration')
    parser.add_argument('--debug', '-d', action='store_true',
                        help='enable debug output')

    sys.exit(main(parser.parse_args()))
