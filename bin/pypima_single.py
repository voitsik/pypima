#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 17:50:20 2013

@author: Petr Voytsik
"""

import argparse
import logging
import os.path
import sys
import tempfile

import psycopg2

from pypima import DB, PimaError, RaExperiment, RaExperimentError


def parse_args():
    """Parse command line arguments"""
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
    parser.add_argument('--fits', nargs='+',
                        help='external FITS-IDI file(s)')
    parser.add_argument('--orbit',
                        help='reconstructed orbit file')
    parser.add_argument('--scan-length', type=float,
                        help='set scan length in seconds')
    parser.add_argument('--frequency-group', type=int, default=1,
                        help='frequency group')
    parser.add_argument('--debug', '-d', action='store_true',
                        help='enable debug output')

    group = parser.add_argument_group("bandpass settings")
    group.add_argument('--ref-sta', metavar='STA',
                       help='reference station')
    group.add_argument('--no-bandpass', action='store_true',
                       help='disable bandpass calibration')
    group.add_argument('--bpas-mode', metavar='MODE',
                       choices=['INIT', 'ACCUM', 'FINE'],
                       help='set bandpass calibration mode')
    group.add_argument('--bpas-use', metavar='BANDPASS_USE', default='PHS',
                       choices=['AMP', 'PHS', 'AMP_PHS', 'NO'],
                       help='set BANDPASS_USE PIMA parameter')
    group.add_argument('--no-ampl-bpas', action='store_true',
                       help='disable amplitude bandpass calibration')
    group.add_argument('--bpas-var', type=int, choices=[0, 1, 2, 3, 4, 5],
                       default=3,
                       help='predefined bandpass parameters')
    group.add_argument('--flag-chann', type=int, default=2, metavar='N',
                       help='flag N edge spectral channels of the bandpass')

    return parser.parse_args()


def main():
    """Main"""
    args = parse_args()

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

    if args.scan_length:
        if args.scan_length <= 0:
            logging.error('scan_length must be positive')
            return 1

        scan_length = args.scan_length
    else:
        scan_length = 1500

    try:
        ra_exp = RaExperiment(exper, band, DB(), gvlbi=args.gvlbi,
                              data_dir=data_dir, uv_fits=args.fits,
                              orbit=args.orbit)
        ra_exp.init_workdir()

        ra_exp.load(update_db=False, force_small=args.force_small,
                    scan_length=scan_length)
        ra_exp.flag_edge_chann(args.flag_chann)

        if not polar:
            if band == 'l':
                polar = 'RR'
            else:
                polar = 'LL'

        ra_exp.pima.set_polar(polar)
        ra_exp.pima.set_frq_grp(args.frequency_group)

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
                    fri = ra_exp.fringe_fitting(bandpass=not args.no_bandpass,
                                                accel=not args.no_accel,
                                                bandpass_mode=args.bpas_mode,
                                                ampl_bandpass=not args.no_ampl_bpas,
                                                bandpass_use=args.bpas_use,
                                                reference_station=args.ref_sta)
                    print('IF #{}'.format(ind+1))
                    print(fri)

                # Restore
                ra_exp.pima.update_cnt({'BEG_FRQ:': str(1),
                                        'END_FRQ:': str(if_num)})

            fri = ra_exp.fringe_fitting(bandpass=not args.no_bandpass,
                                        accel=not args.no_accel,
                                        bandpass_mode=args.bpas_mode,
                                        ampl_bandpass=not args.no_ampl_bpas,
                                        bandpass_var=args.bpas_var,
                                        bandpass_use=args.bpas_use,
                                        reference_station=args.ref_sta)
            print(fri)

            max_scan_len = fri.max_scan_length()
            logging.debug('DEBUG: max_scan_len = %s', max_scan_len)
            if args.split:
                for aver in (False, True):
                    ra_exp.split(average=aver)

                    # Copy final UV-FITS files to the system tmp directory
                    ra_exp.copy_uvfits(tempfile.gettempdir())

    except PimaError as err:
        return 1
    except RaExperimentError as err:
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
    sys.exit(main())
