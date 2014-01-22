#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 15:21:54 2014

@author: Petr Voytsik
"""

import os
#import sys
from shutil import move
import pima as pypima
from optparse import OptionParser


def source_names(name, source_names_file):
    """Return 3 names (IVS-name, J2000-name, B1950-name) for given source"""
    with open(source_names_file) as fil:
        for line in fil:
            if len(line) < 30 or line.startswith('#'):
                continue
            toks = line.split()
            if name in toks[:3]:
                return toks[:3]

    print('Error: could not find source {} in {}'.format(name,
          source_names_file))


def pima_split(pima, aver=0.0, nocl=False, source=None, polar=None,
               _pima_opts=None):
    """PIMA SPLIT: do not averaging if aver is None"""
    pima_opts = []
    if _pima_opts is not None:
        pima_opts.extend(_pima_opts)

    ap_min, ap_max = pima.ap_minmax()
    if ap_min != ap_max:
        print('Warning: different accummulation periods during \
        experiment! Use min.')
    mseg = 1
    if aver > 0:
        mseg = int(aver / ap_min)
    else:
        aver = ap_min

    if polar:
        if polar != pima.cnt_params['POLAR:']:
            print('Info: Polarization differs from previous one: recalibrate \
data')
            pima.update_cnt({'POLAR:': polar, 'SPLT.POLAR:': polar})
            try:
#                pima.load()
                pima.coarse()
                pima.bpas()
                pima.fine()
#                antab_file = '{}/{}{}.antab'.format(pima.work_dir, pima.exper,
#                             pima.band)
#                pima.load_gains(antab_file)
#                pima.load_tsys(antab_file)
            except Exception as ex:
                print('PIMA Error: {}'.format(ex))
                return

    if source:
        pima_opts.extend(['SPLT.SOU_NAME:', source])
        so_names_file = pima.get_cnt_params(['SOU_NAMES:'])['SOU_NAMES:']
        ivs_name, j2000_name, b1950_name = source_names(source, so_names_file)

    if nocl:
        pima_opts.extend(['PIMAVAR_SPLT_SUB_SRT:', 'YES'])

    pima_opts.extend(['SPLT.WEIGHT_TYPE:', 'RMS'])
    pima.split(tim_mseg=mseg, params=pima_opts)

    if source:
        params = pima.get_cnt_params(['EXPER_DIR:', 'SESS_CODE:', 'BAND:',
                                      'POLAR:'])
        #Copy FITS file from pima_scr to final distination
        out_dir = '/home/voitsik/RadioAstron/VLBI'
        pima_fits_path = '{}/{}_uvs/{}_{}_uva.fits'.format(
                         params['EXPER_DIR:'],
                         params['SESS_CODE:'], j2000_name, params['BAND:'])
        if not os.path.isfile(pima_fits_path):
            print("Error: can't find pima output FITS {}".format(
                pima_fits_path))
            return

        final_fits_name = '{}_{}_{}_{}_{}s_uva.fits'.format(b1950_name,
                          pima.exper, pima.band.upper(), params['POLAR:'],
                          int(aver))
        if nocl:
            final_fits_name = final_fits_name.replace('_uva', '_noclosure_uva')
        final_fits_dir = '{}/{}'.format(out_dir, b1950_name)
        if not os.path.isdir(final_fits_dir):
            os.mkdir(final_fits_dir)
        final_fits_path = '{}/{}'.format(final_fits_dir, final_fits_name)
        print('Info: move output FITS file to {}'.format(final_fits_path))
        move(pima_fits_path, final_fits_path)
        pypima.fits_to_txt(final_fits_path)


def main():
    """Main: cmd line arguments parsing"""
    usage = "usage: %prog [-n] <exper> <band> [aver[,aver]]"
    parser = OptionParser(usage=usage, version="%prog 0.1")

    parser.set_defaults(noclosure=False)
    parser.add_option('-n', '--no-closure', action="store_true",
                      dest='noclosure',
                      help='Refer aposteriori VLBI model to the subarray \
                      reference time')

    parser.add_option("-o", "--pima_opts", action="store",
                      dest="pima_opts",
                      metavar='"KEY: VAL"',
                      default="",
                      help="Addtional pima options")

    parser.add_option("-s", "--source", action="store",
                      dest="source",
                      metavar='source_name',
                      default='',
                      help='Source name (by default ALL)')

    parser.add_option("-p", "--polar", action="store",
                      dest="polar",
                      metavar='POL',
                      default='',
                      help='Polarization')

    opts, args = parser.parse_args()

    if len(args) < 2:
        parser.error("incorrect number of arguments")

    exper = args[0].lower()
    band = args[1].lower()
    noclosure = opts.noclosure
    pima_opts = opts.pima_opts.split()
    source = opts.source
    polar = opts.polar

    if len(pima_opts) % 2:
        print('Error: The number of pima arguments should be even, but you \
               specified {}'.format(len(pima_opts)))

    exp_dir = os.getenv('pima_exp_dir')
    work_dir = '{}/{}'.format(exp_dir, exper)
    os.chdir(work_dir)
    pim = pypima.Pima(exper, band, work_dir)

    # Averaging in seconds. If None do not aver
    aver_list = [0]
    if len(args) == 3:
        aver_list = args[2].split(',')

    for aver in aver_list:
        pima_split(pim, float(aver), noclosure, source, polar, pima_opts)

if __name__ == "__main__":
    main()
