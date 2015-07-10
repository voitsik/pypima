#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 18:22:19 2014

@author: Petr Voytsik
"""

import argparse
import math
import matplotlib.pyplot as plt
import numpy as np
import os.path
from scipy.optimize import curve_fit
import sys
import tempfile

PATH = os.path.normpath(os.path.join(os.path.dirname(sys.argv[0]), '..'))
sys.path.insert(0, PATH)
from pypima.pima import Pima
from pypima.fri import Fri


def func(x, A):
    """
    SNR = A * sqrt(time)

    """
    return A * np.sqrt(x)


def main(exper, band, obs):
    """Main"""
    pim = Pima(exper, band)

    if obs <= 0 or obs > pim.obs_number():
        print('Incorrect observation number {}, must be in range [{} {}]'.
              format(obs, 1, pim.obs_number()))
        return 1

    sta1 = sta2 = 'TEST'

    durs = []
    amps = []
    snrs = []

    if os.path.isfile('test1'):
        durs, amps, snrs = np.loadtxt('test1', unpack=True)
    else:
        orig_fri = pim.cnt_params['FRINGE_FILE:']
        orig_frr = pim.cnt_params['FRIRES_FILE:']
        tmp_fri = tempfile.NamedTemporaryFile(suffix='.fri')
        tmp_frr = tempfile.NamedTemporaryFile(suffix='.frr')
        pim.update_cnt({'FRINGE_FILE:': tmp_fri.name,
                        'FRIRES_FILE:': tmp_frr.name})

        fri_file = pim.fine(['FRIB.OBS:', str(obs),
                             'SCAN_LEN_SKIP:', '0.0',
                             'SCAN_LEN_USED:', '1200.0'])
        fri = Fri(fri_file)
        full_duration = fri[0]['duration']
        sta1 = fri[0]['sta1']
        sta2 = fri[0]['sta2']

        for delim in [6, 5, 4, 3, 2, 1.7, 1.5, 1.3, 1.1, 1]:
            dur = full_duration / (delim)
            for j in range(math.ceil(delim)):
                skip = j * dur
                fri_file = pim.fine(['FRIB.OBS:', str(obs),
                                     'SCAN_LEN_SKIP:', str(skip),
                                     'SCAN_LEN_USED:', str(dur)])
                fri = Fri(fri_file)
                if len(fri) == 0:
                    continue

                if fri[0]['SNR'] < 5.7:
                    continue

                durs.append(fri[0]['duration'])
                amps.append(fri[0]['ampl_lsq'])
                snrs.append(fri[0]['SNR'])

                print('{} {} {}'.format(durs[-1], amps[-1], snrs[-1]))

        # Restore original fri and frr files
        pim.update_cnt({'FRINGE_FILE:': orig_fri,
                        'FRIRES_FILE:': orig_frr})
        tmp_fri.close()
        tmp_frr.close()

    durs = np.array(durs)
    amps = np.array(amps) * 1e6  # Amplitude in micro values
    snrs = np.array(snrs)

    p, pcov = curve_fit(func, durs[durs < 400], snrs[durs < 400])
    # err = np.sqrt(np.diag(pcov))

    durs0 = np.linspace(0, durs.max())
    sqrt_snr = func(durs0, *p)
    const_amp = np.ones(durs0.shape) * amps[durs < 400].mean()

    # Plotting
    f, (ax1, ax2) = plt.subplots(2, sharex=True)

    ax1.set_ymargin(0.2)
    ax1.plot(durs, amps, 'o', durs0, const_amp, 'r-')
    ax1.set_ylabel('Amplitude (micro values)')
    title = '{}({}): {} / {}'.format(exper.lower(), band.upper(), sta1, sta2)
    ax1.set_title(title)
#    ax1.set_ylim(ymin=0)

    ax2.plot(durs, snrs, 'o', durs0, sqrt_snr, 'r-')
    ax2.set_xlabel('Integration time (s)')
    ax2.set_ylabel('SNR')

    pic_file_name = '{}_{}_{:02d}_coher.pdf'.format(exper, band, obs)
    plt.savefig(pic_file_name, format='pdf')
#    plt.show()
    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('exper', help='experiment code')
    parser.add_argument('band', help='frequency band')
    parser.add_argument('obs_list',
                        help='coma separated list of observation numbers')

    args = parser.parse_args()

    for obs in args.obs_list.split(','):
        main(args.exper, args.band, int(obs))
