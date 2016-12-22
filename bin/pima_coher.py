#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 18:22:19 2014

@author: Petr Voytsik
"""

import argparse
import logging
import math
import matplotlib.pyplot as plt
import numpy as np
import os.path
from scipy.optimize import curve_fit
import sys
from tempfile import NamedTemporaryFile

PATH = os.path.normpath(os.path.join(os.path.dirname(sys.argv[0]), '..'))
sys.path.insert(0, PATH)
from pypima.pima import Pima
from pypima.pima import Error as PIMAError
from pypima.fri import Fri


def plot(obs, durs, amps, snrs, exper, band, sta1, sta2):
    """
    Plot
    """
    if np.mean(amps) < 1e-3:
        ylabel = 'Amplitude (micro values)'
        amps *= 1e6
    else:
        ylabel = 'Amplitude'

    def sqrt_func(dur, amp):
        return amp * np.sqrt(dur)

    # Plotting
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    ax1.set_ymargin(0.2)
    ax2.set_ymargin(0.1)
    ax1.plot(durs, amps, 'o')
    ax2.plot(durs, snrs, 'o')

    durs_400 = durs[durs < 400]

    if len(durs_400) > 1:
        snrs_400 = snrs[durs < 400]
        params, pcov = curve_fit(sqrt_func, durs_400, snrs_400)
        # err = np.sqrt(np.diag(pcov))

        durs0 = np.linspace(0, durs.max())
        sqrt_snr = sqrt_func(durs0, *params)
        const_amp = np.ones(durs0.shape) * amps[durs < 400].mean()

        ax1.plot(durs0, const_amp, 'r-')
        ax2.plot(durs0, sqrt_snr, 'r-')

    ax1.set_ylabel(ylabel)
    title = '{}({}) obs #{}: {} / {}'.format(exper.lower(), band.upper(),
                                             obs, sta1, sta2)
    ax1.set_title(title)
#    ax1.set_ylim(ymin=0)

    ax2.set_xlabel('Integration time (s)')
    ax2.set_ylabel('SNR')

    pic_file_name = '{}_{}_{:02d}_coher.pdf'.format(exper, band, obs)
    plt.savefig(pic_file_name, format='pdf')


def proc_obs(exper, band, obs, max_dur):
    """Main"""
    try:
        pim = Pima(exper, band)
    except PIMAError as err:
        logging.error('PIMA Error: %s', err)
        return 1
    except OSError as err:
        logging.error('OSError: %s', err)
        return 1

    if obs <= 0 or obs > pim.obs_number():
        logging.error('Incorrect observation number %s, must be in range [%s %s]',
                      obs, 1, pim.obs_number())
        return 1

    # Prepare temporary fri and frr files
    with NamedTemporaryFile(suffix='.fri', delete=False) as tmp:
        tmp_fri = tmp.name
    with NamedTemporaryFile(suffix='.frr', delete=False) as tmp:
        tmp_frr = tmp.name

    fri_file = pim.fine(['FRIB.OBS:', str(obs),
                         'SCAN_LEN_SKIP:', '0.0',
                         'SCAN_LEN_USED:', str(max_dur),
                         'FRINGE_FILE:', tmp_fri,
                         'FRIRES_FILE:', tmp_frr])

    fri = Fri(fri_file)
    if fri[0]['SNR'] < 7.0:
        logging.warning('SNR of obs #%s is too low, skip it.', obs)
        return 1

    full_duration = fri[0]['duration']
    sta1 = fri[0]['sta1']
    sta2 = fri[0]['sta2']

    durs = []
    amps = []
    snrs = []

    for delim in [6, 5, 4, 3, 2, 1.7, 1.5, 1.3, 1.1, 1]:
        dur = full_duration / (delim)
        logging.info('Set scan length to %s', dur)

        for j in range(math.ceil(delim)):
            skip = j * dur
            fri_file = pim.fine(['FRIB.OBS:', str(obs),
                                 'SCAN_LEN_SKIP:', str(skip),
                                 'SCAN_LEN_USED:', str(dur),
                                 'FRINGE_FILE:', tmp_fri,
                                 'FRIRES_FILE:', tmp_frr])
            fri = Fri(fri_file)
            if len(fri) == 0:
                continue

            if fri[0]['SNR'] < 5.7:
                continue

            durs.append(fri[0]['duration'])
            amps.append(fri[0]['ampl_lsq'])
            snrs.append(fri[0]['SNR'])

            print('{} {} {}'.format(durs[-1], amps[-1], snrs[-1]))

    # Delete temporary files
    if os.path.isfile(tmp_fri):
        os.remove(tmp_fri)
    if os.path.isfile(tmp_frr):
        os.remove(tmp_frr)

    if durs:
        plot(obs, np.array(durs), np.array(amps), np.array(snrs),
             exper, band, sta1, sta2)
    else:
        logging.warning('Nothing to plot...')

    return 0


def main(args):
    """Main"""
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s',
                        level=logging.INFO)

    for obs_id in args.obs_list.split(','):
        proc_obs(args.exper, args.band, int(obs_id), args.scan_length)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('exper', help='experiment code')
    parser.add_argument('band', help='frequency band')
    parser.add_argument('obs_list',
                        help='coma separated list of observation numbers')

    parser.add_argument('--scan-length', type=float, default=1200.,
                        help='full scan length')

    main(parser.parse_args())
