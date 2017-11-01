#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 18:22:19 2014

@author: Petr Voytsik
"""

import argparse
from datetime import datetime
import logging
import math
import os.path
from tempfile import NamedTemporaryFile

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from pypima.pima import Pima
from pypima.pima import Error as PIMAError
from pypima.fri import Fri


def plot_theo_curves(ax_ampl, ax_snr, dur_arr, ampl_arr, snr_arr,
                     max_coher_dur=400):
    """
    Plot theoretical curves for amplitude and SNR.

    """
    # Mask to select values with int time less then coher time
    coher_mask = dur_arr < max_coher_dur

    # At least 2 points to fit
    if np.count_nonzero(coher_mask) <= 1:
        return

    def sqrt_func(time, ampl):
        """
        Return A * sqrt(T)
        """
        return ampl * np.sqrt(time)

    fit_params, _ = curve_fit(sqrt_func,
                              dur_arr[coher_mask],
                              snr_arr[coher_mask])

    dur_theo_arr = np.linspace(0, dur_arr.max(), num=200)
    snr_theo_arr = sqrt_func(dur_theo_arr, *fit_params)
    ampl_theo_arr = np.full_like(dur_theo_arr, ampl_arr[coher_mask].mean())

    ax_ampl.plot(dur_theo_arr, ampl_theo_arr, 'r-')
    ax_snr.plot(dur_theo_arr, snr_theo_arr, 'r-')


def plot(obs_id, dur_arr, ampl_arr, snr_arr, exper, band, sta1, sta2,
         start_time):
    """
    Plot
    """
    dur_arr = np.asarray(dur_arr)
    ampl_arr = np.asarray(ampl_arr)
    snr_arr = np.asarray(snr_arr)

    if np.mean(ampl_arr) < 1e-3:
        ylabel = 'Amplitude (micro values)'
        ampl_arr *= 1e6
    else:
        ylabel = 'Amplitude'

    # Plotting
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    ax1.set_ymargin(0.2)
    ax2.set_ymargin(0.1)
    ax1.plot(dur_arr, ampl_arr, 'o')
    ax2.plot(dur_arr, snr_arr, 'o')

    plot_theo_curves(ax1, ax2, dur_arr, ampl_arr, snr_arr)

    ax1.set_ylabel(ylabel)
    title = '{}({}) obs #{}: {} / {}'.format(exper.lower(), band.upper(),
                                             obs_id, sta1, sta2)
    ax1.set_title(title)
    ax1.set_ylim(ymin=0)
    ax1.grid(True)

    ax2.set_xlabel('Integration time (s)')
    ax2.set_ylabel('SNR')
    ax2.grid(True)

    ax2.text(0.95, 0.03,
             'Generated at {:%Y-%m-%d %H:%M}'.format(datetime.now()),
             color='gray', size='xx-small', ha='right',
             transform=ax2.transAxes)

    plot_file_name = '{}_{}_{:02d}_coher.pdf'.format(exper, band, obs_id)
    plt.savefig(plot_file_name, format='pdf')


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
        logging.error(
            'Incorrect observation number %s, must be in range [%s %s]',
            obs, 1, pim.obs_number())
        return 1

    # Prepare temporary fri and frr files
    with NamedTemporaryFile(suffix='.fri', delete=False) as tmp:
        tmp_fri = tmp.name
    with NamedTemporaryFile(suffix='.frr', delete=False) as tmp:
        tmp_frr = tmp.name

    try:
        fri_file = pim.fine(['FRIB.OBS:', str(obs),
                             'SCAN_LEN_SKIP:', '0.0',
                             'SCAN_LEN_USED:', str(max_dur),
                             'FRINGE_FILE:', tmp_fri,
                             'FRIRES_FILE:', tmp_frr])
    except PIMAError as err:
        logging.error('PIMA Error: %s', err)
        return 1

    fri = Fri(fri_file)

    if not fri:
        logging.error('PIMA fri-file is empty')
        return 1

    if fri[0]['SNR'] < 7.0:
        logging.warning('SNR of obs #%s is too low, skip it.', obs)
        return 1

    full_duration = fri[0]['duration']
    sta1 = fri[0]['sta1']
    sta2 = fri[0]['sta2']
    start_time = fri[0]['start_time']

    dur_arr = []
    ampl_arr = []
    snr_arr = []

    for divisor in [6, 5, 4, 3, 2, 1.7, 1.5, 1.3, 1.1, 1]:
        dur = round(full_duration / divisor)
        logging.debug('Set SCAN_LEN_USED %s', dur)

        for j in range(math.ceil(divisor)):
            if divisor < 2 and j == 1:
                skip = full_duration - dur
            else:
                skip = j * dur

            logging.debug('Set SCAN_LEN_SKIP to %s', skip)
            fri_file = pim.fine(['FRIB.OBS:', str(obs),
                                 'SCAN_LEN_SKIP:', str(skip),
                                 'SCAN_LEN_USED:', str(dur),
                                 'FRINGE_FILE:', tmp_fri,
                                 'FRIRES_FILE:', tmp_frr])
            fri = Fri(fri_file)
            if not fri:
                logging.warning('Skip empty fri-file')
                continue

            logging.debug('SNR = %s', fri[0]['SNR'])

            if fri[0]['SNR'] < 5.7:
                continue

            dur_arr.append(fri[0]['duration'])
            ampl_arr.append(fri[0]['ampl_lsq'])
            snr_arr.append(fri[0]['SNR'])

            print('{} {} {}'.format(dur_arr[-1], ampl_arr[-1], snr_arr[-1]))

    # Delete temporary files
    if os.path.isfile(tmp_fri):
        os.remove(tmp_fri)
    if os.path.isfile(tmp_frr):
        os.remove(tmp_frr)

    if dur_arr:
        plot(obs, dur_arr, ampl_arr, snr_arr, exper, band, sta1, sta2,
             start_time)
    else:
        logging.warning('Nothing to plot...')

    return 0


def main():
    """Main"""
    parser = argparse.ArgumentParser()
    parser.add_argument('exper', help='experiment code')
    parser.add_argument('band', help='frequency band')
    parser.add_argument('obs_list',
                        help='coma separated list of observation numbers')

    parser.add_argument('--scan-length', type=float, default=1200.,
                        help='full scan length')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='be more verbose')

    args = parser.parse_args()

    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s',
                        level=logging.INFO)

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    for obs_id in args.obs_list.split(','):
        proc_obs(args.exper, args.band, int(obs_id), args.scan_length)


if __name__ == '__main__':
    main()
