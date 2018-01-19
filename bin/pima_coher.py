#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 18:22:19 2014

@author: Petr Voytsik
"""

import argparse
from collections import namedtuple
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

ObsInfo = namedtuple('ObsInfo', ['id', 'exper', 'band', 'sta1', 'sta2',
                                 'start_time'])


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


def save_txt(obs_info, dur_arr, ampl_arr, snr_arr):
    """
    Save data as a text table.

    """
    data = np.vstack((dur_arr, ampl_arr, snr_arr)).T

    file_name = '{}_{}_{:02d}_coher.{}'.format(obs_info.exper,
                                               obs_info.band,
                                               obs_info.id,
                                               'txt')

    header = '{:^4} {:^10} {:^8}'.format('solint', 'ampl', 'SNR')

    np.savetxt(file_name, data, fmt=['%6.1f', '%10.3e', '%8.1f'],
               header=header)


def plot(obs_info, dur_arr, ampl_arr, snr_arr, out_format='pdf'):
    """
    Plot
    """
    dur_arr = np.asarray(dur_arr)
    ampl_arr = np.asarray(ampl_arr)
    snr_arr = np.asarray(snr_arr)

    # Plotting
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)

    ax1.yaxis.get_major_formatter().set_powerlimits((-4, 6))

    ax1.plot(dur_arr, ampl_arr, '.')
    ax1.plot(0, 0, alpha=0)  # Dirty hack
    ax2.plot(dur_arr, snr_arr, '.')

    plot_theo_curves(ax1, ax2, dur_arr, ampl_arr, snr_arr)

    ax1.set_ylabel('Amplitude')
    ax1.grid(True)

    ax2.set_xlabel('Integration time (s)')
    ax2.set_ylabel('SNR')
    ax2.grid(True)

    ax2.text(0.95, 0.01,
             'Generated at {:%Y-%m-%d %H:%M}'.format(datetime.now()),
             color='gray', size='xx-small', ha='right',
             transform=ax2.transAxes)

    title = '{:8}({:1}) {:8} / {:8} [{:%Y-%m-%d %H:%M} UTC]'.format(
        obs_info.exper.lower(), obs_info.band.upper(), obs_info.sta1,
        obs_info.sta2, obs_info.start_time)
    fig.suptitle(title)

    plot_file_name = '{}_{}_{:02d}_coher.{}'.format(obs_info.exper,
                                                    obs_info.band,
                                                    obs_info.id,
                                                    out_format)
    if out_format == 'png':
        dpi = 300
    else:
        dpi = None

    plt.savefig(plot_file_name, bbox_inches='tight', pad_inches=0.1,
                dpi=dpi)


def proc_obs(exper, band, obs, max_dur, plot_format='pdf'):
    """
    Process observation.

    """
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
    start_time = fri[0]['start_time'] + pim.exper_info['utc_minus_tai']

    obs_info = ObsInfo(obs, exper, band, fri[0]['sta1'], fri[0]['sta2'],
                       start_time)

    dur_arr = []
    ampl_arr = []
    snr_arr = []

    for divisor in [20, 15, 12, 10, 8, 6, 5, 4, 3, 2, 1.7, 1.5, 1.3, 1.1, 1]:
        dur = round(full_duration / divisor)

#        if dur < 60:
#            continue

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
        save_txt(obs_info, dur_arr, ampl_arr, snr_arr)
        plot(obs_info, dur_arr, ampl_arr, snr_arr, out_format=plot_format)
    else:
        logging.warning('Nothing to plot...')

    return 0


def valid_obs_list(string):
    """
    Convert a string with comma-separated values to list of integers.

    """
    try:
        return [int(obs) for obs in string.split(',')]
    except ValueError:
        msg = "Not a valid comma-separated list of integers: '{0}'.".\
            format(string)
        raise argparse.ArgumentTypeError(msg)


def main():
    """Main"""
    parser = argparse.ArgumentParser()
    parser.add_argument('exper', help='experiment code')
    parser.add_argument('band', help='frequency band')
    parser.add_argument('obs_list', type=valid_obs_list,
                        help='comma-separated list of observation numbers')

    parser.add_argument('--scan-length', type=float, default=1200.,
                        help='full scan length')
    parser.add_argument('--plot-format', default='pdf',
                        help='output plot file format (PDF by default)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='be more verbose')

    args = parser.parse_args()

    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s',
                        level=logging.INFO)

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    if args.plot_format not in \
            plt.gcf().canvas.get_supported_filetypes().keys():
        logging.error('Format %s does not supported by matplotlib',
                      args.plot_format)

        return 1

    for obs_id in args.obs_list:
        proc_obs(args.exper, args.band, obs_id, args.scan_length,
                 plot_format=args.plot_format)


if __name__ == '__main__':
    main()
