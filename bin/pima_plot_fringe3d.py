#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
# import os.path
import sys
import tempfile

import pypima
from pypima.pima import Pima
from pypima.fri import Fri


def sta2sta(sta):
    if sta == 'RADIO-AS':
        sta = 'RadioAstron'
    elif sta == 'GBT-VLBA':
        sta = 'GBT'
    elif sta == 'EVPTRIYA':
        sta = 'Evpatoria'
    else:
        sta = sta.title()

    return sta


def sou2sou(sou):
    if sou == '3C274':
        sou = 'M87'

    return sou


def plot(file_name, fri, format, title=False):
    axis1_num_point = 0
    axis2_num_point = 0
    axis3_num_point = 0
    plot_title = ''
    subtitle = ''
    exper, band = fri['session_code'].split('_')
    obs = fri['obs']

    with open(file_name, 'r') as fil:
        for line in fil:
            key, sep, val = line.partition(':')
            if key == 'PLOT_TITLE':
                plot_title = val.strip()
            elif key == 'SUBTITLE':
                subtitle = val.strip()
            elif key == 'AXIS1_NUM_POI':
                axis1_num_point = int(val.strip())
            elif key == 'AXIS2_NUM_POI':
                axis2_num_point = int(val.strip())
            elif key == 'AXIS3_NUM_POI':
                axis3_num_point = int(val.strip())
            elif key == 'POINT':
                break

        i, j, x, y, z = val.replace('D', 'e').split()
        X = np.zeros([axis1_num_point, axis2_num_point])
        Y = np.zeros([axis1_num_point, axis2_num_point])
        Z = np.zeros([axis1_num_point, axis2_num_point])
        i = int(i)-1
        j = int(j)-1
        X[i, j] = float(x)
        Y[i, j] = float(y)
        Z[i, j] = float(z)

        for line in fil:
            key, sep, val = line.partition(':')
            if key == 'POINT':
                i, j, x, y, z = val.replace('D', 'e').split()
                i = int(i)-1
                j = int(j)-1
                X[i, j] = float(x)
                Y[i, j] = float(y)
                Z[i, j] = float(z)

    mean = np.mean(Z)
    std = np.std(Z)
    mean_corr = np.mean(Z[Z < mean + 2*std])

    print("Mean = {:.3e}".format(mean))
    print("Std = {:.3e}".format(std))
    print("Peak = {:.3e}".format(Z.max()))
    print("Peak/mean = {:.1f}" .format(Z.max() / mean_corr))

    plt.style.use('seaborn-paper')
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    freq = fri['ref_freq']

    ax.plot_surface(X*1e9, freq*Y*1e3, Z / mean_corr, rstride=1, cstride=1,
                    cmap=cm.jet, linewidth=0)
    # ax.set_title(subtitle, fontsize='small')
    ax.set_ylabel('Fringe rate (mHz)')
    ax.set_xlabel('Delay (ns)')
    ax.set_zlabel('Signal-to-Noise Ratio')
    # ax.text2D(0.05, 0.85, "Date: 14 Mar 2012", transform=ax.transAxes)
    # ax.text2D(0.05, 0.80, "Time of integration: 65s", transform=ax.transAxes)
    # ax.grid(b=False)
    # ax.w_xaxis.set_pane_color((0, 0, 0, 0))
    # ax.set_zticklabels((), alpha=0)
    # ax.view_init(elev=90, azim=0)
    source = sou2sou(fri['source'])

    if title:
        date = fri['start_time'].strftime('%d.%m.%Y')
        band2wl = {'k': 1.35, 'c': 6, 'l': 18}
        title = '{} ({})\n'.format(exper, date)
        title += '{}, {} cm, {}-{}'.format(source, band2wl[band],
                                           sta2sta(fri['sta1']),
                                           sta2sta(fri['sta2']))
        title += ' ({:.1f} Earth Diameters)'.format(fri['uv_rad_ed'])
        # plt.title(title)
        # fig.suptitle(suptitle)
        ax.set_title(title)

    file_name = '{}_{}_{}_{:02d}_fringe3D'.format(source, exper, band, obs)
    file_name = '{}.{}'.format(file_name, format)

    if format == 'eps':
        ax.w_yaxis.set_pane_color((1, 1, 1, 0))

    if format == 'png':
        dpi = 300
    else:
        dpi = None
    plt.savefig(file_name, format=format, bbox_inches='tight', pad_inches=0.1,
                dpi=dpi)


def main(args):
    """Main"""
    exper = args.exper.lower()
    band = args.band.lower()
    obs = args.obs

    try:
        pim = Pima(exper, band)
    except pypima.pima.Error as err:
        print('PIMA Error: ', err, file=sys.stderr)
        return 1

    if obs <= 0 or obs > pim.obs_number:
        print('Incorrect observation number {} must be in range [{} {}]'.
              format(obs, 1, pim.obs_number))
        return 1

    params = ['FRIB.OBS:', str(obs), 'FRIB.2D_FRINGE_PLOT:', 'TXT']

    if args.delay_window:
        delay_window = args.delay_window * 1e-6
    else:
        delay_window = 500e-9

    if args.rate_window:
        rate_window = args.rate_window
    elif band == 'k':
        rate_window = 2e-12
    elif band == 'c':
        rate_window = 4e-12
    else:
        rate_window = 8e-12

    delay_str = '{:E}'.format(delay_window)
    rate_str = '{:E}'.format(rate_window)

    params.extend(['FRIB.PLOT_DELAY_WINDOW_WIDTH:', delay_str,
                   'FRIB.PLOT_RATE_WINDOW_WIDTH:', rate_str])

    with tempfile.NamedTemporaryFile(suffix='.fri') as tmp_fri:
        params.extend(['FRINGE_FILE:', tmp_fri.name])
        try:
            fri = Fri(pim.fine(params))
        except pypima.pima.Error as err:
            print('PIMA Error: ', err, file=sys.stderr)
            return 1

    time_code = fri[0]['time_code']
    sta1 = fri[0]['sta1'].lower()
    sta1 = sta1 + '_' * (8 - len(sta1))
    sta2 = fri[0]['sta2'].lower()
    # sta2 = sta2 + '_' * (8 - len(sta2))
    path = '{}/{}_fpl/fr2d_{}__{}_{}_{}.txt'.format(
        pim.cnt_params['EXPER_DIR:'],
        pim.cnt_params['SESS_CODE:'],
        time_code, band, sta1, sta2)

    # print('DEBUG: path =', path)

    try:
        plot(path, fri[0], args.format, title=args.title)
    except OSError as err:
        print('OSError: ', err, file=sys.stderr)
        return 1

    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('exper', help='experiment code')
    parser.add_argument('band', help='frequency band code')
    parser.add_argument('obs', type=int, help='observation number')

    parser.add_argument('--title', action='store_true',
                        help='add title to plot')
    parser.add_argument('--format', default='pdf',
                        help='plot file format (pdf, ps, png, ...)')

    parser.add_argument('-d', '--delay-window', type=float,
                        help='delay window width in microseconds')
    parser.add_argument('-r', '--rate-window', type=float,
                        help='rate window width')

    sys.exit(main(parser.parse_args()))
