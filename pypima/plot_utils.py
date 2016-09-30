# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 14:01:07 2016

@author: Petr Voytsik
"""

from datetime import timedelta
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np
import os
import pypima.pima


def generate_autospectra(pima, out_dir):
    """
    Plot autospectrum for each station for each scan.

    Parameters
    ----------
    pima : ``Pima`` object
        Instance of ``Pima`` object
    out_dir : str
        Output directory

    """
    # Sometimes PIMA crashes on `acta` task
    try:
        file_list = pima.acta()
    except pypima.pima.Error:
        # Remove core dump file.
        if os.path.isfile('core'):
            os.remove('core')

        return

    out_dir = os.path.join(out_dir, '{}_{}'.format(pima.exper, pima.band))
    os.makedirs(out_dir, exist_ok=True)

    polar = pima.cnt_params['POLAR:']

    central_freq = 1e-6 * pima.frequencies()[-1]['freq']
    utc_tai = timedelta(seconds=pima.exper_info['utc_minus_tai'])
    out_format = 'pdf'
    fig = Figure()
    FigureCanvas(fig)
    ax = fig.add_subplot(111)

    for txt_path in file_list:
        acta_file = pypima.pima.ActaFile(txt_path)

        # Convert TAI to UTC
        date = acta_file.header['start_date'] + utc_tai

        date_str = date.strftime('%Y-%m-%d %H:%M:%S')
        ax.cla()
        ax.set_title('{} - {}({}) - {}'.format(acta_file.header['station'],
                                               pima.exper,
                                               pima.band.upper(),
                                               date_str))
        ax.set_xlabel('Frequency (MHz)')
        ax.set_ylabel('Amplitude')
        ax.grid(True)
        ax.plot(np.asarray(acta_file.freq) - central_freq,
                acta_file.ampl, marker='o', ms=2)

        date_str = date.strftime('%Y%m%dT%H%M')
        out_file = \
            'AUTOSPEC_{}_{}_{}_{}_{}.{}'.format(date_str,
                                                pima.exper,
                                                pima.band.upper(),
                                                polar,
                                                acta_file.header['station'],
                                                out_format)

        out_file = os.path.join(out_dir, out_file)
        fig.savefig(out_file, format=out_format)
