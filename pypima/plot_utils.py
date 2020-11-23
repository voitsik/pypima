# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 14:01:07 2016

@author: Petr Voytsik
"""

import logging
import os

import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

# import pypima.pima


def plot_autospectra(acta_file_list, out_dir_base):
    """
    Plot each autocorrelation spectrum in `acta_file_list`.

    Parameters
    ----------
    acta_file_list : list of ``ActaFile`` objects
        List of autospectra in ``ActaFile`` format.
    out_dir_base : str
        Base output directory.

    """
    out_format = "pdf"
    fig = Figure()
    FigureCanvas(fig)
    ax = fig.add_subplot(111)

    for acta_file in acta_file_list:
        exper, band = acta_file.header["experiment"].split("_")
        out_dir = os.path.join(out_dir_base, "{}_{}".format(exper, band))
        os.makedirs(out_dir, exist_ok=True)

        date = acta_file.header["start_date"]
        date_str = date.strftime("%Y-%m-%d %H:%M:%S")

        ax.cla()
        ax.set_ymargin(0.05)
        ax.xaxis.set_ticks(np.arange(-16, 16 + 1, 4))

        ax.set_title(
            "{} - {}({}) - {}".format(
                acta_file.header["station"], exper, band.upper(), date_str
            )
        )
        ax.set_xlabel("Frequency (MHz)")
        ax.set_ylabel("Amplitude")
        ax.grid(True)
        central_freq = np.mean(acta_file.freq)
        logging.debug("plot_autospectra: central_freq = %s", central_freq)
        freq = np.asarray(acta_file.freq) - central_freq
        ax.plot(freq, acta_file.ampl, marker="o", ms=2)
        ax.set_xlim(-16, 16)

        date_str = date.strftime("%Y%m%dT%H%M")
        polar = acta_file.header["polar"]
        out_file = "AUTOSPEC_{}_{}_{}_{}_{}.{}".format(
            date_str,
            exper,
            band.upper(),
            polar,
            acta_file.header["station"],
            out_format,
        )

        out_file = os.path.join(out_dir, out_file)
        fig.savefig(out_file, format=out_format)


# def plot_ampl_phas(fri, out_dir):
#    """
#    """
#    if not fri:
#        return
#
#    # Check if Session code in RA AGN survey format
#    if '_' not in fri[0]['session_code']:
#        return
#
#    exper, band = fri[0]['session_code'].split('_')
#    out_dir = os.path.join(out_dir, '{}_{}'.format(exper, band))
#    os.makedirs(out_dir, exist_ok=True)
#
#    out_format = 'pdf'
#    fig = Figure()
#    FigureCanvas(fig)
#    ax_phas = fig.add_subplot(211)
#    ax_ampl = fig.add_subplot(212, sharex=ax_phas)
#
#    for rec in fri:
#        polar = rec['polar']
#        time_code = rec['time_code']
#        sta1 = rec['sta1']
#        sta2 = rec['sta2']
#
#        file_name_base = '{:_<8}_{}_{:_<8}_{:_<8}'.format(time_code, band,
#                                                          sta1.lower(),
#                                                          sta2.lower())
