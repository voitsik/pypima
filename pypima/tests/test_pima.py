#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 16:09:11 2018

@author: Petr Voytsik
"""

import pytest

from collections import namedtuple
import os.path
import shutil

from .. import Pima, Fri
from .data import SAMPLE_RA_FITS, SAMPLE_RA_CNT, SAMPLE_SCF, SAMPLE_FRI


LoadResult = namedtuple('LoadResult', ['sp_chann_num',
                                       'ap_minmax',
                                       'number_of_deselected_points',
                                       'station_list_ivs',
                                       'station_list_corr',
                                       'obs_number',
                                       ])


@pytest.fixture(scope='module')
def work_dir(tmpdir_factory):
    """
    Create and prepare a temporary working directory.
    """
    wdir = tmpdir_factory.mktemp('work_dir')
    sdir = tmpdir_factory.mktemp('pima_scr')

    for cnt_file in SAMPLE_RA_CNT.values():
        shutil.copy(cnt_file, str(wdir))

    return (str(wdir), str(sdir))


def fri_rec_equivalent(rec1, rec2):
    """
    Check two records of fri-file for equivalence.
    """
    return ((rec1['exper_code'] == rec2['exper_code']) and
            (rec1['session_code'] == rec2['session_code']) and
            (rec1['source'] == rec2['source']) and
            (rec1['sta1'] == rec2['sta1']) and
            (rec1['sat2'] == rec2['sta2']))


class TestPima:
    load_results = \
        {('raes03eo', 'l'): LoadResult(64,
                                       (1.0, 1.0),
                                       0,
                                       ['RADIO-AS', 'ZELENCHK'],
                                       ['RA', 'ZC'],
                                       1,
                                       ),
         ('raks16nq', 'c'): LoadResult(64,
                                       (0.5, 0.5),
                                       0,
                                       ['EFLSBERG', 'IRBENE16', 'RADIO-AS',
                                        'TORUN'],
                                       ['EF', 'IB', 'RA', 'TR'],
                                       14,
                                       ),
         }

    def test_environ(self):
        """
        Test environment variables.
        """
        pima_dir = os.getenv('PIMA_DIR')
        assert pima_dir is not None
        pima_exec = os.path.join(pima_dir, 'bin', 'pima')
        assert os.path.isfile(pima_exec)

    @pytest.mark.parametrize(('exper', 'band'), [('raes03eo', 'l'),
                                                 ('raks16nq', 'c')])
    def test_pima_init(self, work_dir, exper, band):
        """
        Test Pima.__init__
        """
        wdir, _ = work_dir
        pima = Pima(exper, band, work_dir=wdir)

        assert pima.exper == exper
        assert pima.band == band
        assert os.path.basename(pima.cnt_file_name) == \
            os.path.basename(SAMPLE_RA_CNT[exper, band])

    @pytest.mark.parametrize(('exper', 'band'), [('raes03eo', 'l'),
                                                 ('raks16nq', 'c')])
    def test_pima_update_cnt(self, work_dir, exper, band):
        """
        Test pima.update_cnt
        """
        wdir, sdir = work_dir
        pima = Pima(exper, band, work_dir=wdir)

        fits_file = SAMPLE_RA_FITS[exper, band]
        scf_file = SAMPLE_SCF[exper]
        pima.update_cnt({'UV_FITS:': fits_file,
                         'STAGING_DIR:': 'NO',
                         'EXPER_DIR:': str(sdir),
                         'EPHEMERIDES_FILE:': scf_file,
                         })

        assert pima.cnt_params['UV_FITS:'] == [fits_file]
        assert pima.cnt_params['STAGING_DIR:'] == 'NO'
        assert pima.cnt_params['EXPER_DIR:'] == str(sdir)
        assert pima.cnt_params['EPHEMERIDES_FILE:'] == scf_file

    @pytest.mark.parametrize(('exper', 'band'), [('raes03eo', 'l'),
                                                 ('raks16nq', 'c')])
    def test_pima_load(self, work_dir, exper, band):
        """
        Test pima.load
        """
        wdir, _ = work_dir
        pima = Pima(exper, band, work_dir=wdir)

        pima.load()
        result = self.load_results[exper, band]

        assert pima.exper_info.exper == exper
        assert pima.exper_info.band == band
        assert pima.exper_info['sp_chann_num'] == result.sp_chann_num

        assert pima.ap_minmax == result.ap_minmax
        assert pima.number_of_deselected_points == \
            result.number_of_deselected_points
        assert pima.station_list(ivs_name=False) == result.station_list_corr
        assert pima.station_list(ivs_name=True) == result.station_list_ivs
        assert pima.chan_number == result.sp_chann_num
        assert pima.obs_number == result.obs_number

    @pytest.mark.parametrize(('exper', 'band'), [('raes03eo', 'l'),
                                                 ('raks16nq', 'c')])
    def test_pima_coarse(self, work_dir, exper, band):
        """
        Test pima.coarse
        """
        wdir, _ = work_dir
        pima = Pima(exper, band, work_dir=wdir)

        fri_file = pima.coarse()

        assert os.path.isfile(fri_file)

        fri = Fri(fri_file)
        ref_fri = Fri(SAMPLE_FRI[exper, band, 'nobps'])

        rec1 = fri[0]
        rec2 = ref_fri[0]

        assert rec1['obs'] == rec2['obs']
        assert rec1['exper_code'] == rec2['exper_code']
        assert rec1['session_code'] == rec2['session_code']
        assert rec1['polar'] == rec2['polar']
        assert rec1['source'] == rec2['source']
        assert rec1['sta1'] == rec2['sta1']
        assert rec1['sta2'] == rec2['sta2']
        assert rec1['U'] == pytest.approx(rec2['U'])
        assert rec1['V'] == pytest.approx(rec2['V'])
        assert rec1['SNR'] == pytest.approx(rec2['SNR'], rel=1e-3)
        assert rec1['ampl_lsq'] == pytest.approx(rec2['ampl_lsq'], rel=1e-3)
        assert rec1['delay'] == pytest.approx(rec2['delay'], rel=1e-3)
        assert rec1['rate'] == pytest.approx(rec2['rate'], rel=1e-3)
        assert rec1['accel'] == pytest.approx(rec2['accel'], rel=1e-3)

    @pytest.mark.parametrize(('exper', 'band'), [('raes03eo', 'l'),
                                                 ('raks16nq', 'c')])
    def test_pima_bpas(self, work_dir, exper, band):
        """
        Test pima.bpas
        """
        wdir, _ = work_dir
        pima = Pima(exper, band, work_dir=wdir)

        mseg = pima.chan_number // 2

        bpas_params = {
            'BPS.MODE:': 'ACCUM',
            'BPS.NOBS_ACCUM:': '6',
            'BPS.MSEG_ACCUM:': mseg,
            'BPS.NOBS_FINE:': '12',
            'BPS.MINOBS_FINE:': '8',
            'BPS.MSEG_FINE:': mseg,
            'BPS.SNR_MIN_ACCUM:': '5.5',
            'BPS.SNR_MIN_FINE:': '5.5',
            'BPS.AMPL_REJECT:': '0.4',
            'BPS.PHAS_REJECT:': '0.2',
            'BPS.INTRP_METHOD:': 'LINEAR',
            'BPS.DEG_AMP:': '0',
            'BPS.DEG_PHS:': '1',
            'BPS.AMP_MIN:': '0.01',
            'BPS.NORML:': 'IF',
            'BPS.SEFD_USE:': 'NO',
            }

        pima.update_cnt(bpas_params)

        pima.bpas()

        bps_file = pima.cnt_params['BANDPASS_FILE:']
        assert os.path.isfile(bps_file)

    @pytest.mark.parametrize(('exper', 'band'), [('raes03eo', 'l'),
                                                 ('raks16nq', 'c')])
    def test_pima_fine(self, work_dir, exper, band):
        """
        Test pima.fine
        """
        wdir, _ = work_dir
        pima = Pima(exper, band, work_dir=wdir)

        polar = pima.cnt_params['POLAR:']
        fri_file = '{}_{}_{}.fri'.format(exper, band, polar)
        fri_file = os.path.join(wdir, fri_file)
        frr_file = '{}_{}_{}.frr'.format(exper, band, polar)
        frr_file = os.path.join(wdir, frr_file)
        pima.update_cnt({'FRINGE_FILE:': fri_file,
                         'FRIRES_FILE:': frr_file})

        fri_file = pima.fine()

        assert os.path.isfile(fri_file)

        ref_fri = Fri(SAMPLE_FRI[exper, band])
        fri = Fri(fri_file)

        rec1 = fri[0]
        rec2 = ref_fri[0]

        assert rec1['obs'] == rec2['obs']
        assert rec1['exper_code'] == rec2['exper_code']
        assert rec1['session_code'] == rec2['session_code']
        assert rec1['polar'] == rec2['polar']
        assert rec1['source'] == rec2['source']
        assert rec1['sta1'] == rec2['sta1']
        assert rec1['sta2'] == rec2['sta2']
        assert rec1['U'] == pytest.approx(rec2['U'])
        assert rec1['V'] == pytest.approx(rec2['V'])
        assert rec1['SNR'] == pytest.approx(rec2['SNR'], rel=1e-3)
        assert rec1['ampl_lsq'] == pytest.approx(rec2['ampl_lsq'], rel=1e-3)
        assert rec1['delay'] == pytest.approx(rec2['delay'], rel=1e-3)
        assert rec1['rate'] == pytest.approx(rec2['rate'], rel=1e-3)
        assert rec1['accel'] == pytest.approx(rec2['accel'], rel=1e-3)
