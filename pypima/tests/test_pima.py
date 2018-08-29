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

from .. import Pima
from .data import SAMPLE_RA_FITS, SAMPLE_RA_CNT, SAMPLE_SCF


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


class TestPima:
    load_results = \
        {('raes03eo', 'l'): LoadResult(64,
                                       (1.0, 1.0),
                                       0,
                                       ['RADIO-AS', 'ZELENCHK'],
                                       ['RA', 'ZC'],
                                       1,
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

    @pytest.mark.parametrize(('exper', 'band'), [('raes03eo', 'l')])
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

    @pytest.mark.parametrize(('exper', 'band'), [('raes03eo', 'l')])
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

    @pytest.mark.parametrize(('exper', 'band'), [('raes03eo', 'l')])
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

    @pytest.mark.parametrize(('exper', 'band'), [('raes03eo', 'l')])
    def test_pima_coarse(self, work_dir, exper, band):
        """
        Test pima.coarse
        """
        wdir, _ = work_dir
        pima = Pima(exper, band, work_dir=wdir)

        fri = pima.coarse()

        assert os.path.isfile(fri)
