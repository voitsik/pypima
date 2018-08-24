#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 16:09:11 2018

@author: Petr Voytsik
"""

import pytest

import os.path
import shutil

from .. import Pima
from ..data import SAMPLE_RA_FITS, SAMPLE_RA_CNT, SAMPLE_SCF


@pytest.fixture(scope='module')
def work_dir(tmpdir_factory):
    """
    Create and prepare a temporary working directory.
    """
    wd = tmpdir_factory.mktemp('raes03eo_l')

    shutil.copy(SAMPLE_RA_CNT, str(wd))

    return str(wd)


class TestPima:
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
        pima = Pima(exper, band, work_dir=work_dir)

        assert pima.exper == exper
        assert pima.band == band
        assert os.path.basename(pima.cnt_file_name) == \
            os.path.basename(SAMPLE_RA_CNT)

    @pytest.mark.parametrize(('exper', 'band'), [('raes03eo', 'l')])
    def test_pima_load(self, tmpdir, work_dir, exper, band):
        pima = Pima(exper, band, work_dir=work_dir)

        pima.update_cnt({'UV_FITS:': SAMPLE_RA_FITS,
                         'STAGING_DIR:': 'NO',
                         'EXPER_DIR:': str(tmpdir),
                         'EPHEMERIDES_FILE:': SAMPLE_SCF,
                        })

        assert pima.cnt_params['UV_FITS:'] == [SAMPLE_RA_FITS]

        pima.load()
        assert pima.exper_info.exper == exper
        assert pima.exper_info.band == band
        assert pima.exper_info['sp_chann_num'] == 64

        assert pima.ap_minmax() == (1.0, 1.0)
        assert pima.number_of_deselected_points() == 0
        assert pima.station_list(ivs_name=False) == ['RA', 'ZC']
        assert pima.station_list(ivs_name=True) == ['RADIO-AS', 'ZELENCHK']
        assert pima.chan_number() == 64
        assert pima.obs_number() == 1
