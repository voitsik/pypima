# -*- coding: utf-8 -*-
"""Sample files"""

# Use private names to avoid inclusion in the sphinx documentation.
from os import path as _path


def _full_path(name, dirname=_path.dirname(_path.abspath(__file__))):
    return _path.join(dirname, name)


SAMPLE_RA_FITS = \
    {('raes03eo', 'l'):
        _full_path('RADIOASTRON_RAES03EO_L_20121016T100000_ASC_V3.idifits'),
     }
"""RadioAstron FITS-IDI sample.
"""

SAMPLE_RA_CNT = {('raes03eo', 'l'): _full_path('raes03eo_l_pima.cnt'),
                 }
"""PIMA control file sample.
"""

SAMPLE_SCF = {'raes03eo': _full_path('RA121015_2100_v02.scf'),
              }
"""Reconstructed orbit file sample.
"""
