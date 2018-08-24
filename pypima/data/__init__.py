# -*- coding: utf-8 -*-
"""Sample files"""

# Use private names to avoid inclusion in the sphinx documentation.
from os import path as _path


def _full_path(name, dirname=_path.dirname(_path.abspath(__file__))):
    return _path.join(dirname, name)


SAMPLE_RA_FITS = \
    _full_path('RADIOASTRON_RAES03EO_L_20121016T100000_ASC_V3.idifits')
"""RadioAstron FITS-IDI sample.
"""
