#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 17:50:20 2013

@author: Petr Voytsik
"""

from __future__ import print_function
import sys
import os.path
path = os.path.normpath(os.path.join(os.path.dirname(sys.argv[0]), '..'))
sys.path.insert(0, path)
import pypima
from pypima.raexperiment import RaExperiment

if len(sys.argv) != 3:
    print('Usage: {} <exper> <band>'.format(sys.argv[0]), file=sys.stderr)
    sys.exit(2)

try:
    p = RaExperiment(sys.argv[1], sys.argv[2])
    p.load()

    for polar in ['RR', 'RL', 'LR', 'LL']:
        p.pima.set_polar(polar)
        p.fringe_fitting(True, True)
        p.fringes2db()

    p.delete_uvfits()
except pypima.pima.Error as err:
    print('PIMA Error: ', err)
    sys.exit(1)
except pypima.raexperiment.Error as err:
    print('RaExperiment Error: ', err)
    sys.exit(1)
except:
    print("Unexpected error: ", sys.exc_info()[0])
    raise
