#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 17:50:20 2013

@author: Petr Voytsik
"""

import sys
sys.path.insert(0, '/home/voitsik/prog/pypima')
import pypima.pima
from pypima.raexperiment import RaExperiment

if len(sys.argv) != 3:
    print('Usage: {} <exper> <band>'.format(sys.argv[0]))
    sys.exit(1)

try:
    p = RaExperiment(sys.argv[1], sys.argv[2])
    p.load()
    #print(p.pima.sta_list())
    p.fringe_fitting(True)
    p.fringes2db()
except pypima.pima.Error as err:
    print('PIMA Error: {}'.format(err))
    sys.exit(1)
except Exception as ex:
    print('Exception: {}'.format(ex))
    sys.exit(2)
