#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 17:50:20 2013

@author: Petr Voytsik
"""

import raexperiment
import sys

if len(sys.argv) != 3:
    print('Usage: {} <exper> <band>'.format(sys.argv[0]))
    sys.exit(1)

try:
    p = raexperiment.RaExperiment(sys.argv[1], sys.argv[2])
    p.load()
    #print(p.pima.sta_list())
    p.fringe_fitting(True)
except Exception as ex:
    print('Exception: {}'.format(ex))
    sys.exit(1)
