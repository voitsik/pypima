#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 17:50:20 2013

@author: Petr Voytsik
"""

import raexperiment

p = raexperiment.RaExperiment('raes03dz', 'c', 'BADARY')
p.load()
p.fringe_fitting()
