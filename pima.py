#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 17:21:29 2013

@author: Petr Voytsik
"""

import os.path


class Pima:
    def __init__(self, exper, band):
        self.exper = exper.lower()
        self.band = band.lower()

        self.exp_dir = os.getenv('pima_exp_dir')
        if not self.exp_dir:
            raise Exception("Environment variable 'pima_exp_dir' is not set")

        work_dir = self.exp_dir + '/' + self.exper + '_auto'
        if not os.path.exists(work_dir):
            os.mkdir(work_dir)

        work_dir_link = self.exp_dir + '/' + self.exper
        if(os.path.isdir(work_dir_link)):
            os.rename(work_dir_link, work_dir_link + '_bak')
        elif(os.path.exists(work_dir_link)):
            os.remove(work_dir_link)
        os.symlink(os.path.basename(work_dir), work_dir_link)

        self.work_dir = work_dir_link
        os.chdir(self.work_dir)

    def load(self):
        return
