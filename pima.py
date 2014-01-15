#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 17:21:29 2013

@author: Petr Voytsik
"""

import datetime
import os.path
import subprocess


class Pima(object):
    """Pima class is analog of pima_fringe.csh script"""
    def __init__(self, experiment_code, band):
        # First, set common variables
        self.exper = experiment_code.lower()
        self.band = band.lower()

        self.pima_dir = os.getenv('PIMA_DIR')
        if not os.path.isdir(self.pima_dir):
            raise Exception('Could not find PIMA directory. \
                Please, set $PIMA_DIR environment variable.')
        self.pima_exec = self.pima_dir + '/bin/pima'
        if not os.path.isfile(self.pima_exec):
            raise Exception('Could not find pima executable')

        # PIMA control file path
        self.cnt_file_name = self.exper + '_' + self.band + '_' + 'pima.cnt'

        if not os.path.isfile(self.cnt_file_name):
            raise Exception('Could not find control file \
                {}'.format(self.cnt_file_name))

    def update_cnt(self, opts):
        """Update pima cnt-file according 'opts' dictionary"""
        if opts is None:
            return

        old_cnt = open(self.cnt_file_name, 'r')
        new_cnt = open(self.cnt_file_name + '.new', 'w')

        for line in old_cnt:
            key = line.split()[0].strip()
            if key in opts.keys():
                line = '{:<20}{}\n'.format(key, opts[key])
            elif line.startswith('# Last update on'):
                line = '# Last update on  \
                    {}\n'.format(str(datetime.datetime.now()))
            new_cnt.write(line)

        old_cnt.close()
        new_cnt.close()
        os.rename(self.cnt_file_name + '.new', self.cnt_file_name)

    def _exec(self, operation, options=None, log_name=None):
        """Execute PIMA binary"""
        cmd_line = [self.pima_exec, self.cnt_file_name, operation]
        if options:
            cmd_line.extend(options)
        print('DEBUG: cmd_line = {}'.format(cmd_line))

        if log_name is None:
            log_name = self.exper + '_' + self.band + '_' + operation + '.log'

        log = open(log_name, 'w')
        log.write(str(datetime.datetime.now()) + '\n\n')
        log.flush()

        ret = subprocess.call(cmd_line, stdout=log, universal_newlines=True)

        log.write('\n' + str(datetime.datetime.now()) + '\n\n')
        log.close()

        return ret

    def load(self):
        """Run pima load"""
        log_name = self.exper + '_' + self.band + '_load.log'
        opts = ['BANDPASS_FILE:', 'NO',
                'POLARCAL_FILE:', 'NO']
        ret = self._exec('load', log_name=log_name, options=opts)
        if ret:
            raise Exception('load failed with code {}'.format(ret))

    def coarse(self):
        """Do coarse fringe fitting"""
        log_name = self.exper + '_' + self.band + '_coarse.log'
        fri_file = '{}_{}_nobps.fri'.format(self.exper, self.band)
        if os.path.isfile(fri_file):
            os.remove(fri_file)
        frr_file = '{}_{}_nobps.frr'.format(self.exper, self.band)
        if os.path.isfile(frr_file):
            os.remove(frr_file)

        opts = ['FRINGE_FILE:', fri_file,
                'FRIRES_FILE:', frr_file,
                'BANDPASS_USE:', 'NO',
                'BANDPASS_FILE:', 'NO',
                'POLARCAL_FILE:', 'NO',
                'FRIB.OVERSAMPLE_MD:', '1',
                'FRIB.OVERSAMPLE_RT:', '1',
                'FRIB.SECONDARY_SNR_MIN:', '0.0',
                'FRIB.SECONDARY_MAX_TRIES:', '0',
                'FRIB.FINE_SEARCH:', 'PAR',
                'MKDB.FRINGE_ALGORITHM:', 'DRF']
        ret = self._exec('frib', opts, log_name)
        if ret:
            raise Exception('coarse failed with code {}'.format(ret))

    def fine(self):
        """Do file fringe fitting"""
        log_name = self.exper + '_' + self.band + '_fine.log'
        fri_file = '{}_{}.fri'.format(self.exper, self.band)
        if os.path.isfile(fri_file):
            os.remove(fri_file)
        frr_file = '{}_{}.frr'.format(self.exper, self.band)
        if os.path.isfile(frr_file):
            os.remove(frr_file)

        opts = ['FRINGE_FILE:', fri_file,
                'FRIRES_FILE:', frr_file]
        ret = self._exec('frib', opts, log_name)
        if ret:
            raise Exception('fine failed with code {}'.format(ret))
