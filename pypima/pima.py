#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 17:21:29 2013

@author: Petr Voytsik
"""

import datetime
import glob
import os.path
import shutil
import subprocess


class Error(Exception):
    """Raised when PIMA error occurs"""
    def __init__(self, exper, band, msg):
        self.exper = exper
        self.band = band
        self.msg = msg

    def __str__(self):
        return '{}({}): {}'.format(self.exper, self.band, self.msg)


class Pima(object):
    """Pima class is analog of pima_fringe.csh script"""
    def __init__(self, experiment_code, band, work_dir=None):
        # First, set common variables
        self.exper = experiment_code.lower()
        self.band = band.lower()
        self.work_dir = work_dir
        if work_dir is None:
            self.work_dir = os.getcwd()

        self.pima_dir = os.getenv('PIMA_DIR')
        if not os.path.isdir(self.pima_dir):
            raise Error(self.exper, self.band, 'Could not find PIMA directory.\
 Please, set $PIMA_DIR environment variable.')
        self.pima_exec = self.pima_dir + '/bin/pima'
        if not os.path.isfile(self.pima_exec):
            raise Error(self.exper, self.band, 'Could not find pima \
executable')

        # PIMA control file path
        self.cnt_file_name = '{}/{}_{}_pima.cnt'.format(self.work_dir,
                             self.exper, self.band)

        if not os.path.isfile(self.cnt_file_name):
            raise Error(self.exper, self.band, 'Could not find control file \
{}'.format(self.cnt_file_name))

        # Dictionary with all parameters from cnt-file
        self.cnt_params = {}
        self._update_cnt_params()

    def _update_cnt_params(self):
        """Read cnt-file and update cnt_params dictionary"""
        with open(self.cnt_file_name, 'r') as cnt_file:
            for line in cnt_file:
                line = line.split('#')[0].strip()
                if len(line) < 8:
                    continue
                key, val = line.split()
                self.cnt_params[key] = val

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
        self._update_cnt_params()

    def get_cnt_params(self, opts):
        """Get parameters name list 'opts' and return dictionary"""
        ret = dict()

        with open(self.cnt_file_name, 'r') as f:
            for line in f:
                if line.startswith('#') or len(line) < 8:
                    continue
                key, val = line.strip().split(None, 1)
                if key in opts:
                    ret[key] = val

        return ret

    def _exec(self, operation, options=None, log_name=None):
        """Execute PIMA binary"""
        cmd_line = [self.pima_exec, self.cnt_file_name, operation]
        if options:
            cmd_line.extend(options)
#        print('DEBUG: cmd_line = {}'.format(cmd_line))

        if log_name is None:
            log_name = self.exper + '_' + self.band + '_' + operation + '.log'

        log = open(log_name, 'w')
        log.write(str(datetime.datetime.now()) + '\n\n')
        log.flush()

        ret = subprocess.call(cmd_line, stdout=log, universal_newlines=True)

        log.write('\n' + str(datetime.datetime.now()) + '\n\n')
        log.close()

        return ret

    def _print_info(self, msg):
        """Print some information"""
        print('Info: {}({}): {}'.format(self.exper, self.band, msg))

    def _error(self, msg):
        """Raise pima.Error exception"""
        raise Error(self.exper, self.band, msg)

    def load(self):
        """Run pima load"""
        # Delete existing PIMA auxiliary files
        auxiliary_files = glob.glob('{}/{}*'.format(
                                    self.cnt_params['EXPER_DIR:'],
                                    self.cnt_params['SESS_CODE:']))
#        print('DEBUG: auxiliary_files = {}'.format(auxiliary_files))
        for aux_file in auxiliary_files:
            if os.path.isfile(aux_file):
                os.remove(aux_file)
            elif os.path.isdir(aux_file):
                shutil.rmtree(aux_file)

        log_name = self.exper + '_' + self.band + '_load.log'
        opts = ['BANDPASS_FILE:', 'NO',
                'POLARCAL_FILE:', 'NO']
        ret = self._exec('load', log_name=log_name, options=opts)
        if ret:
            self._error('load failed with code {}'.format(ret))

        self._print_info('load ok')

    def coarse(self, params=None):
        """Do coarse fringe fitting"""
        log_name = '{}/{}_{}_coarse.log'.format(self.work_dir, self.exper,
                   self.band)
        fri_file = '{}/{}_{}_nobps.fri'.format(self.work_dir, self.exper,
                   self.band)
        if os.path.isfile(fri_file):
            os.remove(fri_file)
        frr_file = '{}/{}_{}_nobps.frr'.format(self.work_dir, self.exper,
                   self.band)
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
        if params:
            opts.extend(params)
        ret = self._exec('frib', opts, log_name)
        if ret:
            self._error('coarse failed with code {}'.format(ret))

        self._print_info('coarse ok')

    def fine(self, params=None):
        """Do file fringe fitting"""
        log_name = '{}/{}_{}_fine.log'.format(self.work_dir, self.exper,
                   self.band)
        fri_file = '{}/{}_{}.fri'.format(self.work_dir, self.exper, self.band)
        if os.path.isfile(fri_file):
            os.remove(fri_file)
        frr_file = '{}/{}_{}.frr'.format(self.work_dir, self.exper, self.band)
        if os.path.isfile(frr_file):
            os.remove(frr_file)

        self.update_cnt({'FRINGE_FILE:': fri_file,
                         'FRIRES_FILE:': frr_file})
        ret = self._exec('frib', params, log_name=log_name)
        if ret:
            self._error('fine failed with code {}'.format(ret))

        self._print_info('fine ok')

    def bpas(self, params=None):
        """Do bandpass calibration"""
        fri_file = '{}/{}_{}_nobps.fri'.format(self.work_dir, self.exper,
                   self.band)
        log_file = '{}/{}_{}_bps.log'.format(self.work_dir, self.exper,
                   self.band)
        exc_obs_file = '{}/{}_{}_bpas_obs.exc'.format(self.work_dir,
                       self.exper, self.band)

        if self.cnt_params['BANDPASS_FILE:'] == 'NO':
            bps_file = '{}/{}_{}.bps'.format(self.work_dir, self.exper,
                       self.band)
            self.update_cnt({'BANDPASS_FILE:': bps_file})

        opts = ['FRINGE_FILE:', fri_file,
                'DEBUG_LEVEL:', '3']

        if os.path.isfile(exc_obs_file):
            opts.extend(['EXCLUDE_OBS_FILE:', exc_obs_file])

        if self.cnt_params['POLAR:'] in ['I', 'RPL']:
            opts.extend(['POLAR:', 'RR'])

        ret = self._exec('bpas', opts, log_file)
        if ret:
            self._error('bpas failed with code {}'.format(ret))

        self._print_info('bpas ok')

    def split(self, tim_mseg=1, params=None):
        """Do SPLIT"""
        opts = ['SPLT.TIM_MSEG:', str(tim_mseg)]
        if params:
            opts.extend(params)
        ret = self._exec('splt', opts)

        if ret:
            self._error('splt failed with code {}'.format(ret))

        self._print_info('split ok')

    def load_gains(self, gain_file, params=None):
        """Load Gains from gain_file"""
        opts = ['evn_gain', gain_file, 'DEBUG_LEVEL:', '6']
        if params:
            opts.extend(params)

        log_file = '{}/{}_{}_gain.log'.format(self.work_dir, self.exper,
                   self.band)

        ret = self._exec('gean', opts, log_file)
        if ret:
            self._error('evn_gain failed with code {}'.format(ret))

    def load_tsys(self, tsys_file, params=None):
        """Load Tsys from tsys_file"""
        opts = ['vlba_log_file', tsys_file, 'DEBUG_LEVEL:', '1']
        if params:
            opts.extend(params)

        log_file = '{}/{}_{}_tsys.log'.format(self.work_dir, self.exper,
                   self.band)

        ret = self._exec('gean', opts, log_file)
        if ret:
            self._error('vlba_log_file failed with code {}'.format(ret))

    # Additional useful utilites
    def ap_minmax(self):
        """Get minimum accummulation period in experiment"""
        ap_min = ap_max = 0
        stt_file = '{}/{}.stt'.format(self.cnt_params['EXPER_DIR:'],
                   self.cnt_params['SESS_CODE:'])

        if os.path.isfile(stt_file):
            with open(stt_file, 'r') as f:
                for line in f:
                    if line.startswith('Accummulation period length_min'):
                        ap_min = float(line.split()[3])
                    elif line.startswith('Accummulation period length_max'):
                        ap_max = float(line.split()[3])

        return ap_min, ap_max

    def number_of_deselected_points(self):
        """Total number of deselected points"""
        num = 9999999
        stt_file = '{}/{}.stt'.format(self.cnt_params['EXPER_DIR:'],
                                      self.cnt_params['SESS_CODE:'])
        if os.path.isfile(stt_file):
            with open(stt_file, 'r') as fil:
                for line in fil:
                    if line.startswith('Total number of deselected points:'):
                        num = int(line.split(':')[1])

        return num

    def sta_list(self):
        """Get station list"""
        sta_file = '{}/{}.sta'.format(self.cnt_params['EXPER_DIR:'],
                                      self.cnt_params['SESS_CODE:'])

        sta_l = []
        if os.path.isfile(sta_file):
            with open(sta_file, 'r') as f:
                for line in f:
                    toks = line.split()
                    sta_l.append(toks[3])

        return sta_l

    def obs_number(self):
        """Get number of observation in experiment"""
        num = 0
        obs_file = '{}/{}.obs'.format(self.cnt_params['EXPER_DIR:'],
                                      self.cnt_params['SESS_CODE:'])

        if os.path.isfile(obs_file):
            with open(obs_file, 'r') as fil:
                for line in fil:
                    if line.startswith('#'):
                        continue
                    num = num + 1

        return num


def fits_to_txt(fits_file):
    """Dump visibilites from 'fits_file' using fits_to_radplot utility"""
    if not os.path.isfile(fits_file):
        raise IOError(2, 'No such file: {}'.format(fits_file))

    txt_file = fits_file.replace('.fits', '.txt')
    cmd_line = ['fits_to_radplot', '-o', txt_file, fits_file]
    out = subprocess.check_output(cmd_line, universal_newlines=True)

    return out
