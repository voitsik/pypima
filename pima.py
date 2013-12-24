#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 17:21:29 2013

@author: Petr Voytsik
"""

import datetime
import os.path
import subprocess
import psycopg2
import urllib.request


def get_orbit_url(exper_name):
    """Returns orbit file url for given experiment"""
    conn = psycopg2.connect(database='ra_results', user='guest', host='odin',
                            password='tufJoorit3')
    conn.autocommit = True
    cur = conn.cursor()

    cur.execute("SELECT scf_files.file_name FROM scf_files, vex_files WHERE \
          vex_files.exper_name = %s AND \
          scf_files.start_time <= vex_files.exper_nominal_start AND \
          scf_files.stop_time >= vex_files.exper_nominal_stop;", (exper_name,))
    orbit = cur.fetchone()
    cur.close()
    conn.close()

    if orbit:
        orbit = orbit[0]

    return orbit


def get_antab(exper_name, band, dist_dir):
    """Download antab-file for the experiment and return path"""
    antab_path = None

    conn = psycopg2.connect(database='ra_results', user='guest', host='odin',
                            password='tufJoorit3')
    conn.autocommit = True
    cur = conn.cursor()

    cur.execute("SELECT to_char(exper_nominal_start, 'YYYY_MM_DD') \
        FROM vex_files WHERE exper_name = %s;", (exper_name,))
    date = cur.fetchone()
    cur.close()
    conn.close()

    if date:
        date = date[0]
        date1 = date[0:7]
        url = 'ftp://radik:Zap069_qq@webinet.asc.rssi.ru/radioastron/ampcal/'+\
            date1 + '/' + date + '_' + exper_name + '/' + exper_name + band + \
            '.antab'
        filename = dist_dir + '/' + exper_name + band + '.antab'
        antab_path, _ = urllib.request.urlretrieve(url, filename)

    return antab_path


class Pima:
    data_base = '/home/voitsik/data/VLBI/RA/ASC_results/'
    ftp_base = 'ftp://quasars:JcSzi5_k@archive.asc.rssi.ru'

    def __init__(self, exper, band, uv_fits=None, url=None, orbit=None):
        # First, set common variables
        self.exper = exper.lower()
        self.band = band.lower()

        self.pima_dir = os.getenv('PIMA_DIR')
        if not self.pima_dir:
            raise Exception('Environment variable $PIMA_DIR is not set')

        self.pima_exec = self.pima_dir + '/bin/pima'

        self.exp_dir = os.getenv('pima_exp_dir')
        if not self.exp_dir:
            raise Exception("Environment variable $pima_exp_dir is not set")

        self.pima_scr = os.getenv('pima_scr_dir')

        # Create working directory and symlink
        work_dir = self.exp_dir + '/' + self.exper + '_auto'
        if not os.path.exists(work_dir):
            os.mkdir(work_dir)

        work_dir_link = self.exp_dir + '/' + self.exper
        if os.path.islink(work_dir_link):
            os.remove(work_dir_link)
        elif os.path.exists(work_dir_link):
            os.rename(work_dir_link, work_dir_link + '_bak')

        os.symlink(os.path.basename(work_dir), work_dir_link)

        self.work_dir = work_dir_link
        os.chdir(self.work_dir)

        # Make directory for antabs
        antab_dir = self.work_dir + '/antab'
        if not os.path.exists(antab_dir):
            os.mkdir(antab_dir)
        self.antab = antab_dir + '/' + self.exper + self.band + '.antab'

        self.cnt_file_name = self.work_dir + '/' + \
            self.exper + '_' + self.band + '_' + 'pima.cnt'

        # Prepare data paths
        self.url = url
        self.uv_fits = uv_fits
        self.orbit = orbit

        self.data_dir = self.data_base + self.exper
        if not self.uv_fits:
            self._set_fits_name()
        if not self.orbit:
            self._get_orbit()
        self._mk_cnt()

    def _set_fits_name(self):
        if not self.url:
            conn = psycopg2.connect(database='ra_results', user='voitsik')
            conn.autocommit = True
            cur = conn.cursor()

            cur.execute('SELECT path FROM fits_files WHERE exper_name = %s \
                         AND band = %s ORDER BY corr_date DESC;',
                       (self.exper, self.band.upper()))
            path = cur.fetchone()
            cur.close()
            conn.close()

            if not path:
                raise Exception('Could not get UV-FITS name from DB for the \
                experiment {}, band {}'.format(self.exper, self.band))
            else:
                path = path[0]
            self.url = self.ftp_base + path

        self.uv_fits = self.data_dir + '/' + os.path.basename(self.url)

    def _mk_cnt(self):
        cnt_templ_name = self.pima_dir + '/share/pima/EXPER_' + self.band + \
            '_pima.cnt'
        cnt_templ = open(cnt_templ_name, 'r')
        cnt_file = open(self.cnt_file_name, 'w')

        for line in cnt_templ:
            if 'CDATE' in line:
                line = line.replace('CDATE', str(datetime.datetime.now()))
            elif line.startswith('SESS_CODE:'):
                line = 'SESS_CODE:      ' + self.exper + '_' + self.band + '\n'
            elif line.startswith('EXPER_DIR:'):
                line = line.replace('SCRDIR', self.pima_scr)
            elif line.startswith('UV_FITS:'):
                line = 'UV_FITS:            ' + self.uv_fits + '\n'
            elif line.startswith('BANDPASS_FILE:') or \
                    line.startswith('BANDPASS_MASK_FILE:') or \
                    line.startswith('FRINGE_FILE:') or \
                    line.startswith('FRIRES_FILE:'):
                line = line.replace('EXP_DIR', self.exp_dir)
                line = line.replace('EXPERN', self.exper)
            elif line.startswith('EPHEMERIDES_FILE:'):
                line = 'EPHEMERIDES_FILE:   ' + self.orbit + '\n'
            cnt_file.write(line)

        cnt_templ.close()
        cnt_file.close()

    def _download_fits(self):
        if not os.path.isdir(self.data_dir):
            os.mkdir(self.data_dir)

        if not os.path.isfile(self.uv_fits):
            print('Info: Start downloading file {}...'.format(self.uv_fits))
            urllib.request.urlretrieve(self.url, filename=self.uv_fits)
            print('Info: Done')
        else:
            print('Info: file {} already exists'.format(self.uv_fits))

    def _get_orbit(self):
        """Download reconstructed orbit file from FTP"""
        orbit = get_orbit_url(self.exper)

        if not orbit:
            raise Exception("Can't find reconstructed orbit for experiment \
            {}".format(self.exper))

        self.orbit = self.work_dir + '/' + orbit

        if not os.path.isfile(self.orbit):
            print('Info: Start downloading orbit file')
            user = 'radik'
            passw = 'Zap069_qq'
            orbit_url = 'ftp://' + user + ':' + passw + \
                        '@webinet.asc.rssi.ru/radioastron/oddata/reconstr/' + \
                        orbit
            orb_request = urllib.request.urlopen(orbit_url)
            orb_data = orb_request.read().decode().replace('\r\n', '\n')
            orb_file = open(self.orbit, 'w')
            if not orb_data.startswith('CCSDS_OEM_VERS'):
                orb_file.write('CCSDS_OEM_VERS = 2.0\n')
            orb_file.write(orb_data)
            orb_file.close()
            print('Info: Done.')
        else:
            print('Info: file {} already exists'.format(self.orbit))

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
        if not os.path.exists(self.uv_fits):
            self._download_fits()

        if not os.path.isfile(self.antab):
            self.antab = get_antab(self.exper, self.band, self.work_dir +
                '/antab')
            if self.antab is None:
                raise Exception('Could not get antab-file')
            print('DEBUG: antab -> {}'.format(self.antab))

        log_name = self.exper + '_' + self.band + '_load.log'
        opts = ['STA_REF:', 'BADARY']
        ret = self._exec('load', log_name=log_name, options=opts)
        if ret:
            raise Exception('load failed with code {}'.format(ret))
