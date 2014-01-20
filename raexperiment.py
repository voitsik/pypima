#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 04:02:35 2013

@author: Petr Voytsik
"""
import datetime
import os.path
import psycopg2
import urllib.request

import pima


class DB:
    """Interface to ra_results DataBase"""
    def __init__(self):
        self.conn = psycopg2.connect(database='ra_results', user='guest',
                                     host='odin', password='tufJoorit3')
        self.conn.autocommit = True

    def get_uvfits_url(self, exper_name, band):
        """Get FITS-file url from DB for given experiment and band"""
        url = None
        url_base = 'ftp://quasars:JcSzi5_k@archive.asc.rssi.ru'

        with self.conn.cursor() as cursor:
            cursor.execute('SELECT path FROM fits_files WHERE exper_name = %s \
                            AND band = %s ORDER BY corr_date DESC;',
                          (exper_name, band.upper()))
            path = cursor.fetchone()

        if path:
            path = path[0]
            url = url_base + path

        return url

    def get_orbit_url(self, exper_name):
        """Returns orbit file url for given experiment"""
        with self.conn.cursor() as cursor:
            cursor.execute("SELECT scf_files.file_name FROM scf_files, \
                            vex_files WHERE vex_files.exper_name = %s AND \
                            scf_files.start_time <= \
                            vex_files.exper_nominal_start AND \
                            scf_files.stop_time >= \
                            vex_files.exper_nominal_stop;", (exper_name,))
            orbit = cursor.fetchone()

        if orbit:
            orbit = orbit[0]

        return orbit

    def get_antab(self, exper_name, band, dist_dir):
        """Download antab-file for the experiment and return path"""
        antab_path = None

        with self.conn.cursor() as cursor:
            cursor.execute("SELECT to_char(exper_nominal_start, 'YYYY_MM_DD') \
                            FROM vex_files WHERE exper_name = %s;",
                          (exper_name,))
            date = cursor.fetchone()

        if date:
            date = date[0]
            date1 = date[0:7]
            url = 'ftp://radik:Zap069_qq@webinet.asc.rssi.ru/radioastron/ampcal/'+\
                date1 + '/' + date + '_' + exper_name + '/' + exper_name + band + \
                '.antab'
            filename = dist_dir + '/' + exper_name + band + '.antab'
            antab_path, _ = urllib.request.urlretrieve(url, filename)

        return antab_path


class RaExperiment(object):
    """This class describe experiment in RadioAstron AGN survey"""
    data_dir_base = '/home/voitsik/data/VLBI/RA/ASC_results/'

    def __init__(self, experiment_code, band, reference_station, uv_fits=None,
                 orbit=None):
                # First, set common variables
        self.exper = experiment_code.lower()
        self.band = band.lower()
        self.sta_ref = reference_station
        self.pima_dir = os.getenv('PIMA_DIR')
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

        # PIMA control file path
        self.cnt_file_name = self.work_dir + '/' + \
            self.exper + '_' + self.band + '_' + 'pima.cnt'

        # Prepare data paths
        self.uv_fits = uv_fits
        self.orbit = orbit

        self.data_dir = self.data_dir_base + self.exper

        # Connect to DB
        self.db = DB()

        # Get fits-file url from DB if don't specified explicitly
        if not self.uv_fits:
            self.uv_fits = self.db.get_uvfits_url(self.exper, self.band)
            if self.uv_fits is None:
                raise Exception('Could not get UV-FITS name from DB for \
                    experiment {}'.format(self.exper))
            local_file_name = '{}/{}'.format(self.data_dir,
                              os.path.basename(self.uv_fits))
            if os.path.isfile(local_file_name):
                self.uv_fits = local_file_name

        if not self.orbit:
            self._get_orbit()
        self._mk_cnt()

        self.pima = pima.Pima(self.exper, self.band, self.work_dir)

    def _mk_cnt(self):
        """Make new cnt-file from template"""
        cnt_templ_name = self.pima_dir + '/share/pima/EXPER_' + self.band + \
            '_pima.cnt'
        cnt_templ = open(cnt_templ_name, 'r')
        cnt_file = open(self.cnt_file_name, 'w')

        for line in cnt_templ:
            if 'CDATE' in line:
                line = line.replace('CDATE', str(datetime.datetime.now()))
            elif line.startswith('SESS_CODE:'):
                line = '{:<20}{}_{}\n'.format('SESS_CODE:', self.exper,
                       self.band)
            elif line.startswith('EXPER_DIR:'):
                line = line.replace('SCRDIR', self.pima_scr)
            elif line.startswith('UV_FITS:'):
                line = '{:<20}{}\n'.format('UV_FITS:', self.uv_fits)
            elif line.startswith('BANDPASS_MASK_FILE:') or \
                    line.startswith('FRINGE_FILE:') or \
                    line.startswith('FRIRES_FILE:'):
                line = line.replace('EXP_DIR', self.exp_dir)
                line = line.replace('EXPERN', self.exper)
            elif line.startswith('STA_REF:'):
                line = '{:<20}{}\n'.format('STA_REF:', self.sta_ref)
            elif line.startswith('BANDPASS_FILE:'):
                bps_file_name = \
                    '{}/{}_{}.bps'.format(self.work_dir, self.exper, self.band)
                line = '{:<20}{}\n'.format('BANDPASS_FILE:', bps_file_name)
            elif line.startswith('EPHEMERIDES_FILE:'):
                line = 'EPHEMERIDES_FILE:   ' + self.orbit + '\n'
            cnt_file.write(line)

        cnt_templ.close()
        cnt_file.close()

    def _download_fits(self):
        """Download FITS-file from remote place"""
        local_file_name = '{}/{}'.format(self.data_dir,
                          os.path.basename(self.uv_fits))

        if not os.path.isfile(self.uv_fits):
            if not os.path.isdir(self.data_dir):
                os.mkdir(self.data_dir)
            print('Info: Start downloading file {}...'.format(self.uv_fits))
            new_file_name, _ = urllib.request.urlretrieve(self.uv_fits,
                               filename=local_file_name)
            print('Info: Done')
            self.uv_fits = new_file_name
        else:
            print('Info: file {} already exists'.format(self.uv_fits))

    def _get_orbit(self):
        """Download reconstructed orbit file from FTP"""
        orbit = self.db.get_orbit_url(self.exper)

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

    def load(self):
        """Download data and run pima load"""
        if not os.path.exists(self.uv_fits):
            self._download_fits()
            self.pima.update_cnt({'UV_FITS:': self.uv_fits})

        if not os.path.isfile(self.antab):
            self.antab = self.db.get_antab(self.exper, self.band,
                                           self.work_dir + '/antab')
            if self.antab is None:
                raise Exception('Could not get antab-file')
            print('DEBUG: antab -> {}'.format(self.antab))

        self.pima.load()

    def fringe_fitting(self, bandpass=False):
        """Do fringe fitting"""
        self.pima.coarse()
        if bandpass:
            self.pima.bpas()
        else:
            self.pima.update_cnt({'BANDPASS_FILE:': 'NO'})
        self.pima.fine()
