#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 04:02:35 2013

@author: Petr Voytsik
"""
import datetime
import netrc
import os.path
import psycopg2
import re
import urllib.request

import pima


class DB:
    """Interface to ra_results DataBase"""
    def __init__(self):
        self.conn = psycopg2.connect(database='ra_results', user='guest',
                                     host='odin', password='tufJoorit3')
        self.conn.autocommit = True

        nrc = netrc.netrc()
        self.web_login = nrc.authenticators('webinet.asc.rssi.ru')[0]
        self.web_passw = nrc.authenticators('webinet.asc.rssi.ru')[2]
        self.arc_login = nrc.authenticators('archive.asc.rssi.ru')[0]
        self.arc_passw = nrc.authenticators('archive.asc.rssi.ru')[2]

    def get_uvfits_url(self, exper_name, band):
        """Get FITS-file url from DB for given experiment and band"""
        url = None
        url_base = 'ftp://{}:{}@archive.asc.rssi.ru'.format(self.arc_login,
                   self.arc_passw)

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
        url = None
        url_base = 'ftp://{}:{}@webinet.asc.rssi.ru/radioastron/\
oddata/reconstr/'.format(self.web_login, self.web_passw)

        with self.conn.cursor() as cursor:
            cursor.execute("SELECT scf_files.file_name FROM scf_files, \
vex_files WHERE vex_files.exper_name = %s AND \
scf_files.start_time <= vex_files.exper_nominal_start AND \
scf_files.stop_time >= vex_files.exper_nominal_stop;", (exper_name,))
            orbit = cursor.fetchone()

        if orbit:
            orbit_file = orbit[0]
            url = url_base + orbit_file

        return url

    def get_antab_url(self, exper, band):
        """Download antab-file for the experiment and return path"""
        url = None
        url_base = 'ftp://{}:{}@webinet.asc.rssi.ru/radioastron/\
ampcal'.format(self.web_login, self.web_passw)

        with self.conn.cursor() as cursor:
            cursor.execute("SELECT to_char(exper_nominal_start, 'YYYY_MM_DD') \
                            FROM vex_files WHERE exper_name = %s;",
                          (exper,))
            date = cursor.fetchone()

        if date:
            date = date[0]
            date1 = date[0:7]
            url = '{0}/{1}/{2}_{3}/{3}{4}.antab'.format(url_base, date1, date,
                  exper, band)

        return url


class RaExperiment(object):
    """This class describe experiment in RadioAstron AGN survey"""

    def __init__(self, experiment_code, band, reference_station=None,
                 uv_fits=None, orbit=None):
                # First, set common variables
        self.exper = experiment_code.lower()
        self.band = band.lower()
        self.sta_ref = reference_station
        if self.sta_ref is None:
            self.sta_ref = 'RADIO-AS'

        if not self.band in ['p', 'l', 'c', 'k']:
            raise Exception('Error: unknown band {}'.format(band))

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

        # Connect to DB
        self.db = DB()

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
            elif line.startswith('UV_FITS:') and self.uv_fits:
                line = '{:<20}{}\n'.format('UV_FITS:', self.uv_fits)
            elif line.startswith('BANDPASS_MASK_FILE:') or \
                    line.startswith('FRINGE_FILE:') or \
                    line.startswith('FRIRES_FILE:'):
                line = line.replace('EXP_DIR', self.exp_dir)
                line = line.replace('EXPERN', self.exper)
            elif line.startswith('STA_REF:') and self.sta_ref:
                line = '{:<20}{}\n'.format('STA_REF:', self.sta_ref)
            elif line.startswith('BANDPASS_FILE:'):
                line = '{:<20}{}\n'.format('BANDPASS_FILE:', 'NO')
            elif line.startswith('EPHEMERIDES_FILE:') and self.orbit:
                line = '{:<20}{}\n'.format('EPHEMERIDES_FILE:', self.orbit)
            cnt_file.write(line)

        cnt_templ.close()
        cnt_file.close()

    def _download_fits(self):
        """Download FITS-file from remote place"""
        data_dir = '/home/voitsik/data/VLBI/RA/ASC_results/' + self.exper
        fits_url = self.db.get_uvfits_url(self.exper, self.band)
        if not fits_url:
            raise Exception('Error: Could not find FITS file name in DB for \
 the experiment {}({})'.format(self.exper, self.band))

        self.uv_fits = '{}/{}'.format(data_dir,
                                      os.path.basename(fits_url))

        if not os.path.isfile(self.uv_fits):
            if not os.path.isdir(data_dir):
                os.mkdir(data_dir)
            print('Info: Start downloading file {}...'.format(fits_url))
            self.uv_fits, _ = urllib.request.urlretrieve(fits_url,
                                                         filename=self.uv_fits)
            print('Info: Done')
        else:
            print('Info: file {} already exists'.format(self.uv_fits))

        self.pima.update_cnt({'UV_FITS:': self.uv_fits})

    def _get_orbit(self):
        """Download reconstructed orbit file from FTP"""
        orbit_url = self.db.get_orbit_url(self.exper)

        if not orbit_url:
            raise Exception("Can't find reconstructed orbit for experiment \
{}".format(self.exper))

        self.orbit = self.work_dir + '/' + os.path.basename(orbit_url)

        if not os.path.isfile(self.orbit):
            print('Info: Start downloading orbit file {} ...'.format(
                orbit_url))
            with urllib.request.urlopen(orbit_url) as orb_request:
                orb_data = orb_request.read().decode().replace('\r\n', '\n')
            with open(self.orbit, 'w') as orb_file:
                if not orb_data.startswith('CCSDS_OEM_VERS'):
                    orb_file.write('CCSDS_OEM_VERS = 2.0\n')
                orb_file.write(orb_data)
            print('Info: Done.')
        else:
            print('Info: file {} already exists'.format(self.orbit))

        self.pima.update_cnt({'EPHEMERIDES_FILE:': self.orbit})

    def load(self):
        """Download data and run pima load"""
        if self.uv_fits is None:
            self._download_fits()

        if self.orbit is None:
            self._get_orbit()

        if not os.path.isfile(self.antab):
            antab_url = self.db.get_antab_url(self.exper, self.band)
            if antab_url is None:
                print('Warning: Could not get antab file url')
            else:
                try:
                    print('Info: Start downloading file {}'.format(antab_url))
                    self.antab, _ = urllib.request.urlretrieve(antab_url,
                        filename=self.antab)
                except urllib.error.URLError as ex:
                    print('Warning: could not download file: {}'.format(
                        ex.reason))
#            print('DEBUG: antab -> {}'.format(self.antab))

        self.pima.load()

    def fringe_fitting(self, bandpass=False):
        """Do fringe fitting"""
        if bandpass:
            self.pima.coarse()

            # Now auto select reference station
            coarse_log_file = '{}/{}_{}_coarse.log'.format(self.work_dir,
                              self.exper, self.band)
            snrs = []
            with open(coarse_log_file, 'r') as f:
                for line in f:
                    if not line.startswith(
                            self.pima.cnt_params['SESS_CODE:']):
                        continue
#                        print('DEBUG: line = {}'.format(line))
                    sta1 = line.split('/')[0].split()[-1]
                    sta2 = line.split('/')[1].split()[0]
                    snr = float(line.split('SNR=')[1])

                    if sta1 == 'RADIO-AS':
                        snrs.append((float(snr), sta2))
                    elif sta2 == 'RADIO-AS':
                        snrs.append((float(snr), sta1))

            self.sta_ref = sorted(snrs, reverse=True)[0][1]
            max_snr = sorted(snrs, reverse=True)[0][0]
            snr = max(5.5, min(10.0, max_snr-0.1))

            if not self.sta_ref:
                raise Exception('Error: could not select reference station')
            self.pima.update_cnt({'STA_REF:': self.sta_ref,
                                  'BPS.SNR_MIN_ACCUM:': str(snr),
                                  'BPS.SNR_MIN_FINE:': str(snr)})

            print('Info: new reference station is {}'.format(self.sta_ref))
            print('Info: set SNR_MIN to {:.1f}'.format(snr))

            if max_snr < 5.5:
                print('Warning: SNR on space baselines is too low for \
bandpass calibration')
            else:
                self.pima.bpas()

        self.pima.fine()
