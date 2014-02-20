#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 04:02:35 2013

@author: Petr Voytsik
"""

from datetime import datetime
from datetime import timedelta
import netrc
import os.path
import psycopg2
from urllib.error import URLError
import urllib.request as urlreq
import pypima.pima
from pypima.fri import Fri
from pypima.pima import Pima
import time


class Error(Exception):
    """Raised when RaExperiment error occurs"""
    def __init__(self, exper, band, msg):
        self.exper = exper
        self.band = band
        self.msg = msg

    def __str__(self):
        return '{}({}): {}'.format(self.exper, self.band, self.msg)


class DB(object):
    """Interface to ra_results DataBase"""
    def __init__(self, exper, band):
        self.exper = exper
        self.band = band

        self.conn = psycopg2.connect(database='ra_results', user='guest',
                                     host='odin')
        self.connw = psycopg2.connect(database='ra_results', user='editor',
                                      host='odin')
        self.conn.autocommit = True

        nrc = netrc.netrc()
        self.web_login = nrc.authenticators('webinet.asc.rssi.ru')[0]
        self.web_passw = nrc.authenticators('webinet.asc.rssi.ru')[2]
        self.arc_login = nrc.authenticators('archive.asc.rssi.ru')[0]
        self.arc_passw = nrc.authenticators('archive.asc.rssi.ru')[2]

    def get_uvfits_url(self):
        """Get FITS-file url from DB for given experiment and band"""
        url = None
        size = 0
        url_base = 'ftp://{}:{}@archive.asc.rssi.ru'.format(self.arc_login,
                   self.arc_passw)

        with self.conn.cursor() as cursor:
            cursor.execute('SELECT path, size FROM fits_files WHERE \
exper_name = %s AND band = %s ORDER BY corr_date DESC, path DESC;',
                          (self.exper, self.band.upper()))
            repl = cursor.fetchone()

        if repl:
            path = repl[0]
            size = repl[1]
            url = url_base + path

        return url, size

    def get_orbit_url(self):
        """Returns orbit file url for given experiment"""
        url = None
        url_base = 'ftp://{}:{}@webinet.asc.rssi.ru/radioastron/\
oddata/reconstr/'.format(self.web_login, self.web_passw)

        with self.conn.cursor() as cursor:
            cursor.execute("SELECT scf_files.file_name FROM scf_files, \
vex_files WHERE vex_files.exper_name = %s AND \
scf_files.start_time <= vex_files.exper_nominal_start AND \
scf_files.stop_time >= vex_files.exper_nominal_stop;", (self.exper,))
            orbit = cursor.fetchone()

        if orbit:
            orbit_file = orbit[0]
            url = url_base + orbit_file

        return url

    def get_antab_url(self):
        """Download antab-file for the experiment and return path"""
        url = None
        url_base = 'ftp://{}:{}@webinet.asc.rssi.ru/radioastron/\
ampcal'.format(self.web_login, self.web_passw)

        with self.conn.cursor() as cursor:
            cursor.execute("SELECT to_char(exper_nominal_start, 'YYYY_MM_DD') \
                            FROM vex_files WHERE exper_name = %s;",
                          (self.exper,))
            date = cursor.fetchone()

        if date:
            date = date[0]
            date1 = date[0:7]
            url = '{0}/{1}/{2}_{3}/{3}{4}.antab'.format(url_base, date1, date,
                  self.exper, self.band)

        return url

    def _check_and_delete(self, exper, band, polar):
        """
        Delete records if any exists.

        """
        with self.connw.cursor() as cursor:
            cursor.execute("SELECT * FROM pima_observations WHERE \
exper_name = %s AND band = %s AND polar = %s", (exper, band, polar))
            test = cursor.fetchone()
            if test:
                cursor.execute("DELETE FROM pima_observations WHERE \
exper_name = %s AND band = %s AND polar = %s", (exper, band, polar))

    def fri2db(self, fri_file):
        """
        Store information from the fri-file to the DB.

        """
        exper, band = fri_file[0]['session_code'].split('_')

        polar = fri_file[0]['polar']

        self._check_and_delete(exper, band, polar)

        print(fri_file)
        with self.connw.cursor() as cur:
            for rec in fri_file:
                obs = rec['obs']
                start_time = rec['start_time']
                stop_time = rec['stop_time']
            #        exper, band = rec['session_code'].split('_')
                source = rec['source']
                sta1 = rec['sta1']
                sta2 = rec['sta2']
                snr = rec['SNR']
                delay = rec['delay']
                rate = rec['rate']
                if rec['FRIB.FINE_SEARCH'] == 'ACC':
                    accel = rec['ph_acc']
                else:
                    accel = rec['accel']
                ampl = rec['ampl_lsq']
                dur = rec['duration']
                u = rec['U']
                v = rec['V']
                uv_rad_ed = rec['uv_rad_ed']
                freq = rec['ref_freq']

                query = "INSERT INTO pima_observations (obs, start_time, \
stop_time, exper_name, band, source, polar, st1, st2, delay, rate, accel, \
snr, ampl, solint, u, v, base_ed, ref_freq) VALUES (%s, %s, %s, %s, %s, %s, \
%s, %s, %s,  %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);"
                cur.execute(query, (obs, start_time, stop_time, exper, band,
                                    source, polar, sta1, sta2, delay, rate,
                                    accel, snr, ampl, dur, u, v, uv_rad_ed,
                                    freq))
        self.connw.commit()

    def exper_info2db(self, exper_info, uv_fits):
        """
        Put experiment info to the DB.

        """
        with self.connw.cursor() as cursor:
            cursor.execute("SELECT exper_name, band FROM pima_experiments \
WHERE exper_name = %s AND band = %s", (self.exper, self.band))
            test = cursor.fetchone()
            if test is None:
                cursor.execute("INSERT INTO pima_experiments (exper_name, \
band) VALUES (%s, %s);", (self.exper, self.band))

        query = 'UPDATE pima_experiments SET fits_idi = %s, sp_chann_num = %s,\
 time_epochs_num = %s, scans_num = %s, obs_num = %s, uv_points_num = %s, \
uv_points_used_num = %s, deselected_points_num = %s, no_auto_points_num = %s, \
accum_length = %s, utc_minus_tai = %s, nominal_start = %s, nominal_end = %s \
WHERE exper_name = %s AND band = %s'
        with self.connw.cursor() as cursor:
            cursor.execute(query, (uv_fits, exper_info.sp_chann_num,
                                   exper_info.time_epochs_num,
                                   exper_info.scans_num,
                                   exper_info.obs_num,
                                   exper_info.uv_points_num,
                                   exper_info.uv_points_used_num,
                                   exper_info.deselected_points_num,
                                   exper_info.no_auto_points_num,
                                   exper_info.accum_length,
                                   timedelta(seconds=exper_info.utc_minus_tai),
                                   exper_info.nominal_start,
                                   exper_info.nominal_end,
                                   self.exper, self.band))

        self.connw.commit()


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
            self._error('unknown band {}'.format(band))

        self.pima_dir = os.getenv('PIMA_DIR')
        self.exp_dir = os.getenv('pima_exp_dir')
        if not self.exp_dir:
            self._error("Environment variable $pima_exp_dir is not set")

        self.pima_scr = os.getenv('pima_scr_dir')

        # Create working directory and symlink
        work_dir = os.path.join(self.exp_dir, self.exper + '_auto')
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
        antab_dir = os.path.join(self.work_dir, 'antab')
        if not os.path.exists(antab_dir):
            os.mkdir(antab_dir)
        self.antab = os.path.join(antab_dir, self.exper + self.band + '.antab')

        # PIMA control file path
        self.cnt_file_name = os.path.join(self.work_dir, '{}_{}_pima.cnt'.
                                          format(self.exper, self.band))

        # Prepare data paths
        self.uv_fits = uv_fits
        self.orbit = orbit

        # Connect to DB
        self.db = DB(self.exper, self.band)

        self._mk_cnt()

        self.pima = Pima(self.exper, self.band, self.work_dir)

    def _mk_cnt(self):
        """Make new cnt-file from template"""
        cnt_templ_name = os.path.join(self.pima_dir, 'share/pima',
                                      'EXPER_{}_pima.cnt'.format(self.band))

        cnt_templ = open(cnt_templ_name, 'r')
        cnt_file = open(self.cnt_file_name, 'w')

        for line in cnt_templ:
            if 'CDATE' in line:
                line = line.replace('CDATE', str(datetime.now()))
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
        fits_url, size = self.db.get_uvfits_url()
        if not fits_url:
            self._error('Could not find FITS file name in DB')

        # Delete spaces in filename
        self.uv_fits = os.path.join(data_dir, os.path.basename(fits_url).
                                    replace(' ', ''))
        lock_file = self.uv_fits + '.lock'

        if os.path.isfile(self.uv_fits) and \
                os.path.getsize(self.uv_fits) == size:
            self._print_info('file {} already exists'.format(self.uv_fits))
        elif os.path.isfile(lock_file):
            self._print_info('file {} is being downloaded now, wait')
            while os.path.isfile(lock_file):
                print('.', end='')
                time.sleep(10)
            print('')
        else:
            if not os.path.isdir(data_dir):
                os.makedirs(data_dir)
            self._print_info('Start downloading file {}...'.format(fits_url))
            lock_ = open(lock_file, 'w')
            lock_.close()
            try:
                self.uv_fits, _ = urlreq.urlretrieve(fits_url,
                                                     filename=self.uv_fits)
            except URLError as ex:
                os.remove(lock_file)
                self._error('Could not download file {}: {}'.format(
                    fits_url, ex.reason))

            os.remove(lock_file)
            self._print_info('Done')

        self.pima.update_cnt({'UV_FITS:': self.uv_fits})

    def _get_orbit(self):
        """Download reconstructed orbit file from FTP"""
        orbit_url = self.db.get_orbit_url()

        if not orbit_url:
            self._error('Could not find reconstructed orbit')

        self.orbit = os.path.join(self.work_dir, os.path.basename(orbit_url))

        if not os.path.isfile(self.orbit):
            self._print_info('Start downloading orbit file {} ...'.format(
                orbit_url))
            with urlreq.urlopen(orbit_url) as orb_request:
                orb_data = orb_request.read().decode().replace('\r\n', '\n')
            with open(self.orbit, 'w') as orb_file:
                if not orb_data.startswith('CCSDS_OEM_VERS'):
                    orb_file.write('CCSDS_OEM_VERS = 2.0\n')
                orb_file.write(orb_data)
            self._print_info('Done')
        else:
            self._print_info('file {} already exists'.format(self.orbit))

        self.pima.update_cnt({'EPHEMERIDES_FILE:': self.orbit})

    def _get_antab(self):
        """
        Download antab-file from FTP server.

        """
        antab_url = self.db.get_antab_url()
        if antab_url is None:
            self._print_warn('Could not get antab file url from DB.')
        else:
            self._print_info('Start downloading file {}'.format(antab_url))
            try:
                self.antab, _ = urlreq.urlretrieve(antab_url,
                                                   filename=self.antab)
            except URLError as ex:
                self._print_warn('Could not download file {}: {}'.format(
                    antab_url, ex.reason))

    def _print_info(self, msg):
        """Print some information"""
        print('Info: {}({}): {}'.format(self.exper, self.band, msg))

    def _print_warn(self, msg):
        """Print warning"""
        print('Warning: {}({}): {}'.format(self.exper, self.band, msg))

    def _error(self, msg):
        """Raise pima.Error exception"""
        raise Error(self.exper, self.band, msg)

    def load(self, download_only=False):
        """
        Download data, run pima load, and do some checks.

        """
        if self.uv_fits is None:
            self._download_fits()

        if self.orbit is None:
            self._get_orbit()

        if not os.path.isfile(self.antab):
            self._get_antab()

        if download_only:
            return

        if self.band == 'p':
            self.pima.update_cnt({'END_FRQ:': '1'})

        self.pima.load()

        self.db.exper_info2db(self.pima.exper_info,
                              os.path.basename(self.uv_fits))

        if self.pima.obs_number() == 0:
            self._error('ZERO observations have been loaded')

        if 'RADIO-AS' not in self.pima.station_list():
            self._print_warn('RADIO-AS is not in station list')
            self.sta_ref = self.pima.station_list()[0]
            self.pima.update_cnt({'STA_REF:': self.sta_ref})

        desel_nam = self.pima.number_of_deselected_points()
        if desel_nam > 10:
            self._print_warn('Total number of deselected points is ' +
                             desel_nam)

    def _select_ref_sta(self, fri_file):
        """
        Select reference station for bandpass calibration.
        """
        fri = Fri(fri_file)
        self.sta_ref = None
        snr = 0
        obs = fri.max_snr('RADIO-AS')

        if len(obs):
            if obs['SNR'] < 5.6:
                self._print_warn('SNR is too low on space baseline for \
bandpass: ' + str(obs['SNR']))
            else:
                if obs['sta1'] == 'RADIO-AS':
                    self.sta_ref = obs['sta2']
                elif obs['sta2'] == 'RADIO-AS':
                    self.sta_ref = obs['sta1']
        else:
            self._print_info('There is no scans with RADIO-AS')

        if self.sta_ref is None:
            obs = fri.max_snr()
            if obs['SNR'] < 5.6:
                self._print_warn('SNR is too low for bandpass: ' +
                                 str(obs['SNR']))
            else:
                good_stations = ['ARECIBO', 'GBT-VLBA', 'EFLSBERG']
                for sta in good_stations:
                    if sta in [obs['sta1'], obs['sta2']]:
                        self.sta_ref = sta
                        break
                if self.sta_ref is None:
                    self.sta_ref = obs['sta1']

        if self.sta_ref:
            snr = min(10.0, obs['SNR']-0.1)
            self.pima.update_cnt({'STA_REF:': self.sta_ref,
                                  'BPS.SNR_MIN_ACCUM:': str(snr),
                                  'BPS.SNR_MIN_FINE:': str(snr)})
            self._print_info('new reference station is {}'.format(
                             self.sta_ref))
            self._print_info('set SNR_MIN to {:.1f}'.format(snr))
            return True
        else:
            return False

    def fringe_fitting(self, bandpass=False, accel=False):
        """
        Do fringe fitting.

        Parameters
        ----------
        bandpass : bool, optional
            If True try to do bandpass calibration. Default is False.

        accel: boot, optional
            If True turn on phase acceleration fitting.

        """
        if accel:
            self.pima.update_cnt({'FRIB.FINE_SEARCH:': 'ACC',
                                  'PHASE_ACCEL_MIN:': '-1D-14',
                                  'PHASE_ACCEL_MAX:': '1D-14'})

        if self.pima.chan_number() > 512:
            self._print_warn('Too many spectral channels for bandpass: {}'.
                             format(self.pima.chan_number()))
            bandpass = False

        if bandpass:
            fri_file = self.pima.coarse()

            # Now auto select reference station
            if self._select_ref_sta(fri_file):
                try:
                    self.pima.bpas()
                except pypima.pima.Error as err:
                    print(err)
                    self._print_info('Try INIT bandpass')
                    self.pima.bpas(['BPS.MODE:', 'INIT'])

        self.pima.fine()

    def split(self, source=None):
        """
        Do SPLIT.

        """
        pass

    def fringes2db(self):
        """Put fringes information to DB"""
        fri_file = self.pima.cnt_params['FRINGE_FILE:']
        if not os.path.isfile(fri_file):
            return

        self.db.fri2db(Fri(fri_file))

    def delete_uvfits(self):
        """Delete UV-FITS file"""
        if os.path.isfile(self.uv_fits):
            os.remove(self.uv_fits)
