# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 04:02:35 2013

@author: Petr Voytsik
"""

from datetime import datetime
from datetime import timedelta
import glob
import netrc
import os.path
import psycopg2
from urllib.error import URLError
import urllib.request as urlreq
import pypima.pima
from pypima.fri import Fri
from pypima.pima import Pima
import sys
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

        self.conn = None
        self.connw = None
        self.connected = False

        nrc = netrc.netrc()
        self.web_login = nrc.authenticators('webinet.asc.rssi.ru')[0]
        self.web_passw = nrc.authenticators('webinet.asc.rssi.ru')[2]
        self.arc_login = nrc.authenticators('archive.asc.rssi.ru')[0]
        self.arc_passw = nrc.authenticators('archive.asc.rssi.ru')[2]

    def _connect(self):
        """Realy connect to DB"""
        if not self.conn:
            self.conn = psycopg2.connect(database='ra_results', user='guest',
                                         host='odin')
            self.conn.autocommit = True
        if not self.connw:
            self.connw = psycopg2.connect(database='ra_results', user='editor',
                                          host='odin')
        self.connected = True

    def close(self):
        """Close connections"""
        if self.connected:
            self.conn.close()
            self.connw.close()

    def get_uvfits_url(self):
        """Get FITS-file url from DB for given experiment and band"""
        if not self.connected:
            self._connect()

        url = None
        size = 0
        url_base = 'ftp://{}:{}@archive.asc.rssi.ru'.format(self.arc_login,
                   self.arc_passw)

        with self.conn.cursor() as cursor:
            cursor.execute('SELECT path, size FROM fits_files WHERE \
LOWER(exper_name) = LOWER(%s) AND LOWER(band) = LOWER(%s) \
ORDER BY corr_date DESC, path DESC;',
                          (self.exper, self.band))
            reply = cursor.fetchone()

        if reply:
            path = reply[0]
            size = reply[1]
            url = url_base + path

        return url, size

    def get_orbit_url(self):
        """Returns orbit file url for given experiment"""
        if not self.connected:
            self._connect()

        url = None
        url_base = 'ftp://{}:{}@webinet.asc.rssi.ru/radioastron/\
oddata/reconstr/'.format(self.web_login, self.web_passw)

        with self.conn.cursor() as cursor:
            cursor.execute("SELECT scf_files.file_name FROM scf_files, \
vex_files WHERE vex_files.exper_name = %s AND \
scf_files.start_time <= vex_files.exper_nominal_start AND \
scf_files.stop_time >= vex_files.exper_nominal_stop;", (self.exper,))
            reply = cursor.fetchone()

        if reply:
            orbit_file = reply[0]
            url = url_base + orbit_file

        return url

    def get_antab_url(self):
        """Download antab-file for the experiment and return path"""
        if not self.connected:
            self._connect()

        url = None
        url_base = 'ftp://{}:{}@webinet.asc.rssi.ru/radioastron/\
ampcal'.format(self.web_login, self.web_passw)

        with self.conn.cursor() as cursor:
            cursor.execute("SELECT to_char(exper_nominal_start, 'YYYY_MM_DD') \
                            FROM vex_files WHERE exper_name = %s;",
                          (self.exper,))
            reply = cursor.fetchone()

        if reply:
            date = reply[0]
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
            reply = cursor.fetchone()
            if reply:
                cursor.execute("DELETE FROM pima_observations WHERE \
exper_name = %s AND band = %s AND polar = %s", (exper, band, polar))

    def fri2db(self, fri_file, exper_info):
        """
        Store information from the fri-file to the DB.

        """
        if not self.connected:
            self._connect()

        polar = fri_file[0]['polar']

        self._check_and_delete(self.exper, self.band, polar)

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
                cur.execute(query, (obs, start_time, stop_time, self.exper,
                                    self.band,
                                    source, polar, sta1, sta2, delay, rate,
                                    accel, snr, ampl, dur, u, v, uv_rad_ed,
                                    freq))

            # Update status of the observations
            query = "UPDATE pima_observations SET status = %s WHERE \
exper_name = %s AND band = %s AND polar = %s AND snr <= %s;"
            cur.execute(query, ('n', self.exper, self.band, polar, 5.3))

            query = "UPDATE pima_observations SET status = %s WHERE \
exper_name = %s AND band = %s AND polar = %s AND snr >= %s;"
            if exper_info.sp_chann_num <= 128:
                cur.execute(query, ('y', self.exper, self.band, polar, 5.7))
            else:
                cur.execute(query, ('y', self.exper, self.band, polar, 7.0))

        self.connw.commit()

    def exper_info2db(self, exper_info, uv_fits):
        """
        Put experiment info to the DB.

        """
        if not self.connected:
            self._connect()

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
accum_length = %s, utc_minus_tai = %s, nominal_start = %s, nominal_end = %s, \
proc_date = %s, last_error = %s \
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
                                   datetime.now(), '',
                                   self.exper, self.band))

        self.connw.commit()

    def set_error_msg(self, msg):
        """
        Put error comment into DB.

        """
        query = 'UPDATE pima_experiments SET last_error = %s \
WHERE exper_name = %s AND band = %s'

        if not self.connected:
            self._connect()

        with self.connw.cursor() as cursor:
            cursor.execute(query, (msg, self.exper, self.band))

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

        #  Select directory for raw data from a correlator
        self.data_dir = os.getenv('PYPIMA_DATA_DIR')
        if not self.data_dir:
            home_dir = os.getenv('HOME')
            self.data_dir = os.path.join(home_dir, 'data/VLBI/RA/ASC_results')

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
        self.calibration_loaded = False

        # PIMA control file path
        self.cnt_file_name = os.path.join(self.work_dir, '{}_{}_pima.cnt'.
                                          format(self.exper, self.band))

        # Prepare data paths
        self.uv_fits = uv_fits
        self.orbit = orbit

        # Connect to DB
        self.db = DB(self.exper, self.band)

        # Create PIMA control file
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
        data_dir = os.path.join(self.data_dir, self.exper)
        fits_url, size = self.db.get_uvfits_url()
        if not fits_url:
            self._error('Could not find FITS file name in DB')

        # Delete spaces in filename
        self.uv_fits = os.path.join(data_dir, os.path.basename(fits_url).
                                    replace(' ', ''))
        lock_file_name = self.uv_fits + '.lock'

        if os.path.isfile(self.uv_fits) and \
                os.path.getsize(self.uv_fits) == size:
            self._print_info('file {} already exists'.format(self.uv_fits))
        elif os.path.isfile(lock_file_name):
            self._print_info('File {} is being downloaded now, wait'.format(
                             self.uv_fits))
            while os.path.isfile(lock_file_name):
                print('.', end='')
                sys.stdout.flush()
                time.sleep(10)
            print('')
        else:
            if not os.path.isdir(data_dir):
                os.makedirs(data_dir)
            self._print_info('Start downloading file {}...'.format(fits_url))
            lock_file = open(lock_file_name, 'w')
            lock_file.close()
            try:
                self.uv_fits, _ = urlreq.urlretrieve(fits_url,
                                                     filename=self.uv_fits)
            except URLError as ex:
                self._error('Could not download file {}: {}'.format(
                    fits_url, ex.reason))
            finally:
                os.remove(lock_file_name)

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
                orb_data = orb_request.read().decode()

            orb_data = orb_data.replace('\r\n', '\n').split('\n')
            with open(self.orbit, 'w') as orb_file:
                if not orb_data[0].startswith('CCSDS_OEM_VERS'):
                    orb_file.write('CCSDS_OEM_VERS = 2.0\n')

                # Fix meta information
                for line in orb_data:
                    if line.startswith('CENTER_NAME'):
                        line = 'CENTER_NAME   = Earth Barycenter'
                    elif line.startswith('OBJECT_NAME'):
                        line = 'OBJECT_NAME   = RADIO-ASTRON'
                    elif line.startswith('CREATION'):
                        line = line.replace('CREATION DATE', 'CREATION_DATE')
                    elif line.startswith('STOP_TIME') and len(line) < 20:
                        for back_line in reversed(orb_data):
                            cols = back_line.split()
                            if len(cols):
                                break
                        line = line.strip() + ' ' + cols[0]

                    orb_file.write(line + '\n')

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
        sys.stdout.flush()

    def _print_warn(self, msg):
        """Print warning"""
        print('Warning: {}({}): {}'.format(self.exper, self.band, msg))
        sys.stdout.flush()

    def _error(self, msg):
        """Raise pima.Error exception"""
        raise Error(self.exper, self.band, msg)

    def _fix_antab(self):
        """
        Fix antab.

        """
        if not os.path.isfile(self.antab):
            return

        new_antab = os.path.join(self.work_dir, os.path.basename(self.antab))

        freq_setup = self.pima.frequencies()
        if len(freq_setup) != 2:
            self._error('Expected 2 IFs but get {}'.format(len(freq_setup)))

        # Should we fix frequency setup?
        fix_freq = False
        if freq_setup[0]['side_band'] != -1:
            fix_freq = True

        sta_list = self.pima.station_list(ivs_name=False)

        with open(self.antab, 'r') as inp, open(new_antab, 'w') as out:
            magic = inp.readline()
            if not magic.startswith('! Produced by: TSM'):
                self._print_warn('antab file {} does NOT have magic in the \
first line'.format(self.antab))
                self.antab = None
                return

            # Do not forget to write a 'magic' line to the output file
            out.write(magic)

            for line in inp:
                toks = line.split()
                if len(toks) == 0:
                    continue

                if fix_freq and len(toks) > 9 and toks[1].isdigit():
                    if toks[6] == 'L':
                        toks[6] = 'U'
                        toks[9] = '{:.2f}MHz'.format(
                            freq_setup[0]['freq'] * 1e-6)

                if toks[0] == 'TSYS' and len(toks) > 4:
                    toks[4] = toks[4].upper()
                    if toks[4] not in sta_list:
                        toks.insert(0, '!')

                out.write(' '.join(toks) + '\n')

        self.antab = new_antab

    def load(self, download_only=False):
        """
        Download data, run pima load, and do some checks.

        """
        os.chdir(self.work_dir)

        if self.uv_fits is None:
            self._download_fits()

        if download_only:
            return

        if self.orbit is None:
            self._get_orbit()

        # Always download antab-file.
        self._get_antab()

        # Only one sideband in P-band
        if self.band == 'p':
            self.pima.update_cnt({'END_FRQ:': '1'})

        # Set maximum scan length to 1200 s
        self.pima.update_cnt({'MAX_SCAN_LEN:': '1200.0',
                              'SCAN_LEN_USED:': '1200.0'})

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
                             str(desel_nam))

        # Save memory by reducing oversampling
        if self.pima.ap_minmax()[0] < 0.1:
            self.pima.update_cnt({'FRIB.OVERSAMPLE_MD:': '2',
                                  'FRIB.OVERSAMPLE_RT:': '2'})

        self._fix_antab()

        # Try to load calibration information from ANTAB
        if os.path.isfile(self.antab):
            try:
                self.pima.load_gains(self.antab)
                self.pima.load_tsys(self.antab)
                self.calibration_loaded = True
            except pypima.pima.Error as err:
                print(err)
                self._print_warn('Could not load calibration information')
                self.calibration_loaded = False

    def _select_ref_sta(self, fri_file):
        """
        Select reference station for bandpass calibration.
        """
        snr_detecton = 5.6
        fri = Fri(fri_file)
        self.sta_ref = None
        snr = 0
        obs = fri.max_snr('RADIO-AS')

        if len(obs):
            if obs['SNR'] < snr_detecton:
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
            if obs['SNR'] < snr_detecton:
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
                                  'BPS.SNR_MIN_FINE:': str(snr),
                                  'FRIB.SNR_DETECTION:': str(snr_detecton)})
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
        else:
            self.pima.update_cnt({'FRIB.FINE_SEARCH:': 'LSQ',
                                  'PHASE_ACCEL_MIN:': '0',
                                  'PHASE_ACCEL_MAX:': '0'})

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
        Do SPLIT and copy UV-FITS files to the final destination.

        Parameters
        ----------
        source : string, optional
            Do split only for given source. By default split all sources in
            the experiment.

        """
        if self.pima.chan_number() > 512:
            self._print_warn('Too many spectral channels to do SPLIT: {}'.
                             format(self.pima.chan_number()))
            return

        if not self.calibration_loaded:
            self._print_warn('Could not do split due to absence of \
calibartion information')
            return

        # Average all spectral channels in each IF
        self.pima.update_cnt({'SPLT.FRQ_MSEG:': str(self.pima.chan_number())})

        if source:
            self.pima.update_cnt({'SPLT.SOU_NAME:': source})
        else:
            self.pima.update_cnt({'SPLT.SOU_NAME:': 'ALL'})

        exper_dir = self.pima.cnt_params['EXPER_DIR:']
        sess_code = self.pima.cnt_params['SESS_CODE:']
        output_dir = os.path.join(exper_dir, sess_code + '_uvs')

        # Delete old uv-fits remained from previous run
        if os.path.isdir(output_dir):
            old_uvfits = glob.glob(output_dir + '/*')
            for fil in old_uvfits:
                os.remove(fil)

        self.pima.split()

    def fringes2db(self):
        """
        Put fringe fitting information to the DB

        """
        fri_file = self.pima.cnt_params['FRINGE_FILE:']
        if not os.path.isfile(fri_file):
            return

        self.db.fri2db(Fri(fri_file), self.pima.exper_info)

    def delete_uvfits(self):
        """Delete UV-FITS file"""
        if os.path.isfile(self.uv_fits):
            os.remove(self.uv_fits)
