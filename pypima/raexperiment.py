# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 04:02:35 2013

@author: Petr Voytsik
"""

from datetime import datetime
import glob
from io import BytesIO
import os.path
import pycurl
import shutil
# import sys
import threading
# import time


import pypima.pima
from pypima.fri import Fri
from pypima.pima import Pima


class Error(Exception):
    """Raised when RaExperiment error occurs"""
    def __init__(self, exper, band, msg):
        self.exper = exper
        self.band = band
        self.msg = msg
        self.time = str(datetime.now())

    def __str__(self):
        return '{}({}): {}'.format(self.exper, self.band, self.msg)


class RaExperiment:
    """This class describe experiment in RadioAstron AGN survey"""

    def __init__(self, experiment_code, band, data_base, data_dir=None,
                 uv_fits=None, orbit=None):
        """
        Parameters
        ----------

        experiment_code : str
            Experiment code.
        band : srt
            One letter frequency bad code.
        data_base : pypima.db.DB
            pypima.db.DB instance.
        data_dir : str, optional
            Directory from FITS-IDI. If `None` working directory of the current
            experiment is used.
        uv_fits : str
            Path to the data file (FITS-IDI). If `None` get a file name from
            data base and download file from the FTP archive.
        orbit : str
            Path to a file with reconstracted orbit. If `None` download it from
            the FTP archive.

        """
        # First, set common variables
        self.exper = experiment_code.lower()
        self.band = band.lower()
        self.db = data_base
        self.sta_ref = 'RADIO-AS'
        self.run_id = 0  # Record id in pima_runs database table
        self.lock = threading.Lock()  # Lock for FITS file downloading control

        if self.band not in ('p', 'l', 'c', 'k'):
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

        work_dir_link = os.path.join(self.exp_dir, self.exper)
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
        self.antab = None
        self.antab_downloaded = False
        self.calibration_loaded = False

        # PIMA control file path
        self.cnt_file_name = os.path.join(self.work_dir, '{}_{}_pima.cnt'.
                                          format(self.exper, self.band))

        # Prepare data paths
        self.uv_fits = uv_fits
        self.orbit = orbit

        #  Select directory for raw data from a correlator
        if data_dir:
            self.data_dir = data_dir
        else:
            self.data_dir = self.work_dir

        # Create PIMA control file
        self._mk_cnt()

        self.pima = Pima(self.exper, self.band, self.work_dir)

        #
        self.split_time_aver = 0

    def _print_info(self, msg):
        """Print some information"""
        now = str(datetime.now())
        print(now, 'Info: {}({}): {}'.format(self.exper, self.band, msg),
              flush=True)

    def _print_warn(self, msg):
        """Print warning"""
        now = str(datetime.now())
        print(now, 'Warning: {}({}): {}'.format(self.exper, self.band, msg),
              flush=True)

    def _error(self, msg):
        """Raise pima.Error exception"""
        raise Error(self.exper, self.band, msg)

    def _mk_cnt(self):
        """
        Make new cnt-file from template.

        """
        cnt_templ_name = os.path.join(self.pima_dir, 'share', 'pima',
                                      'TEMPLATE_pima.cnt')

        cnt_templ = open(cnt_templ_name, 'r')
        cnt_file = open(self.cnt_file_name, 'w')

        sess_code = '{}_{}'.format(self.exper, self.band)
        fringe_file = os.path.join(self.work_dir, sess_code + '.fri')
        frires_file = os.path.join(self.work_dir, sess_code + '.frr')

        if self.band in ('l', 'p'):
            polar = 'RR'
        else:
            polar = 'LL'

        for line in cnt_templ:
            if '@CDATE@' in line:
                line = line.replace('@CDATE@', str(datetime.now()))
            elif line.startswith('SESS_CODE:'):
                line = line.replace('@sess_code@', sess_code)
            elif line.startswith('BAND:'):
                line = line.replace('@band@', self.band.upper())
            elif line.startswith('EXPER_DIR:'):
                line = line.replace('@exper_dir@', self.pima_scr)
            elif line.startswith('UV_FITS:') and self.uv_fits:
                line = line.replace('@uv_fits@', self.uv_fits)
            elif line.startswith('FRINGE_FILE:'):
                line = line.replace('@fringe_file@', fringe_file)
            elif line.startswith('FRIRES_FILE:'):
                line = line.replace('@frires_file@', frires_file)
            elif line.startswith('STA_REF:') and self.sta_ref:
                line = '{:<20}{}\n'.format('STA_REF:', self.sta_ref)
            elif line.startswith('EPHEMERIDES_FILE:') and self.orbit:
                line = line.replace('@ephemerides_file@', self.orbit)
            elif line.startswith('POLAR:') or line.startswith('SPLT.POLAR:'):
                line = line.replace('@polar@', polar)
            cnt_file.write(line)

        cnt_templ.close()
        cnt_file.close()

    def _download_fits(self, gvlbi=False):
        """
        Download FITS-file from remote place

        """
        data_dir = os.path.join(self.data_dir, self.exper)
        fits_url, size = self.db.get_uvfits_url(self.exper, self.band, gvlbi)

        if not fits_url:
            self._error('Could not find FITS file name in DB')

        # Delete spaces in filename
        uv_fits = os.path.join(data_dir, os.path.basename(fits_url).
                               replace(' ', ''))
        lock_file_name = uv_fits + '.lock'

        if os.path.isfile(uv_fits) and os.path.getsize(uv_fits) == size:
            self._print_info('File {} already exists'.format(uv_fits))
        elif os.path.isfile(lock_file_name):
            self._error('Lock! This never should happen!!!')
#            self._print_info('File {} is being downloaded now, wait'.
#                             format(uv_fits))
#            while os.path.isfile(lock_file_name):
#                print('.', end='', flush=True)
#                time.sleep(10)
#            print('')
        else:
            if not os.path.isdir(data_dir):
                os.makedirs(data_dir)
            self._print_info('Start downloading file {}...'.format(fits_url))
            lock_file = open(lock_file_name, 'w')
            lock_file.close()
            try:
                with open(uv_fits, 'wb') as fil:
                    _download_it(fits_url, fil, 2)
            except pycurl.error as err:
                self._error('Could not download file {}: {}'.
                            format(fits_url, err))
            finally:
                os.remove(lock_file_name)

            self._print_info('Downloading is complete')

        self.pima.update_cnt({'UV_FITS:': uv_fits})

        # We use self.uv_fits as a flag of FITS file existence, so set it at
        # the end of this function
        self.uv_fits = uv_fits

    def _get_orbit(self):
        """Download reconstructed orbit file from FTP"""
        orbit_url = self.db.get_orbit_url(self.exper)

        if not orbit_url:
            self._error('Could not find reconstructed orbit')

        self.orbit = os.path.join(self.work_dir, os.path.basename(orbit_url))

        self._print_info('Start downloading orbit file {} ...'.
                         format(orbit_url))

        buffer = BytesIO()
        try:
            _download_it(orbit_url, buffer)
        except pycurl.error as err:
            self._error('Could not download file {}: {}'.
                        format(orbit_url, err))

        orb_data = buffer.getvalue().decode().replace('\r\n', '\n').split('\n')

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

        self._print_info('Downloading is complete')
        self.pima.update_cnt({'EPHEMERIDES_FILE:': self.orbit})

    def _get_antab(self):
        """
        Download antab-file from FTP server.

        """
        antab_url = self.db.get_antab_url(self.exper, self.band)

        if not antab_url:
            self._print_warn('Could not get antab file url from DB.')
        else:
            antab_file = os.path.join(self.work_dir, 'antab',
                                      os.path.basename(antab_url) + '.orig')
            self._print_info('Start downloading file {}'.format(antab_url))
            try:
                with open(antab_file, 'wb') as fil:
                    _download_it(antab_url, fil)

                self.antab_downloaded = True
                self._print_info('Downloading is complete.')
                self.antab = self._fix_antab(antab_file)
            except pycurl.error as err:
                self.antab_downloaded = False
                self._print_warn('Could not download file {}: {}'.
                                 format(antab_url, err))

    def _fix_antab(self, antab):
        """
        Fix antab.

        """
        if not antab or not os.path.isfile(antab):
            return None

        new_antab = antab.replace('.orig', '')

        # ANTAB file already exists and prepared
        if antab == new_antab:
            return new_antab

        freq_setup = self.pima.frequencies()

        # Should we fix frequency setup?
        fix_freq = False
        if self.band != 'p' and freq_setup[0]['side_band'] != -1:
            fix_freq = True

        sta_list = self.pima.station_list(ivs_name=False)

        with open(antab, 'r') as inp, open(new_antab, 'w') as out:
            magic = inp.readline()
            if not magic.startswith('! Produced by: TSM'):
                self._print_warn('antab file {} does NOT have magic in the \
first line'.format(antab))
                return None

            # Do not forget to write a 'magic' line to the output file
            out.write(magic)

            for line in inp:
                line = line.strip()

                if line.startswith('POLY') and line.endswith('/'):
                    line = line.replace('/', ' /')
                elif line.startswith('GAIN WB'):
                    line = 'GAIN WB EQUAT DPFU=1.0,1.0 FREQ=1000,5000'
                elif line.startswith('DPFU=1.0'):
                    line = ''
                elif line.startswith('/') and len(line) > 1:
                    line = line.replace('/', '/ ', 1)
                elif not line.startswith('!') and '!' in line:
                    line = line.split('!')[0].strip()

                # VLA (raes11a and friends)
                if 'YY' in sta_list and 'Y27' in line:
                    line = line.replace('Y27', 'YY')

                toks = line.split()

                if len(toks) == 0:
                    continue

                # Fix EF C-band channels table
                if len(toks) == 10 and toks[0] == '!' and toks[1].isdigit():
                    toks.insert(2, "6cm")

                if fix_freq and len(toks) > 9 and toks[1].isdigit():
                    if toks[6] == 'L':
                        toks[6] = 'U'
                        toks[9] = '{:.2f}MHz'.format(
                            freq_setup[0]['freq'] * 1e-6)

                # Deselect stations
                if toks[0] == 'TSYS' and len(toks) > 4:
                    toks[4] = toks[4].upper()
                    if toks[4] not in sta_list:
                        toks.insert(0, '!')
                # EF, L-band GAINs
                elif toks[0] == 'GAIN' and toks[-1] == '/':
                    for ind in range(len(toks)):
                        if toks[ind].startswith('POLY'):
                            toks[ind] = "\n" + toks[ind]
                            break
                elif toks[0] == '/' and len(toks) > 1:
                    toks[1] = "\n" + toks[1]

                out.write(' '.join(toks) + '\n')

        return new_antab

    def load(self, download_only=False, update_db=False,
             scan_length=1200, scan_part=1, gvlbi=False):
        """
        Download data, run pima load, and do some checks.

        Parameters
        ----------
        download_only : bool, optional
            If True, download FITS-file and return.
        update_db : bool, optional
            If True, update database with experiment information.
        scan_length : float, optional
            Set maximum length of scan. Default is 20 min.
        scan_part : int
            1 is full scan, 2 is half of scan. In general `scan_part` can be
            used as run index.
        gvlbi : bool
            Process ground only (GVLBI) part of the experiment

        """
        os.chdir(self.work_dir)

        # If self.uv_fits is not None assume FITS file already exists
        with self.lock:
            if self.uv_fits is None:
                self._download_fits(gvlbi)

        if download_only:
            return

        if self.orbit is None:
            self._get_orbit()

        # Only one sideband at P-band
        if self.band == 'p':
            self.pima.update_cnt({'END_FRQ:': '1'})

        # Set maximum scan length
        self._print_info('Set maximum scan length to {} s'.format(scan_length))
        self.pima.update_cnt({'MAX_SCAN_LEN:': str(scan_length),
                              'SCAN_LEN_USED:': str(scan_length)})

        staging_dir = os.getenv('PYPIMA_STAGING_DIR', default='NO')
        self.pima.update_cnt({'STAGING_DIR:': staging_dir})

        if update_db:
            self.run_id = self.db.add_exper_info(self.exper, self.band,
                                                 os.path.basename(self.uv_fits),
                                                 scan_part)
        self.pima.load()
        if update_db:
            self.db.update_exper_info(self.pima.exper_info, self.run_id)
            if scan_part == 1:
                self.db.model2db(self.run_id, self.pima.clock_model())

        #
        # Various checks and setups
        #
        if self.pima.obs_number() == 0:
            self._error('ZERO observations have been loaded')

        sou_dist = self.pima.source_dist()
        for source, distance in sou_dist.items():
            if distance > 1.0:
                self._error('Dist = {} arcsec for source {}'.
                            format(distance, source))

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

        # Average all spectral channels in each IF when splitting.
        self.pima.update_cnt({'SPLT.FRQ_MSEG:': str(self.pima.chan_number())})

        # Always download antab-file.
        if not self.antab_downloaded:
            self._get_antab()

        # Try to load calibration information from ANTAB
        if self.antab and os.path.isfile(self.antab):
            try:
                self.pima.load_gains(self.antab)
                self.pima.load_tsys(self.antab)
                self.calibration_loaded = True
            except pypima.pima.Error as err:
                print(err)
                self._print_warn('Could not load calibration information')
                self.calibration_loaded = False

    def _select_ref_sta(self, fri):
        """
        Select reference station for bandpass calibration.

        Parameters
        ----------
        fri : pypima.fri.Fri

        """
        snr_detecton = 5.6
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
            self._print_info('new reference station is {}'.
                             format(self.sta_ref))
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

        accel : boot, optional
            If True turn on phase acceleration fitting.

        Returns
        -------
        fri_file : str
            Name of the PIMA fri-file.

        """
        if accel:
            self.pima.update_cnt({'FRIB.FINE_SEARCH:': 'ACC',
                                  'PHASE_ACCEL_MIN:': '-1.D-14',
                                  'PHASE_ACCEL_MAX:': '1.D-14'})
        else:
            self.pima.update_cnt({'FRIB.FINE_SEARCH:': 'LSQ',
                                  'PHASE_ACCEL_MIN:': '0',
                                  'PHASE_ACCEL_MAX:': '0'})

        if bandpass and self.pima.chan_number() > 512:
            self._print_warn('Too many spectral channels for bandpass: {}'.
                             format(self.pima.chan_number()))
            bandpass = False

        if bandpass:
            fri_file = self.pima.coarse()
            fri = Fri(fri_file)
            if not fri:
                self._error('fri-file is empty after coarse.')

            # Now auto select reference station
            if self._select_ref_sta(fri):
                try:
                    self.pima.bpas()
                except pypima.pima.Error as err:
                    print(err)
                    self._print_info('Try INIT bandpass')
                    try:
                        self.pima.bpas(['BPS.MODE:', 'INIT'])
                    except pypima.pima.Error as err:
                        print(err)
                        self._print_info('Continue without bandpass')
                        self.pima.update_cnt({'BANDPASS_FILE:': 'NO'})

            # Set reasonable SNR detection limit
            self.pima.update_cnt({'FRIB.SNR_DETECTION:': 5.7})

        return self.pima.fine()

    def split(self, source=None, average=0):
        """
        Do SPLIT.

        Parameters
        ----------
        source : string, optional
            Do split only for given source. By default split all sources in
            the experiment.
        average : float, optional
            Number of seconds to average data when splitting. `average` <= 0
            disables time averaging.

        """
        if not self.calibration_loaded:
            self._print_warn('Could not do split due to absence of \
calibartion information')
            return

        split_params = []
        time_segments = 1

        if source:
            split_params.extend(('SPLT.SOU_NAME:', source))
        else:
            split_params.extend(('SPLT.SOU_NAME:', 'ALL'))

        ap = self.pima.ap_minmax()[0]
        if average > 0:
            time_segments = round(average / ap)

        self.split_time_aver = time_segments * ap

        exper_dir = self.pima.cnt_params['EXPER_DIR:']
        sess_code = self.pima.cnt_params['SESS_CODE:']
        pima_fits_dir = os.path.join(exper_dir, sess_code + '_uvs')

        # Delete old uv-fits remained from previous run
        if os.path.isdir(pima_fits_dir):
            old_uvfits = glob.glob(pima_fits_dir + '/*')
            for fil in old_uvfits:
                os.remove(fil)

        self.pima.split(tim_mseg=time_segments, params=split_params)

    def copy_uvfits(self, out_dir):
        """
        Copy splited uv-fits from pima scratch dir to out_dir

        """
        exper_dir = self.pima.cnt_params['EXPER_DIR:']
        sess_code = self.pima.cnt_params['SESS_CODE:']
        band = self.pima.cnt_params['BAND:']
        polar = self.pima.cnt_params['POLAR:']

        pima_fits_dir = os.path.join(exper_dir, sess_code + '_uvs')

        sources = self.pima.source_list()
        splt_sou_name = self.pima.cnt_params['SPLT.SOU_NAME:']

        for source_names in sources:
            if splt_sou_name != 'ALL' and splt_sou_name not in source_names:
                continue

            pima_fits_name = '{}_{}_uva.fits'.format(source_names[1], band)
            pima_fits_path = os.path.join(pima_fits_dir, pima_fits_name)

            if not os.path.isfile(pima_fits_path):
                self._print_warn('UV-FITS "{}" does not exists.'.format(
                    pima_fits_path))
                continue

            # Use B1950 name for output directory
            out_fits_dir = os.path.join(out_dir, source_names[2])
            if not os.path.isdir(out_fits_dir):
                os.mkdir(out_fits_dir)

            out_fits_name = \
                '{}_{}_{}_{}_{:04d}s_uva.fits'.format(source_names[2],
                                                      self.exper,
                                                      self.band.upper(),
                                                      polar,
                                                      round(self.split_time_aver))
            out_fits_path = os.path.join(out_fits_dir, out_fits_name)

            self._print_info('Copy {} to {}'.format(pima_fits_path,
                                                    out_fits_path))
            shutil.copy(pima_fits_path, out_fits_path)

            # Run `fits_to_radplot`
            pypima.pima.fits_to_txt(out_fits_path)

    def fringes2db(self):
        """
        Put fringe fitting information to the database.

        """
        if self.run_id <= 0:
            return

        fri_file = self.pima.cnt_params['FRINGE_FILE:']
        if not os.path.isfile(fri_file):
            return

        self.db.fri2db(Fri(fri_file), self.pima.exper_info, self.run_id)

    def delete_uvfits(self):
        """
        Delete UV-FITS file.

        """
        if os.path.isfile(self.uv_fits):
            os.remove(self.uv_fits)

    def generate_autospectra(self, out_dir):
        """
        """
        # Sometimes PIMA crashes on `acta` task
        try:
            file_list = self.pima.acta()
        except pypima.pima.Error as err:
            print('PIMA Error: ', err)
            return

        polar = self.pima.cnt_params['POLAR:']
        plot_dir = '{}_{}_{}'.format(self.exper, self.band, polar)
        tmp_plot_dir = os.path.join('/tmp', plot_dir)

        if os.path.exists(tmp_plot_dir):
            shutil.rmtree(tmp_plot_dir)

        os.mkdir(tmp_plot_dir)

        for file_name in file_list:
            out_plot_name = os.path.basename(file_name).replace('.txt', '.gif')
            tmp_plot_path = os.path.join(tmp_plot_dir, out_plot_name)
            pypima.pima.acta_plot(file_name, tmp_plot_path)

        # Check destination
        final_dir = os.path.join(out_dir, plot_dir)
        if os.path.exists(final_dir):
            shutil.rmtree(final_dir)

        shutil.move(tmp_plot_dir, out_dir)


def _download_it(url, buffer, max_retries=0):
    """
    Download data from `url` and write it to `buffer` using pycurl.

    Parameters
    ----------
    url : str
        URL
    buffer : object
        Object with `write` function. For inctance, BytesIO or file descriptor.
    max_retries : int
        Number of attempts to download.

    """
    done = False
    retries = 0

    curl = pycurl.Curl()
    curl.setopt(pycurl.URL, url)
    curl.setopt(pycurl.CONNECTTIMEOUT, 30)
    curl.setopt(pycurl.LOW_SPEED_LIMIT, 10000)
    curl.setopt(pycurl.LOW_SPEED_TIME, 60)
    curl.setopt(pycurl.WRITEDATA, buffer)

    while not done:
        try:
            curl.perform()
        except pycurl.error as err:
            errno, errstr = err.args
            if errno == 28 and retries < max_retries:
                retries += 1
                buffer.seek(0)
                buffer.truncate(0)
            else:
                raise
        else:
            done = True

    curl.close()
