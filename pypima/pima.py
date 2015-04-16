# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 17:21:29 2013

@author: Petr Voytsik
"""

from datetime import datetime
import glob
import os.path
import shutil
import subprocess
import sys


class Error(Exception):
    """Raised when PIMA error occurs"""
    def __init__(self, exper, band, msg):
        self.exper = exper
        self.band = band
        self.msg = msg
        self.time = str(datetime.now())

    def __str__(self):
        return '{}({}): {}'.format(self.exper, self.band, self.msg)


class ExperInfo:
    """Some experiment information"""

    def __init__(self, exper_name, band, stt_file=None):
        self.exper = exper_name
        self.band = band

        self.sp_chann_num = None
        self.time_epochs_num = None
        self.scans_num = None
        self.obs_num = None
        self.uv_points_num = None
        self.uv_points_used_num = None
        self.deselected_points_num = None
        self.no_auto_points_num = None
        self.accum_length = None
        self.utc_minus_tai = None
        self.nominal_start = None
        self.nominal_end = None

        if stt_file:
            self.update(stt_file)

    def update(self, stt_file):
        """
        Fill ExperInfo with info from PIMA stt-file
        """

        with open(stt_file, 'r') as fil:
            for line in fil:
                if line.startswith('Number of spectral channels:'):
                    self.sp_chann_num = int(line.split()[4])
                elif line.startswith('Number of time epochs:'):
                    self.time_epochs_num = int(line.split()[4])
                elif line.startswith('Number of scans:'):
                    self.scans_num = int(line.split()[3])
                elif line.startswith('Number of observations:'):
                    self.obs_num = int(line.split()[3])
                elif line.startswith('Total number of UV points:'):
                    self.uv_points_num = int(line.split()[5])
                elif line.startswith('Total number of used UV points:'):
                    self.uv_points_used_num = int(line.split()[6])
                elif line.startswith('Total number of deselected points:'):
                    self.deselected_points_num = int(line.split(':')[1])
                elif line.startswith('Number of cross-correl NO_AUTO_1 \
deselected points:'):
                    no_auto1 = int(line.split(':')[1])
                elif line.startswith('Number of cross-correl NO_AUTO_2 \
deselected points:'):
                    no_auto2 = int(line.split(':')[1])
                elif line.startswith('Accummulation period length_min:'):
                    acc_min = float(line.split(':')[1])
                elif line.startswith('Accummulation period length_max:'):
                    acc_max = float(line.split(':')[1])
                elif line.startswith('UTC_minus_TAI:'):
                    self.utc_minus_tai = float(line.split(':')[1])
                elif line.startswith('Experiment nominal start:'):
                    self.nominal_start = datetime.strptime(
                        line.split(':', 1)[1].strip()[:23],
                        '%Y.%m.%d-%H:%M:%S.%f')
                elif line.startswith('Experiment nominal end:'):
                    self.nominal_end = datetime.strptime(
                        line.split(':', 1)[1].strip()[:23],
                        '%Y.%m.%d-%H:%M:%S.%f')

        self.accum_length = (acc_min + acc_max) / 2.
        self.no_auto_points_num = no_auto1 + no_auto2


class Pima(object):
    """
    Pima class is analog of pima_fringe.csh script.

    """
    def __init__(self, experiment_code, band, work_dir=None):
        # First, set common variables
        self.exper = experiment_code.lower()
        self.band = band.lower()
        self.work_dir = work_dir

        if work_dir is None:
            self.work_dir = os.getcwd()

        self.pima_dir = os.getenv('PIMA_DIR')

        if self.pima_dir is None or not os.path.isdir(self.pima_dir):
            raise Error(self.exper, self.band, 'Could not find PIMA directory.\
 Please, set $PIMA_DIR environment variable.')

        self.pima_exec = os.path.join(self.pima_dir, 'bin', 'pima')

        if not os.path.isfile(self.pima_exec):
            raise Error(self.exper, self.band, 'Could not find pima \
executable. Check your PIMA installation!')

        # PIMA control file path
        self.cnt_file_name = '{}_{}_pima.cnt'.format(self.exper, self.band)
        self.cnt_file_name = os.path.join(self.work_dir, self.cnt_file_name)

        if not os.path.isfile(self.cnt_file_name):
            raise Error(self.exper, self.band, 'Could not find control file \
{}'.format(self.cnt_file_name))

        # Dictionary with all parameters from cnt-file
        self.cnt_params = {}
        self._update_cnt_params()
        self.exper_info = ExperInfo(self.exper, self.band)
        # Update exper_info if experiment already loaded.
        stt_file = os.path.join(self.cnt_params['EXPER_DIR:'],
                                self.cnt_params['SESS_CODE:'] + '.stt')
        if os.path.isfile(stt_file):
            self.exper_info.update(stt_file)

    def _update_cnt_params(self):
        """Read cnt-file and fill cnt_params dictionary."""
        self.cnt_params.clear()
        self.cnt_params['UV_FITS:'] = list()

        with open(self.cnt_file_name, 'r') as cnt_file:
            for line in cnt_file:
                line = line.split('#')[0].strip()
                if len(line) < 8:
                    continue
                key, val = line.split(None, 1)

                # Special case
                if key == 'UV_FITS:':
                    self.cnt_params[key].append(val)
                else:
                    self.cnt_params[key] = val

    def update_cnt(self, opts):
        """
        Update pima cnt-file according 'opts' dictionary.

        """
        if not isinstance(opts, dict):
            return

        uv_fits = False

        old_cnt = open(self.cnt_file_name, 'r')
        new_cnt = open(self.cnt_file_name + '.new', 'w')

        for line in old_cnt:
            key = line.split()[0].strip()

            if key in opts.keys():
                val = opts[key]

                # Skip empty values
                if not val:
                    continue

                # Modify UV_FITS only once
                if key == 'UV_FITS:':
                    if uv_fits:
                        continue
                    else:
                        uv_fits = True

                # In case of many FITS-files
                if key == 'UV_FITS:' and isinstance(val, list):
                    line = ''
                    for items in val:
                        line += '{:<20}{}\n'.format(key, items)
                else:
                    line = '{:<20}{}\n'.format(key, val)
            elif line.startswith('# Last update on'):
                line = '# Last update on  {}\n'.format(str(datetime.now()))

            new_cnt.write(line)

        old_cnt.close()
        new_cnt.close()
        os.rename(self.cnt_file_name + '.new', self.cnt_file_name)
        self._update_cnt_params()

    def _exec(self, operation, options=None, log_name=None):
        """Execute PIMA binary."""
        cmd_line = ['nice', '-n', '19', self.pima_exec, self.cnt_file_name,
                    operation]
        if options:
            cmd_line.extend(options)

        if log_name is None:
            log_name = os.path.join(self.work_dir,
                                    '{}_{}_{}.log'.format(self.exper,
                                                          self.band,
                                                          operation))

        log = open(log_name, 'w')
        log.write(str(datetime.now()) + '\n\n')
        log.flush()

        ret = subprocess.call(cmd_line, stdout=log, universal_newlines=True)

        log.write('\n' + str(datetime.now()) + '\n\n')
        log.close()

        return ret

    def _print_info(self, msg):
        """Print some information"""
        now = str(datetime.now())
        print(now, 'Info: {}({}): {}'.format(self.exper, self.band, msg))
        sys.stdout.flush()

    def _error(self, msg):
        """Raise pima.Error exception"""
        raise Error(self.exper, self.band, msg)

    def load(self):
        """
        Run pima load.

        """
        # Delete existing PIMA auxiliary files
        auxiliary_files = os.path.join(self.cnt_params['EXPER_DIR:'],
                                       self.cnt_params['SESS_CODE:'])
        auxiliary_files = glob.glob(auxiliary_files + '*')

        for aux_file in auxiliary_files:
            if os.path.isfile(aux_file):
                os.remove(aux_file)
            elif os.path.isdir(aux_file):
                shutil.rmtree(aux_file)

        opts = ['BANDPASS_FILE:', 'NO',
                'POLARCAL_FILE:', 'NO']
        ret = self._exec('load', options=opts)
        if ret:
            self._error('load failed with code {}'.format(ret))

        stt_file = os.path.join(self.cnt_params['EXPER_DIR:'],
                                self.cnt_params['SESS_CODE:'] + '.stt')
        self.exper_info.update(stt_file)

        self._print_info('load ok')

    def coarse(self, params=None):
        """
        Do coarse fringe fitting.

        This function runs 'pima frib' with spetial parameters:

        1. Disable bandpass

        2. Use fast algorithm for fringe search

        3. Set fri-file name

        Parameters
        ----------
        params : list, optional
            List of the optional pima parameters. Must have even number of
            the elements

        Returns
        -------
        fri_file : str
            Name of the fri-file

        """
        log_name = '{}_{}_coarse.log'.format(self.exper, self.band)
        log_name = os.path.join(self.work_dir, log_name)
        fri_file = '{}_{}_nobps.fri'.format(self.exper, self.band)
        fri_file = os.path.join(self.work_dir, fri_file)
        frr_file = '{}_{}_nobps.frr'.format(self.exper, self.band)
        frr_file = os.path.join(self.work_dir, frr_file)

        if os.path.isfile(fri_file):
            os.remove(fri_file)

        if os.path.isfile(frr_file):
            os.remove(frr_file)

        opts = ['FRINGE_FILE:', fri_file,
                'FRIRES_FILE:', frr_file,
                'BANDPASS_USE:', 'NO',
                'BANDPASS_FILE:', 'NO',
                'POLARCAL_FILE:', 'NO',
                'FRIB.OVERSAMPLE_MD:', '4',
                'FRIB.OVERSAMPLE_RT:', '4',
                'FRIB.SECONDARY_SNR_MIN:', '0.0',
                'FRIB.SECONDARY_MAX_TRIES:', '0',
                'FRIB.FINE_SEARCH:', 'PAR',
                'MKDB.FRINGE_ALGORITHM:', 'DRF',
                'PHASE_ACCEL_MIN:', '0',
                'PHASE_ACCEL_MAX:', '0']

        if params:
            opts.extend(params)

        ret = self._exec('frib', opts, log_name)

        if ret:
            self._error('coarse failed with code {}'.format(ret))

        self._print_info('coarse ok')

        return fri_file

    def fine(self, params=None):
        """
        Do file fringe fitting.

        This function runs 'pima frib' with parameters from cnt-file plus user
        defined params.

        Parameters
        ----------
        params : list, optional
            List of the optional pima parameters. Must have even number of
            elements.

        Returns
        -------
        fri_file : str
            Name of the fri-file.

        Notes
        -----
        Pima.fine() uses fri-file name from the cnt-file. Set up it before run.

        """
        log_name = '{}_{}_fine.log'.format(self.exper, self.band)
        log_name = os.path.join(self.work_dir, log_name)
        fri_file = self.cnt_params['FRINGE_FILE:']
        frr_file = self.cnt_params['FRIRES_FILE:']

        if os.path.isfile(fri_file):
            os.remove(fri_file)

        if os.path.isfile(frr_file):
            os.remove(frr_file)

        ret = self._exec('frib', params, log_name=log_name)

        if ret:
            self._error('fine failed with code {}'.format(ret))

        self._print_info('fine ok')

        return fri_file

    def bpas(self, params=None):
        """
        Do bandpass calibration.

        This function runs 'pima bpas'.

        Parameters
        ----------
        params : list, optional
            List of the optional pima parameters. Must have even number of
            elements.

        """
        fri_file = '{}_{}_nobps.fri'.format(self.exper, self.band)
        fri_file = os.path.join(self.work_dir, fri_file)
        log_file = '{}_{}_bps.log'.format(self.exper, self.band)
        log_file = os.path.join(self.work_dir, log_file)
        exc_obs_file = '{}_{}_bpas_obs.exc'.format(self.exper, self.band)
        exc_obs_file = os.path.join(self.work_dir, exc_obs_file)

        if self.cnt_params['BANDPASS_FILE:'] == 'NO':
            bps_file = '{}_{}.bps'.format(self.exper, self.band)
            bps_file = os.path.join(self.work_dir, bps_file)
            self.update_cnt({'BANDPASS_FILE:': bps_file})

        opts = ['FRINGE_FILE:', fri_file,
                'DEBUG_LEVEL:', '3']

        if os.path.isfile(exc_obs_file):
            opts.extend(['EXCLUDE_OBS_FILE:', exc_obs_file])

        if self.cnt_params['POLAR:'] in ['I', 'RPL']:
            opts.extend(['POLAR:', 'RR'])

        if params:
            opts.extend(params)

        ret = self._exec('bpas', opts, log_file)
        if ret:
            self._error('bpas failed with code {}'.format(ret))

        self._print_info('bpas ok')

    def split(self, tim_mseg=1, params=None):
        """
        Do SPLIT.

        This function runs 'pima splt' with given integration time.

        Parameters
        ----------
        tim_mseg : int
            Number of time segments to integrate.

        params : list, optional
            List of the optional pima parameters. Must have even number of
            elements.

        """
        opts = ['SPLT.TIM_MSEG:', str(tim_mseg)]
        if params:
            opts.extend(params)
        ret = self._exec('splt', opts)

        if ret:
            self._error('splt failed with code {}'.format(ret))

        self._print_info('split ok')

    def load_gains(self, gain_file, params=None):
        """
        Load antenna gains from the `gain_file`.

        """
        opts = ['evn_gain', gain_file, 'DEBUG_LEVEL:', '6']
        if params:
            opts.extend(params)

        log_file = '{}_{}_gain.log'.format(self.exper, self.band)
        log_file = os.path.join(self.work_dir, log_file)

        ret = self._exec('gean', opts, log_file)

        if ret:
            self._error('evn_gain failed with code {}'.format(ret))

        self._print_info('evn_gain ok')

    def load_tsys(self, tsys_file, params=None):
        """
        Load Tsys from the `tsys_file`.

        """
        opts = ['vlba_log_file', tsys_file, 'DEBUG_LEVEL:', '1']
        if params:
            opts.extend(params)

        log_file = '{}_{}_tsys.log'.format(self.exper, self.band)
        log_file = os.path.join(self.work_dir, log_file)

        ret = self._exec('gean', opts, log_file)

        if ret:
            self._error('vlba_log_file failed with code {}'.format(ret))

        self._print_info('vlba_log_file ok')

    def acta(self, params=None):
        """
        Run `acta` pima task for autocorrelation spectrum generation.

        Returns
        -------
        file_list : list
            List of the names of the generated files.

        """
        opts = ['DEBUG_LEVEL:', '2']
        if params:
            opts.extend(params)

        log_file = '{}_{}_acta.log'.format(self.exper, self.band)
        log_file = os.path.join(self.work_dir, log_file)

        ret = self._exec('acta', opts, log_file)

        if ret:
            self._error('acta failed with code {}'.format(ret))
        else:
            self._print_info('acta ok')

        file_list = []
        with open(log_file, 'r') as fil:
            for line in fil:
                if line.startswith('PIMA_ACTA created file'):
                    file_name = line.split(':')[1].strip()
                    file_list.append(file_name)

        return file_list

    # Additional useful utilites
    def set_polar(self, polar):
        """
        Set polarization in control file.

        """
        polar = polar.upper()

        if polar not in ['RR', 'RL', 'LR', 'LL', 'I']:
            self._error('Wrong polarization: ' + polar)

        self._print_info('Set polarization to ' + polar)
        self.update_cnt(({'POLAR:': polar, 'SPLT.POLAR:': polar}))

    def ap_minmax(self):
        """
        Return minimum and maximum accummulation periods in experiment.

        """
        ap_min = ap_max = 0
        stt_file = os.path.join(self.cnt_params['EXPER_DIR:'],
                                self.cnt_params['SESS_CODE:'] + '.stt')

        if os.path.isfile(stt_file):
            with open(stt_file, 'r') as fil:
                for line in fil:
                    if line.startswith('Accummulation period length_min'):
                        ap_min = float(line.split()[3])
                    elif line.startswith('Accummulation period length_max'):
                        ap_max = float(line.split()[3])

        return ap_min, ap_max

    def number_of_deselected_points(self):
        """
        Return total number of deselected points

        """
        return self.exper_info.deselected_points_num

    def station_list(self, ivs_name=True):
        """
        Return a list of the station names.

        This function returns list of station names which participated in the
        experiment.

        Parameters
        ----------
        ivs_name : bool, optional
            If `ivs_name` is False this function returns list of 2-letter
            station codes instead of IVS names.

        Returns
        -------
        list
            List of the station names.

        """
        sta_l = []
        sta_file = os.path.join(self.cnt_params['EXPER_DIR:'],
                                self.cnt_params['SESS_CODE:'] + '.sta')

        if os.path.isfile(sta_file):
            with open(sta_file, 'r') as fil:
                for line in fil:
                    toks = line.split()
                    if ivs_name:
                        sta_l.append(toks[3])
                    else:
                        sta_l.append(toks[5])

        return sta_l

    def source_list(self):
        """
        Return a list of the source names.

        Returns
        -------
        out : list
            Each item in the out is a tuple of 3 names: IVS, J2000, B1950.

        """
        sou_list = []
        sou_file = os.path.join(self.cnt_params['EXPER_DIR:'],
                                self.cnt_params['SESS_CODE:'] + '.sou')

        if os.path.isfile(sou_file):
            with open(sou_file) as fil:
                for line in fil:
                    toks = line.split()
                    sou_list.append((toks[2], toks[3], toks[4]))

        return sou_list

    def source_dist(self):
        """
        Return distance between correlator phase center and source position
        from catalog.

        """
        dist = {}
        sou_file = os.path.join(self.cnt_params['EXPER_DIR:'],
                                self.cnt_params['SESS_CODE:'] + '.sou')

        if os.path.isfile(sou_file):
            with open(sou_file) as fil:
                for line in fil:
                    toks = line.split()
                    try:
                        dist[toks[2]] = float(toks[12])
                    except ValueError:
                        dist[toks[2]] = 9999999.9

        return dist

    def obs_number(self):
        """
        Return number of the observations in the experiment.

        """
        return self.exper_info.obs_num

    def chan_number(self):
        """
        Return number of the spectral channels in uv-data.

        """
        return self.exper_info.sp_chann_num

    def frequencies(self):
        """
        Return list of frequencies used in the experiment.

        Returns
        -------
        freqs : list
            The function returns a list of dictionaries.

        """
        freqs = []
        frq_file = os.path.join(self.cnt_params['EXPER_DIR:'],
                                self.cnt_params['SESS_CODE:'] + '.frq')

        if os.path.isfile(frq_file):
            with open(frq_file) as fil:
                for line in fil:
                    if len(line) < 16:
                        break
                    if line.startswith('Ind_grp:'):
                        toks = line.split()
                        freq = dict()
                        freq['freq'] = float(toks[6])
                        freq['band_width'] = float(toks[8])
                        freq['chan_width'] = float(toks[10])
                        freq['side_band'] = int(toks[12])
                        freqs.append(freq)

        return freqs

    def clock_model(self):
        """
        Read clock model components from the PIMA mdc file and retrun
        they as list of records.

        Returns
        -------
        clock_model : list


        """
        clock_model = []
        mdc_file = os.path.join(self.cnt_params['EXPER_DIR:'],
                                self.cnt_params['SESS_CODE:'] + '.mdc')

        with open(mdc_file, 'r') as fil:
            for line in fil:
                if not line.startswith('CLOCK_MODEL'):
                    continue

                cols = line.strip().split()
                if len(cols) != 19:
                    continue

                sta = cols[2]
                time = datetime.strptime(cols[6], '%Y.%m.%d-%H:%M:%S.%f')
                clock_offset = float(cols[8].replace('D', 'e'))
                clock_rate = float(cols[10].replace('D', 'e'))
                group_delay = float(cols[16].replace('D', 'e'))
                delay_rate = float(cols[18].replace('D', 'e'))

                clock_model.append((sta, time, clock_offset, clock_rate,
                                    group_delay, delay_rate))

        return clock_model


def fits_to_txt(fits_file):
    """
    Dump visibilites from the uvf file using the fits_to_radplot utility.

    fits_to_radplot reads visibility data from the input file and write
    amplitude, phase, and weight for each baseline to the text file.

    Parameters
    ----------
    fits_file : str
        Name of the input FITS file. FITS file must be in uvf format.

    Returns
    -------
    out : str
        Returns the stdout of the fits_to_radplot utility.

    """

    if not os.path.isfile(fits_file):
        raise IOError(2, 'No such file: {}'.format(fits_file))

    txt_file = fits_file.replace('.fits', '.txt')

    # If fits_file has extention other than .fits
    if txt_file == fits_file:
        txt_file = fits_file + '.txt'

    cmd_line = ['fits_to_radplot', '-o', txt_file, fits_file]
    out = subprocess.check_output(cmd_line, universal_newlines=True)

    return out


def acta_plot(input_file, output_file):
    """
    """
    if not os.path.isfile(input_file):
        raise IOError(2, 'No such file: {}'.format(input_file))

    cmd_line = ['acta_plot', input_file, output_file]
    out = subprocess.check_output(cmd_line, universal_newlines=True)

    return out
