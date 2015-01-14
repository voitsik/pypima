# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 19:52:28 2014

@author: Petr Voytsik
"""

from datetime import datetime
import math

ED = 2. * 6371e5  # Earth diameter
C = 2.99792458e10  # Speed of light


# PIMA fri-file obs fields
#  0 -> Obs #
#  1 -> Sca #
#  2 -> Time code %j-%H%M
#  3 -> Source name
#  4 -> Station 1
#  5 -> Station 2
#  7 -> SNR
# 11 -> Scan start time
# 12 -> Scan stop time
# 15 -> Ampl_lsq
# 23 -> Gr_del_lsq value
# 31 -> Ph_rat_lsq value
# 37 -> Ph_acc value
# 79 -> Duration of scan
# 84 -> U
# 85 -> V
# 94 -> Reference frequency
class Fri(object):
    """
    This class represents PIMA fringe file.

    """
    ver100 = '# PIMA Fringe results  v  1.00  Format version of 2010.04.05'
    ver101 = '# PIMA Fringe results  v  1.01  Format version of 2014.02.08'
    ver102 = '# PIMA Fringe results  v  1.02  Format version of 2014.12.24'

    def __init__(self, file_name=None):
        self.records = []
        self.index = 0

        if file_name:
            self.parse_file(file_name)

    def parse_file(self, file_name):
        """
        Parse PIMA fri-file.

        """
        started = False
        header = {}
        self.records.clear()

        with open(file_name, 'r') as fil:
            line = fil.readline()
            if not (line.startswith(self.ver100) or
                    line.startswith(self.ver101) or
                    line.startswith(self.ver102)):
                raise Exception('{} is not PIMA fri-file'.format(file_name))
            for line in fil:
                if line.startswith('# PIMA_FRINGE started'):
                    started = True
                    header.clear()
                    continue
                if line.startswith('# PIMA_FRINGE ended'):
                    started = False
                    continue
                if not started:
                    continue
                toks = line.split()
                if len(toks) < 3:
                    continue
                if toks[0] == '#':
                    if toks[1] == 'Control':
                        header['cnt_file'] = toks[-1]
                    elif toks[1] == 'Experiment':
                        header['exper_code'] = toks[-1]
                    elif toks[1] == 'Session':
                        header['session_code'] = toks[-1]
                    elif toks[1] == 'FRINGE_FITTING_STYLE:':
                        header['FRINGE_FITTING_STYLE'] = toks[-1]
                    elif toks[1] == 'FRIB.POLAR:':
                        header['polar'] = toks[2]
                    elif toks[1] == 'FRIB.FINE_SEARCH:':
                        header['FRIB.FINE_SEARCH'] = toks[-1]
                    elif toks[1] == 'PHASE_ACCELERATION:':
                        header['accel'] = float(toks[2].replace('D', 'e'))
                elif toks[6] != 'FAILURE':
                    self.records.append({})
                    self.records[-1].update(header)
                    self.records[-1]['obs'] = int(toks[0])
                    self.records[-1]['scan'] = int(toks[1])
                    self.records[-1]['time_code'] = toks[2]
                    self.records[-1]['source'] = toks[3]
                    self.records[-1]['sta1'] = toks[4]
                    self.records[-1]['sta2'] = toks[5]
                    try:
                        self.records[-1]['SNR'] = float(toks[7])
                    except ValueError:
                        self.records[-1]['SNR'] = 0
                    self.records[-1]['start_time'] = datetime.strptime(
                        toks[11], '%Y.%m.%d-%H:%M:%S.%f,')
                    self.records[-1]['stop_time'] = datetime.strptime(
                        toks[12], '%Y.%m.%d-%H:%M:%S.%f')
                    self.records[-1]['ampl_lsq'] = \
                        float(toks[15].replace('D', 'e'))
                    self.records[-1]['delay'] = \
                        float(toks[23].replace('D', 'e'))
                    self.records[-1]['rate'] = \
                        float(toks[31].replace('D', 'e'))
                    self.records[-1]['ph_acc'] = \
                        float(toks[37].replace('D', 'e'))
                    self.records[-1]['duration'] = \
                        float(toks[79].replace('D', 'e'))
                    self.records[-1]['U'] = \
                        float(toks[84].replace('D', 'e'))
                    self.records[-1]['V'] = \
                        float(toks[85].replace('D', 'e'))
                    self.records[-1]['ref_freq'] = \
                        float(toks[94].replace('D', 'e'))
                    # Calculated parameters
                    # UV-radius in lambda
                    self.records[-1]['uv_rad'] = math.hypot(
                        self.records[-1]['U'], self.records[-1]['V'])
                    # Position angle
                    self.records[-1]['PA'] = math.degrees(math.atan2(
                        self.records[-1]['U'], self.records[-1]['V']))
                    wave_len = C / self.records[-1]['ref_freq']
                    # UV-radius in Earth diameters
                    self.records[-1]['uv_rad_ed'] = \
                        self.records[-1]['uv_rad'] * wave_len / ED

    def max_snr(self, station=None):
        """
        Return observation record with maximum SNR.

        Parameters
        ----------
        station : str, optional
            If station name is provided select observations with this station
            only.

        Returns
        -------
        result : dict
            Returns fri-record -- dictionary with parameters of the selected
            observation.

        """
        result = {}

        # Select observations with 'station'
        if station:
            records = [rec for rec in self.records if station in (rec['sta1'],
                                                                  rec['sta2'])]
        else:
            records = self.records

        # Sort records by SNR
        records = sorted(records, key=lambda rec: rec['SNR'], reverse=True)
        if records:
            result = records[0]

        return result

    def append(self, rec):
        """
        Append fri-file record to the end of the list.

        """
        self.records.append(rec)

    def __str__(self):
        out = "#Obs Timecode   Source      Sta1/Sta2         SNR    Delay    \
Rate      Accel      Base   Base\n"
        for rec in self.records:
            accel = rec['accel']
            if rec['FRIB.FINE_SEARCH'] == 'ACC':
                accel = rec['ph_acc']

            line = '{:>3}{:>10} {:>8} {:>8}/{:>8} {:8.2f} \
{:8.3f} {:10.3e} {:9.2e} {:8.2f} {:5.1f}\n'.format(rec['obs'],
                                                   rec['time_code'],
                                                   rec['source'],
                                                   rec['sta1'],
                                                   rec['sta2'],
                                                   rec['SNR'],
                                                   1e6 * rec['delay'],
                                                   rec['rate'],
                                                   accel,
                                                   1e-6 * rec['uv_rad'],
                                                   rec['uv_rad_ed'])
            out = out + line

        return out

    def __iter__(self):
        return self

    def __getitem__(self, ind):
        return self.records[ind]

    def __next__(self):
        if self.index == len(self.records):
            raise StopIteration
        self.index = self.index + 1
        return self.records[self.index - 1]

    def __len__(self):
        return len(self.records)
