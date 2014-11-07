# -*- coding: utf-8 -*-
"""
Created on Thu Oct 30 14:23:23 2014

@author: Petr Voytsik
"""

from datetime import datetime
from datetime import timedelta
import netrc
import psycopg2


class DB:
    """Interface to ra_results DataBase"""

    def __init__(self):
        nrc = netrc.netrc()
        self.web_login = nrc.authenticators('webinet.asc.rssi.ru')[0]
        self.web_passw = nrc.authenticators('webinet.asc.rssi.ru')[2]
        self.arc_login = nrc.authenticators('archive.asc.rssi.ru')[0]
        self.arc_passw = nrc.authenticators('archive.asc.rssi.ru')[2]

        self.conn = psycopg2.connect(database='ra_results', user='guest',
                                     host='odin')
        self.connw = psycopg2.connect(database='ra_results', user='editor',
                                      host='odin')

    def get_uvfits_url(self, exper, band):
        """Get FITS-file url from DB for given experiment and band"""

        url = None
        size = 0
        url_base = 'ftp://{}:{}@archive.asc.rssi.ru'.format(self.arc_login,
                                                            self.arc_passw)

        with self.conn.cursor() as cursor:
            cursor.execute('SELECT path, size FROM fits_files WHERE \
LOWER(exper_name) = LOWER(%s) AND LOWER(band) = LOWER(%s) \
ORDER BY corr_date DESC, path DESC;',
                           (exper, band))
            reply = cursor.fetchone()

        if reply:
            path = reply[0]
            size = reply[1]
            url = url_base + path

        return url, size

    def get_orbit_url(self, exper):
        """Returns orbit file url for given experiment"""

        url = None
        url_base = 'ftp://{}:{}@webinet.asc.rssi.ru/radioastron/\
oddata/reconstr/'.format(self.web_login, self.web_passw)

        with self.conn.cursor() as cursor:
            cursor.execute("SELECT scf_files.file_name FROM scf_files, \
vex_files WHERE vex_files.exper_name = %s AND \
scf_files.start_time <= vex_files.exper_nominal_start AND \
scf_files.stop_time >= vex_files.exper_nominal_stop;", (exper,))
            reply = cursor.fetchone()

        if reply:
            orbit_file = reply[0]
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
            reply = cursor.fetchone()

        if reply:
            date = reply[0]
            date1 = date[0:7]
            url = '{0}/{1}/{2}_{3}/{3}{4}.antab'.format(url_base, date1, date,
                                                        exper, band)

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

        polar = fri_file[0]['polar']
        exper, band = fri_file[0]['session_code'].split('_', 1)

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

            # Update status of the observations
            query = "UPDATE pima_observations SET status = %s WHERE \
exper_name = %s AND band = %s AND polar = %s AND snr <= %s;"
            cur.execute(query, ('n', exper, band, polar, 5.3))

            query = "UPDATE pima_observations SET status = %s WHERE \
exper_name = %s AND band = %s AND polar = %s AND snr >= %s;"
            if exper_info.sp_chann_num <= 128:
                cur.execute(query, ('y', exper, band, polar, 5.7))
            else:
                cur.execute(query, ('y', exper, band, polar, 7.0))

        self.connw.commit()

    def exper_info2db(self, exper_info, uv_fits):
        """
        Put experiment info to the DB.

        """

        with self.connw.cursor() as cursor:
            cursor.execute("SELECT exper_name, band FROM pima_experiments \
WHERE exper_name = %s AND band = %s", (exper_info.exper, exper_info.band))
            test = cursor.fetchone()
            if test is None:
                cursor.execute("INSERT INTO pima_experiments (exper_name, \
band) VALUES (%s, %s);", (exper_info.exper, exper_info.band))

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
                                   exper_info.exper, exper_info.band))

        self.connw.commit()

    def set_error_msg(self, exper, band, msg):
        """
        Put error comment into DB.

        """
        query = 'UPDATE pima_experiments SET last_error = %s \
WHERE exper_name = %s AND band = %s'

        with self.connw.cursor() as cursor:
            cursor.execute(query, (msg, exper, band))

        self.connw.commit()
