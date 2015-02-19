"""
Created on Thu Oct 30 14:23:23 2014

@author: Petr Voytsik
"""

from datetime import datetime, timedelta
import netrc
import psycopg2


class DB:
    """
    Interface to the "ra_results" database.

    """
    def __init__(self):
        """
        Connect to the database and load logins and passwords from ``.netrc``.

        """
        nrc = netrc.netrc()
        self.web_login = nrc.authenticators('webinet.asc.rssi.ru')[0]
        self.web_passw = nrc.authenticators('webinet.asc.rssi.ru')[2]
        self.arc_login = nrc.authenticators('archive.asc.rssi.ru')[0]
        self.arc_passw = nrc.authenticators('archive.asc.rssi.ru')[2]

        self.conn = psycopg2.connect(database='ra_results', user='guest',
                                     host='odin')
        self.conn.autocommit = True
        self.connw = psycopg2.connect(database='ra_results', user='editor',
                                      host='odin')

    def get_uvfits_url(self, exper, band):
        """
        Get FITS-file URL and size from the database for the given experiment
        and band.

        Parameters
        ----------
        exper : str
            Experiment name.
        band : str
            Frequency band.

        Returns
        -------
        url, size : (str, int)
            Tuple of file URL and size. Returns (None, 0) if the database reply
            is empty.

        Notes
        -----
        If database contains more than one record for given experiment and
        band, this function selects the most recent FITS-file.

        """
        url = None
        size = 0
        url_base = 'ftp://{}:{}@archive.asc.rssi.ru'.format(self.arc_login,
                                                            self.arc_passw)

        # TODO: Check fits-file versions
        query = 'SELECT path, size FROM fits_files WHERE \
LOWER(exper_name) = LOWER(%s) AND LOWER(band) = LOWER(%s) AND path LIKE %s \
ORDER BY corr_date DESC, path DESC;'

        with self.conn.cursor() as cursor:
            cursor.execute(query, (exper, band, '%RADIOASTRON%'))
            reply = cursor.fetchone()

        if reply:
            path = reply[0]
            size = reply[1]
            url = url_base + path

        return url, size

    def get_orbit_url(self, exper):
        """
        Return reconstructed orbit file URL for the given experiment.

        Parameters
        ----------
        exper : str
            Experiment name.

        Returns
        -------
        url : str
            File URL on FTP server. Returns None the database reply is empty.

        """
        url = None
        url_base = 'ftp://{}:{}@webinet.asc.rssi.ru/radioastron/\
oddata/reconstr/'.format(self.web_login, self.web_passw)

        query = 'SELECT scf_files.file_name FROM scf_files, vex_files \
WHERE vex_files.exper_name = %s AND \
scf_files.start_time <= vex_files.exper_nominal_start AND \
scf_files.stop_time >= vex_files.exper_nominal_stop;'

        with self.conn.cursor() as cursor:
            cursor.execute(query, (exper,))
            reply = cursor.fetchone()

        if reply:
            orbit_file = reply[0]
            url = url_base + orbit_file

        return url

    def get_antab_url(self, exper, band):
        """
        Return ANTAB-file URL for the given experiment and band.

        Parameters
        ----------
        exper : str
            Experiment name.
        band : str
            Frequency band.

        Returns
        -------
        url : str
            File URL on FTP server. Returns None the database reply is empty.

        """
        url = None
        url_base = 'ftp://{}:{}@webinet.asc.rssi.ru/radioastron/\
ampcal'.format(self.web_login, self.web_passw)

        query = "SELECT to_char(exper_nominal_start, 'YYYY_MM_DD') \
                 FROM vex_files WHERE exper_name = %s;"

        with self.conn.cursor() as cursor:
            cursor.execute(query, (exper,))
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
            cursor.execute("SELECT * FROM pima_obs WHERE \
exper_name = %s AND band = %s AND polar = %s", (exper, band, polar))
            reply = cursor.fetchone()
            if reply:
                cursor.execute("DELETE FROM pima_obs WHERE \
exper_name = %s AND band = %s AND polar = %s", (exper, band, polar))

    def fri2db(self, fri_file, exper_info, run_id):
        """
        Store information from the PIMA fri-file to the database.

        """
        exper = exper_info.exper
        band = exper_info.band
        polar = fri_file[0]['polar']

#        self._check_and_delete(exper, band, polar)

        query_insert = 'INSERT INTO pima_obs (obs, start_time, \
stop_time, exper_name, band, source, polar, st1, st2, delay, rate, accel, \
snr, ampl, solint, u, v, base_ed, ref_freq, scan_name, run_id) \
VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, \
%s, %s, %s, %s);'

        with self.connw.cursor() as cur:
            for rec in fri_file:
                if rec['FRIB.FINE_SEARCH'] == 'ACC':
                    accel = rec['ph_acc']
                else:
                    accel = rec['accel']

                cur.execute(query_insert, (rec['obs'],
                                           rec['start_time'],
                                           rec['stop_time'],
                                           exper,
                                           band,
                                           rec['source'],
                                           polar,
                                           rec['sta1'],
                                           rec['sta2'],
                                           rec['delay'],
                                           rec['rate'],
                                           accel,
                                           rec['SNR'],
                                           rec['ampl_lsq'],
                                           rec['duration'],
                                           rec['U'],
                                           rec['V'],
                                           rec['uv_rad_ed'],
                                           rec['ref_freq'],
                                           rec['time_code'],
                                           run_id))

            # Update status of the observations
            query_update = 'UPDATE pima_obs SET status = %s WHERE \
run_id = %s AND polar = %s AND snr <= %s;'
            cur.execute(query_update, ('n', run_id, polar, 5.3))

            query_update = 'UPDATE pima_obs SET status = %s WHERE \
run_id = %s AND polar = %s AND snr > %s;'
            if exper_info.sp_chann_num <= 128:
                cur.execute(query_update, ('y', run_id, polar, 5.7))
            else:
                cur.execute(query_update, ('y', run_id, polar, 7.0))

        self.connw.commit()

    def add_exper_info(self, exper, band, uv_fits, scan_part):
        """
        Add experiment record to the DB.

        Returns
        -------
        rec_id : int
            Id of inserted record.

        Notes
        -----
        This function deletes previous record.

        """
        with self.connw.cursor() as cursor:
            # Delete old before add new
            query = 'DELETE FROM pima_runs WHERE \
exper_name = %s AND band = %s AND fits_idi = %s AND scan_part = %s;'
            cursor.execute(query, (exper, band, uv_fits, scan_part))

            query = 'INSERT INTO pima_runs (exper_name, band, \
proc_date, fits_idi, scan_part) VALUES (%s, %s, %s, %s, %s) RETURNING id;'
            cursor.execute(query, (exper, band, datetime.now(), uv_fits,
                                   scan_part))
            rec_id = cursor.fetchone()[0]

        self.connw.commit()

        return rec_id

    def update_exper_info(self, exper_info, rec_id):
        """
        Add extended experiment information to the DB.

        """
        query = 'UPDATE pima_runs SET sp_chann_num = %s,\
 time_epochs_num = %s, scans_num = %s, obs_num = %s, uv_points_num = %s, \
uv_points_used_num = %s, deselected_points_num = %s, no_auto_points_num = %s, \
accum_length = %s, utc_minus_tai = %s, nominal_start = %s, nominal_end = %s \
WHERE id = %s'
        with self.connw.cursor() as cursor:
            cursor.execute(query, (exper_info.sp_chann_num,
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
                                   rec_id))

        self.connw.commit()

    def set_error_msg(self, run_id, msg):
        """
        Put error comment into DB.

        """
        query = 'UPDATE pima_runs SET last_error = %s WHERE id = %s'

        with self.connw.cursor() as cursor:
            cursor.execute(query, (msg, run_id))

        self.connw.commit()
