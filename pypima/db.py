"""
Created on Thu Oct 30 14:23:23 2014

@author: Petr Voytsik
"""

from datetime import datetime, timedelta
import os.path
import psycopg2


class DB:
    """
    Interface to the "ra_results" database.

    """
    def __init__(self):
        """
        Connect to the database.

        """
        self.conn = psycopg2.connect(database='ra_results', user='guest',
                                     host='odin')
        self.conn.autocommit = True
        self.connw = psycopg2.connect(database='ra_results', user='editor',
                                      host='odin')

    def get_uvfits_url(self, exper, band, gvlbi=False, small=False):
        """
        Get FITS-file URL and size from the database for the given experiment
        and band.

        Parameters
        ----------
        exper : str
            Experiment name.
        band : str
            Frequency band.
        gvlbi : bool
            If True select ground only (GVLBI) FITS-file.

        Returns
        -------
        url, size, ftp_user : (str, int, str)
            Tuple of the file URL, size and FTP user. Returns (None, 0, None)
            if the database reply is empty.

        Notes
        -----
        If database contains more than one record for given experiment and
        band, this function selects the most recent FITS-file.

        """
        url = None
        size = 0
        ftp_user = None
        url_base = 'ftp://archive.asc.rssi.ru'

        query = """
        SELECT path, size, ftp_user FROM fits_files WHERE
        LOWER(exper_name) = LOWER(%s) AND LOWER(band) = LOWER(%s) AND
        path LIKE %s #EXT#
        ORDER BY corr_date DESC, path DESC;
        """

        if small:
            query = query.replace('#EXT#', 'AND ch_num = 64')
        else:
            query = query.replace('#EXT#', '')

        with self.conn.cursor() as cursor:
            if gvlbi:
                cursor.execute(query, (exper, band, '%GVLBI%'))
            else:
                cursor.execute(query, (exper, band, '%RADIOASTRON%'))
            reply = cursor.fetchone()

        if reply:
            path = reply[0]
            size = reply[1]
            ftp_user = reply[2]
            url = url_base + path

        return url, size, ftp_user

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
        url_base = 'ftp://webinet.asc.rssi.ru/radioastron/oddata/reconstr/'

        query = """
        SELECT scf_files.file_name FROM scf_files, vex_files
        WHERE vex_files.exper_name = %s AND
        scf_files.start_time <= vex_files.exper_nominal_start AND
        scf_files.stop_time >= vex_files.exper_nominal_stop AND
        step = 1
        ORDER BY creation_date DESC, file_name DESC;
        """

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
        url_base = 'ftp://webinet.asc.rssi.ru/radioastron/ampcal'

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

    def fri2db(self, fri, exper_info, run_id):
        """
        Store information from the PIMA fri-file to the database.

        Parameters
        ----------
        fri : fri.Fri object
            PIMA fringe fitting results.
        exper_info : pima.ExperInfo object
            Experiment general information.
        run_id : int
            Id of the record in the ``pima_runs`` database table.

        """
        if not fri:
            return

        exper = exper_info.exper
        band = exper_info.band
        polar = fri.records[0]['polar']

        query = """INSERT INTO pima_obs (obs, start_time, stop_time,
exper_name, band, source, polar, st1, st2, delay, rate, accel, snr, ampl,
solint, u, v, base_ed, ref_freq, scan_name, run_id, if_id, status, elevation)
VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s,
%s, %s, %s, %s, %s, %s, %s);"""

        rec_list = []
        for rec in fri.records:
            if rec['FRIB.FINE_SEARCH'] == 'ACC':
                accel = rec['ph_acc']
            else:
                accel = rec['accel']

            # Set IF id to 0 if fringe fitting done over multiple IFs
            if rec['beg_ifrq'] != rec['end_ifrq']:
                if_id = 0
            else:
                if_id = rec['beg_ifrq']

            rec_list.append((rec['obs'],
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
                             run_id,
                             if_id,
                             rec['status'],
                             rec['elevation']))

        with self.connw.cursor() as cur:
            cur.executemany(query, rec_list)

        self.connw.commit()

    def add_exper_info(self, exper, band, uv_fits, scan_part):
        """
        Add experiment record to the DB.

        Returns
        -------
        run_id : int
            Id of inserted record.

        Notes
        -----
        This function deletes previous record.

        """
        with self.connw.cursor() as cursor:
            # Delete old before add new
            query = 'DELETE FROM pima_runs WHERE \
exper_name = %s AND band = %s AND fits_idi LIKE %s AND scan_part = %s;'
            fits_name_base = uv_fits.split('_')[0] + '%'
            cursor.execute(query, (exper, band, fits_name_base, scan_part))

            query = 'INSERT INTO pima_runs (exper_name, band, \
proc_date, fits_idi, scan_part) VALUES (%s, %s, %s, %s, %s) RETURNING id;'
            cursor.execute(query, (exper, band, datetime.now(), uv_fits,
                                   scan_part))
            run_id = cursor.fetchone()[0]

        self.connw.commit()

        return run_id

    def update_exper_info(self, exper_info, run_id):
        """
        Add extended experiment information to the DB.

        """
        query = '''UPDATE pima_runs SET
        sp_chann_num = %s,
        time_epochs_num = %s,
        scans_num = %s,
        obs_num = %s,
        uv_points_num = %s,
        uv_points_used_num = %s,
        deselected_points_num = %s,
        no_auto_points_num = %s,
        accum_length = %s,
        utc_minus_tai = %s,
        nominal_start = %s,
        nominal_end = %s,
        hostname = %s,
        pima_version = %s
        WHERE id = %s
        '''

        utc_tai = timedelta(seconds=exper_info['utc_minus_tai'])
        with self.connw.cursor() as cursor:
            cursor.execute(query, (exper_info['sp_chann_num'],
                                   exper_info['time_epochs_num'],
                                   exper_info['scans_num'],
                                   exper_info['obs_num'],
                                   exper_info['uv_points_num'],
                                   exper_info['uv_points_used_num'],
                                   exper_info['deselected_points_num'],
                                   exper_info['no_auto_points_num'],
                                   exper_info['accum_length'],
                                   utc_tai,
                                   exper_info['nominal_start'],
                                   exper_info['nominal_end'],
                                   exper_info['hostname'],
                                   exper_info['pima_version'],
                                   run_id))

        self.connw.commit()

    def set_error_msg(self, run_id, msg):
        """
        Put error comment into DB.

        """
        query = 'UPDATE pima_runs SET last_error = %s WHERE id = %s'

        with self.connw.cursor() as cursor:
            cursor.execute(query, (msg, run_id))

        self.connw.commit()

    def model2db(self, run_id, clock_model):
        """
        Parameters
        ----------
        run_id : int
            Id of record in `pima_runs` table.

        clock_model : list
            List of clock model components.

        """
        query = 'INSERT INTO clock_models \
(sta, time, clock_offset, clock_rate, group_delay, delay_rate, run_id) \
VALUES (%s, %s, %s, %s, %s, %s, %s);'

        with self.connw.cursor() as cursor:
            for rec in clock_model:
                parameters = list(rec)
                parameters.append(run_id)
                cursor.execute(query, parameters)

        self.connw.commit()

    def uvfits2db(self, fits_file, b1950_name, run_id):
        """
        Parameters
        ----------
        fits_file : UVFits object
            UV-FITS file object
        b1950_name : str
            B1950 source name

        """
        data = []
        file_name = os.path.basename(fits_file.file_name)

        for ind in range(fits_file.gcount):
            time = fits_file.dates[ind]  # As datetime.datetime object
            ant1_name, ant2_name = fits_file.get_ant_by_ind(ind)
            inttime = float(fits_file.inttime[ind])

            for if_ind in range(fits_file.no_if):
                flux = float(fits_file.amplitudes[ind, if_ind])
                weight = float(fits_file.weights[ind, if_ind])

                if weight <= 0:
                    continue

                freq = float(fits_file.freq +
                             fits_file.freq_table[if_ind]['if_freq'])
                uu = float(fits_file.u_raw[ind])
                vv = float(fits_file.v_raw[ind])

                row = (ind+1, time, if_ind+1, b1950_name, fits_file.exper_name,
                       fits_file.band, fits_file.stokes, ant1_name, ant2_name,
                       uu, vv, freq, flux, weight, inttime, file_name, run_id)
                data.append(row)

        query = """INSERT INTO ra_uvfits (ind, time, if_id, source, exper_name,
band, polar, sta1, sta2, u, v, freq, ampl, weight, inttime, file_name, run_id)
VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);"""

        if data:
            with self.connw.cursor() as cur:
                cur.execute('DELETE FROM ra_uvfits WHERE file_name = %s;',
                            (file_name, ))
                cur.executemany(query, data)
            self.connw.commit()
