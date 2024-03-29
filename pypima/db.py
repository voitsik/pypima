"""
Created on Thu Oct 30 14:23:23 2014

@author: Petr Voytsik
"""

import os.path
from datetime import datetime

import psycopg2
from psycopg2.extras import execute_values


class DB:
    """Interface to the "ra_results" database."""

    def __init__(self):
        """Connect to the database."""
        self.conn = psycopg2.connect(database="ra_results", user="guest", host="odin")
        self.conn.set_session(readonly=True, autocommit=True)
        self.connw = psycopg2.connect(database="ra_results", user="editor", host="odin")

    def get_uvfits_path(self, exper, band, gvlbi=False, small=False):
        """
        Retrieve path on FTP server and size of FITS-IDI file from the database.

        Parameters
        ----------
        exper : str
            Experiment name.
        band : str
            Frequency band.
        gvlbi : bool
            If ``True`` select ground only (GVLBI) FITS-file.
        small : bool
            If ``True`` select FITS files with 64 channenls only.

        Returns
        -------
        path, host, size : (str, str, int)
            Tuple of the file path, ftp hostname and size.
            Returns (None, None, 0) if the database reply is empty.

        Notes
        -----
        If database contains more than one record for given experiment and
        band, this function selects the most recent FITS-file.

        """
        path = None
        host = None
        size = 0

        query = """SELECT path, size, ftp_user FROM fits_files
        WHERE exper_name = %(exper)s AND band = %(band)s AND
        split_part(basename, '_', 1) = %(base)s AND
        LOWER(split_part(basename, '_', 2)) = %(exper)s #EXT#
        ORDER BY corr_date DESC, path DESC;"""

        params = {"exper": exper, "band": band}
        where_arr = []

        if gvlbi:
            params["base"] = "GVLBI"
        else:
            params["base"] = "RADIOASTRON"

            # Exclude rubidium
            where_arr.append("basename NOT LIKE %(rub1)s")
            params["rub1"] = "%_RUB.idifits"
            where_arr.append("basename NOT LIKE %(rub2)s")
            params["rub2"] = "%_RB.idifits"

        if small:
            where_arr.append("ch_num = %(ch_num)s")
            params["ch_num"] = 64

        if where_arr:
            query = query.replace("#EXT#", "AND {}".format(" AND ".join(where_arr)))
        else:
            query = query.replace("#EXT#", "")

        with self.conn.cursor() as cursor:
            cursor.execute(query, params)
            reply = cursor.fetchone()

        if reply:
            path = reply[0]
            size = reply[1]
            host = reply[2]

        return path, host, size

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
        url_base = "ftp://webinet.asc.rssi.ru/radioastron/oddata/reconstr/"

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
            File URL on FTP server. Returns ``None`` if database reply is empty.

        """
        url = None
        url_base = "ftp://webinet.asc.rssi.ru/radioastron/ampcal"

        query = "SELECT to_char(exper_nominal_start, 'YYYY_MM_DD') \
                 FROM vex_files WHERE exper_name = %s;"

        with self.conn.cursor() as cursor:
            cursor.execute(query, (exper,))
            reply = cursor.fetchone()

        if reply:
            date = reply[0]
            date1 = date[0:7]
            url = "{0}/{1}/{2}_{3}/{3}{4}.antab2".format(
                url_base, date1, date, exper, band
            )

        return url

    def fri2db(self, fri, exper_info, run_id, nobps=False):
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
        nobps : bool
            If ``True`` store data to the ``pima_obs_nobps`` table instead of
            ``pima_obs`` one.

        """
        if not fri:
            return

        exper = exper_info.exper
        band = exper_info.band

        if nobps:
            table = "pima_obs_nobps"
        else:
            table = "pima_obs"

        query = """INSERT INTO {} (obs, start_time, stop_time,
exper_name, band, source, polar, st1, st2, delay, rate, accel, snr, ampl,
solint, u, v, base_ed, ref_freq, scan_name, run_id, if_id, status, elevation,
bandpass, pfd) VALUES %s;""".format(
            table
        )

        rec_list = []
        for rec in fri.records:
            if rec["FRIB.FINE_SEARCH"] == "ACC":
                accel = rec["ph_acc"]
            else:
                accel = rec["accel"]

            # Set IF id to 0 if fringe fitting done over multiple IFs
            if rec["beg_ifrq"] != rec["end_ifrq"]:
                if_id = 0
            else:
                if_id = rec["beg_ifrq"]

            rec_list.append(
                (
                    rec["obs"],
                    rec["start_time"],
                    rec["stop_time"],
                    exper,
                    band,
                    rec["source"],
                    rec["polar"],
                    rec["sta1"],
                    rec["sta2"],
                    rec["delay"],
                    rec["rate"],
                    accel,
                    rec["SNR"],
                    rec["ampl_lsq"],
                    rec["duration"],
                    rec["U"],
                    rec["V"],
                    rec["uv_rad_ed"],
                    rec["ref_freq"],
                    rec["time_code"],
                    run_id,
                    if_id,
                    rec["status"],
                    rec["elevation"],
                    fri.aux["bandpass"],
                    rec["pfd"],
                )
            )

        with self.connw.cursor() as cursor:
            execute_values(cursor, query, rec_list)

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
            query = """DELETE FROM pima_runs
            WHERE exper_name = %s AND band = %s
            AND fits_idi LIKE %s AND scan_part = %s;"""
            uv_fits_toks = uv_fits.split("_")
            fits_name_base = uv_fits_toks[0] + "%"

            if len(uv_fits_toks) >= 5:
                correlator = uv_fits_toks[4]
                fits_name_base += correlator + "%"

            cursor.execute(query, (exper, band, fits_name_base, scan_part))

            query = """INSERT INTO pima_runs (exper_name, band, proc_date, fits_idi,
            scan_part)
            VALUES (%s, %s, %s, %s, %s) RETURNING id;"""
            cursor.execute(query, (exper, band, datetime.now(), uv_fits, scan_part))
            run_id = cursor.fetchone()[0]

        self.connw.commit()

        return run_id

    def update_exper_info(self, exper_info, run_id):
        """Add extended experiment information to the DB."""
        query = """UPDATE pima_runs SET
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
        pima_version = %s,
        correlator_name = %s
        WHERE id = %s
        """

        with self.connw.cursor() as cursor:
            cursor.execute(
                query,
                (
                    exper_info.sp_chann_num,
                    exper_info.time_epochs_num,
                    exper_info.scans_num,
                    exper_info.obs_num,
                    exper_info.uv_points_num,
                    exper_info.uv_points_used_num,
                    exper_info.deselected_points_num,
                    exper_info.no_auto_points_num,
                    exper_info.accum_length,
                    exper_info.utc_minus_tai,
                    exper_info.nominal_start,
                    exper_info.nominal_end,
                    exper_info.hostname,
                    exper_info.pima_version,
                    exper_info.correlator_name,
                    run_id,
                ),
            )

        self.connw.commit()

    def set_error_msg(self, run_id, msg):
        """Put error comment into DB."""
        query = "UPDATE pima_runs SET last_error = %s WHERE id = %s"

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
        query = """INSERT INTO clock_models
        (sta, time, clock_offset, clock_rate, group_delay, delay_rate, run_id)
        VALUES %s;"""

        data = []

        for rec in clock_model:
            row = list(rec)
            row.append(run_id)
            data.append(row)

        with self.connw.cursor() as cursor:
            execute_values(cursor, query, data)

        self.connw.commit()

    def uvfits2db(self, fits_file, b1950_name, run_id):
        """
        Parameters
        ----------
        fits_file : UVFits object
            UV-FITS file object.
        b1950_name : str
            B1950 source name.
        run_id : int
            Id of the record in ``pima_runs`` table.

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

                freq = float(fits_file.freq + fits_file.freq_table[if_ind]["if_freq"])
                uu = float(fits_file.u_raw[ind])
                vv = float(fits_file.v_raw[ind])

                row = (
                    ind + 1,
                    time,
                    if_ind + 1,
                    b1950_name,
                    fits_file.exper_name,
                    fits_file.band,
                    fits_file.stokes,
                    ant1_name,
                    ant2_name,
                    uu,
                    vv,
                    freq,
                    flux,
                    weight,
                    inttime,
                    file_name,
                    run_id,
                )
                data.append(row)

        query = """INSERT INTO ra_uvfits (ind, time, if_id, source, exper_name,
        band, polar, sta1, sta2, u, v, freq, ampl, weight, inttime, file_name, run_id)
        VALUES %s;"""

        if data:
            with self.connw.cursor() as cursor:
                # cursor.execute('DELETE FROM ra_uvfits WHERE file_name = %s;',
                #                (file_name, ))
                # cur.executemany(query, data)
                execute_values(cursor, query, data)
            self.connw.commit()

    def autospec2db(self, acta_file):
        """Store autocorrelation spectrum to the database."""
        exper, band = acta_file.header["experiment"].split("_")
        polar = acta_file.header["polar"]
        sta = acta_file.header["station"]
        start_date = acta_file.header["start_date"]
        stop_date = acta_file.header["stop_date"]
        obs = acta_file.header["obs"]
        scan_name = acta_file.header["scan_name"]

        delete_query = """DELETE FROM autospec_info
        WHERE exper_name = %s AND band = %s AND polar = %s AND sta = %s AND
        scan_name = %s;"""
        query_info = """INSERT INTO autospec_info
        (exper_name, band, polar, sta, start_date, stop_date, obs, scan_name)
        VALUES %s RETURNING id;"""
        query_data = """INSERT INTO autospec
        (if_num, chann_num, freq, ampl, info_id) VALUES %s;"""
        data = []

        with self.connw.cursor() as cursor:
            cursor.execute(delete_query, (exper, band, polar, sta, scan_name))
            cursor.execute(
                query_info,
                [(exper, band, polar, sta, start_date, stop_date, obs, scan_name)],
            )
            info_id = cursor.fetchone()[0]

            for ind in range(acta_file.header["num_of_points"]):
                row = (
                    acta_file.if_num[ind],
                    acta_file.channel[ind],
                    acta_file.freq[ind],
                    acta_file.ampl[ind],
                    info_id,
                )
                data.append(row)
            execute_values(cursor, query_data, data, page_size=2048)

        self.connw.commit()

    def close(self):
        """Close connections."""
        self.connw.close()
        self.conn.close()
