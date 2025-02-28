"""
Interact with the RadioAstron database.

Created on Thu Oct 30 14:23:23 2014

@author: Petr Voytsik
"""

import logging
import os.path
from datetime import datetime
from typing import Optional

from sqlalchemy import MetaData, Table, create_engine, func, select

from .fri import Fri  # for type hints
from .pima import ActaFile, ExperInfo  # for type hints
from .uvfits import UVFits  # for type hints

logger = logging.getLogger(__name__)


class DataBase:
    """Interface to the "ra_results" database."""

    def __init__(self):
        """Connect to the database."""
        self.engine = create_engine(
            "postgresql://editor@localhost/ra_results",
            insertmanyvalues_page_size=2048,
        )
        self.metadata = MetaData()

        self.fits_files = Table("fits_files", self.metadata, autoload_with=self.engine)
        self.pima_runs = Table("pima_runs", self.metadata, autoload_with=self.engine)
        self.pima_obs = Table("pima_obs", self.metadata, autoload_with=self.engine)
        self.pima_obs_nobps = Table(
            "pima_obs_nobps", self.metadata, autoload_with=self.engine
        )
        self.clock_models = Table(
            "clock_models", self.metadata, autoload_with=self.engine
        )
        self.ra_uvfits = Table("ra_uvfits", self.metadata, autoload_with=self.engine)
        self.autospec_info = Table(
            "autospec_info", self.metadata, autoload_with=self.engine
        )
        self.autospec = Table("autospec", self.metadata, autoload_with=self.engine)
        self.scf_files = Table("scf_files", self.metadata, autoload_with=self.engine)
        self.vex_files = Table("vex_files", self.metadata, autoload_with=self.engine)

    def get_uvfits_path(
        self, exper: str, band: str, gvlbi: bool = False, small: bool = False
    ) -> tuple[str, str, int]:
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
            Returns ("", "", 0) if the database reply is empty.

        Notes
        -----
        If database contains more than one record for given experiment and
        band, this function selects the most recent FITS-file.

        """
        path = ""
        host = ""
        size = 0

        query = select(
            self.fits_files.c.path, self.fits_files.c.size, self.fits_files.c.ftp_user
        ).where(
            self.fits_files.c.exper_name == exper,
            self.fits_files.c.band == band,
            func.split_part(self.fits_files.c.basename, "_", 1)
            == ("GVLBI" if gvlbi else "RADIOASTRON"),
            func.lower(func.split_part(self.fits_files.c.basename, "_", 2))
            == exper.lower(),
        )

        if not gvlbi:
            # Exclude rubidium
            query = query.where(
                ~self.fits_files.c.basename.like("%_RUB.idifits"),
                ~self.fits_files.c.basename.like("%_RB.idifits"),
            )

        if small:
            query = query.where(self.fits_files.c.ch_num == 64)

        query = query.order_by(
            self.fits_files.c.corr_date.desc(), self.fits_files.c.path.desc()
        )

        logger.debug("%s: query:\n%s", "get_uvfits_path", query)

        with self.engine.connect() as conn:
            row = conn.execute(query).first()

        if row:
            path, size, host = row

        return path, host, size

    def get_orbit_url(self, exper: str) -> Optional[str]:
        """
        Return reconstructed orbit file URL for the given experiment `exper`.

        Parameters
        ----------
        exper : str
            Experiment name.

        Returns
        -------
        url : str
            File URL on FTP server. Returns ``None`` the database reply is empty.

        """
        url = None
        url_base = "ftp://webinet.asc.rssi.ru/radioastron/oddata/reconstr"

        query = (
            select(self.scf_files.c.file_name)
            .join(self.vex_files, self.vex_files.c.exper_name == exper)
            .where(
                self.scf_files.c.start_time <= self.vex_files.c.exper_nominal_start,
                self.scf_files.c.stop_time >= self.vex_files.c.exper_nominal_stop,
                self.scf_files.c.step == 1,
            )
            .order_by(
                self.scf_files.c.creation_date.desc(), self.scf_files.c.file_name.desc()
            )
        )

        logger.debug("%s: query:\n%s", "get_orbit_url", query)

        with self.engine.connect() as conn:
            reply = conn.execute(query).scalar()

        if reply:
            url = f"{url_base}/{reply}"

        return url

    def get_antab_url(self, exper: str, band: str) -> Optional[str]:
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

        query = select(
            func.to_char(self.vex_files.c.exper_nominal_start, "YYYY_MM_DD").label(
                "date"
            )
        ).where(self.vex_files.c.exper_name == exper)

        logger.debug("%s: query:\n%s", "get_antab_url", query)

        with self.engine.connect() as conn:
            reply = conn.execute(query).scalar()

        if reply:
            date_year_month_day = reply
            date_year_month = date_year_month_day[0:7]
            url = "{0}/{1}/{2}_{3}/{3}{4}.antab2".format(
                url_base, date_year_month, date_year_month_day, exper, band
            )

        return url

    def fri2db(
        self, fri: Fri, exper_info: ExperInfo, run_id: int, nobps: bool = False
    ) -> None:
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

        table = self.pima_obs_nobps if nobps else self.pima_obs

        data = []
        for rec in fri.records:
            if rec["FRIB.FINE_SEARCH"] == "ACC":
                accel = rec["ph_acc"]
            else:
                accel = rec["accel"]

            # Set IF id to 0 if fringe fitting done over multiple IFs
            if_id = 0 if rec["beg_ifrq"] != rec["end_ifrq"] else rec["beg_ifrq"]

            data.append(
                {
                    "obs": rec["obs"],
                    "scan_name": rec["time_code"],
                    "start_time": rec["start_time"],
                    "stop_time": rec["stop_time"],
                    "exper_name": exper_info.exper,
                    "band": exper_info.band,
                    "source": rec["source"],
                    "polar": rec["polar"],
                    "st1": rec["sta1"],
                    "st2": rec["sta2"],
                    "delay": rec["delay"],
                    "rate": rec["rate"],
                    "accel": accel,
                    "snr": rec["SNR"],
                    "ampl": rec["ampl_lsq"],
                    "solint": rec["duration"],
                    "u": rec["U"],
                    "v": rec["V"],
                    "base_ed": rec["uv_rad_ed"],
                    "ref_freq": rec["ref_freq"],
                    "run_id": run_id,
                    "if_id": if_id,
                    "status": rec["status"],
                    "elevation": rec["elevation"],
                    "bandpass": fri.aux["bandpass"],
                    "pfd": rec["pfd"],
                }
            )

        if data:
            with self.engine.begin() as conn:
                conn.execute(table.insert(), data)

    def add_exper_info(
        self, exper: str, band: str, uv_fits: str, scan_part: int
    ) -> int:
        """
        Add experiment record to the database.

        Parameters
        ----------
        exper : str
            Experiment name.
        band : str
            Frequency band.
        uv_fits : str
            Name of the UV-FITS file.
        scan_part : int
            Number of PIMA run (1 = full scan, 2 = half of scan, etc.).

        Returns
        -------
        run_id : int
            Id of the inserted record.

        Notes
        -----
        This function deletes previous record.

        """
        uv_fits_toks = uv_fits.split("_")
        fits_name_base = uv_fits_toks[0] + "%"  # RADIOASTRON% or GVLBI%

        if len(uv_fits_toks) >= 5:
            correlator = uv_fits_toks[4]  # ASC or DIFX
            fits_name_base += correlator + "%"

        # Delete query
        delete_stmt = self.pima_runs.delete().where(
            self.pima_runs.c.exper_name == exper,
            self.pima_runs.c.band == band,
            self.pima_runs.c.fits_idi.like(fits_name_base),
            self.pima_runs.c.scan_part == scan_part,
        )
        logger.debug("add_exper_info: delete_stmt:\n%s", delete_stmt)

        # Insert query
        insert_stmt = (
            self.pima_runs.insert()
            .values(
                exper_name=exper,
                band=band,
                proc_date=datetime.now(),
                fits_idi=uv_fits,
                scan_part=scan_part,
            )
            .returning(self.pima_runs.c.id)
        )
        logger.debug("add_exper_info: insert_stmt:\n%s", insert_stmt)

        with self.engine.begin() as conn:
            conn.execute(delete_stmt)
            result = conn.execute(insert_stmt)
            run_id = result.scalar()

        return run_id

    def update_exper_info(self, exper_info: ExperInfo, run_id: int) -> None:
        """Add extended experiment information to the database."""
        update_stmt = (
            self.pima_runs.update()
            .where(self.pima_runs.c.id == run_id)
            .values(
                sp_chann_num=exper_info.sp_chann_num,
                time_epochs_num=exper_info.time_epochs_num,
                scans_num=exper_info.scans_num,
                obs_num=exper_info.obs_num,
                uv_points_num=exper_info.uv_points_num,
                uv_points_used_num=exper_info.uv_points_used_num,
                deselected_points_num=exper_info.deselected_points_num,
                no_auto_points_num=exper_info.no_auto_points_num,
                accum_length=exper_info.accum_length,
                utc_minus_tai=exper_info.utc_minus_tai,
                nominal_start=exper_info.nominal_start,
                nominal_end=exper_info.nominal_end,
                hostname=exper_info.hostname,
                pima_version=exper_info.pima_version,
                correlator_name=exper_info.correlator_name,
            )
        )
        logger.debug("update_exper_info: update_stmt:\n%s", update_stmt)

        with self.engine.begin() as conn:
            conn.execute(update_stmt)

    def set_error_msg(self, run_id: int, msg: str) -> None:
        """Write error `msg` to the database."""
        update_stmt = (
            self.pima_runs.update()
            .where(self.pima_runs.c.id == run_id)
            .values(last_error=msg)
        )
        logger.debug("set_error_msg: update_stmt:\n%s", update_stmt)

        with self.engine.begin() as conn:
            conn.execute(update_stmt)

    def model2db(self, run_id: int, clock_model: list) -> None:
        """
        Store clock model to the database.

        Parameters
        ----------
        run_id : int
            Id of record in `pima_runs` table.

        clock_model : list
            List of clock model components.

        """
        data = []
        for rec in clock_model:
            data.append(
                {
                    "sta": rec.sta,
                    "time": rec.time,
                    "clock_offset": rec.clock_offset,
                    "clock_rate": rec.clock_rate,
                    "group_delay": rec.group_delay,
                    "delay_rate": rec.delay_rate,
                    "run_id": run_id,
                }
            )

        if data:
            with self.engine.begin() as conn:
                conn.execute(self.clock_models.insert(), data)

    def uvfits2db(self, fits_file: UVFits, b1950_name: str, run_id: int) -> None:
        """
        Store UV-data from `fits_file` to the database.

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
                ampl = float(fits_file.amplitudes[ind, if_ind])
                weight = float(fits_file.weights[ind, if_ind])

                if weight <= 0:
                    continue

                freq = float(fits_file.freq + fits_file.freq_table[if_ind]["if_freq"])
                uu = float(fits_file.u_raw[ind])
                vv = float(fits_file.v_raw[ind])

                data.append(
                    {
                        "ind": ind + 1,
                        "time": time,
                        "if_id": if_ind + 1,
                        "source": b1950_name,
                        "exper_name": fits_file.exper_name,
                        "band": fits_file.band,
                        "polar": fits_file.stokes,
                        "sta1": ant1_name,
                        "sta2": ant2_name,
                        "u": uu,
                        "v": vv,
                        "freq": freq,
                        "ampl": ampl,
                        "weight": weight,
                        "inttime": inttime,
                        "file_name": file_name,
                        "run_id": run_id,
                    }
                )

        if data:
            with self.engine.begin() as conn:
                conn.execute(self.ra_uvfits.insert(), data)

    def autospec2db(self, acta_file: ActaFile) -> None:
        """Store autocorrelation spectrum to the database."""
        exper, band = acta_file.header["experiment"].split("_")
        polar = acta_file.header["polar"]
        sta = acta_file.header["station"]
        start_date = acta_file.header["start_date"]
        stop_date = acta_file.header["stop_date"]
        obs = acta_file.header["obs"]
        scan_name = acta_file.header["scan_name"]

        # Delete existing records
        delete_stmt = self.autospec_info.delete().where(
            self.autospec_info.c.exper_name == exper,
            self.autospec_info.c.band == band,
            self.autospec_info.c.polar == polar,
            self.autospec_info.c.sta == sta,
            self.autospec_info.c.scan_name == scan_name,
        )

        # Insert new info record
        insert_info_stmt = (
            self.autospec_info.insert()
            .values(
                exper_name=exper,
                band=band,
                polar=polar,
                sta=sta,
                start_date=start_date,
                stop_date=stop_date,
                obs=obs,
                scan_name=scan_name,
            )
            .returning(self.autospec_info.c.id)
        )

        with self.engine.begin() as conn:
            conn.execute(delete_stmt)
            result = conn.execute(insert_info_stmt)
            info_id = result.scalar()

            # Insert autospec data
            data = [
                {
                    "if_num": acta_file.if_num[ind],
                    "chann_num": acta_file.channel[ind],
                    "freq": acta_file.freq[ind],
                    "ampl": acta_file.ampl[ind],
                    "info_id": info_id,
                }
                for ind in range(acta_file.header["num_of_points"])
            ]

            if data:
                conn.execute(self.autospec.insert(), data)

    def close(self):
        """Close connections."""
        self.engine.dispose()
