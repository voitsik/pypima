#!/usr/bin/env python3
"""
Created on Fri Aug 24 16:09:11 2018

@author: Petr Voytsik
"""

import os.path
import shutil
import re
from typing import NamedTuple, List, Tuple

import pytest

from .. import Fri, Pima, PimaError
from .data import SAMPLE_FRI, SAMPLE_RA_CNT, SAMPLE_RA_FITS, SAMPLE_SCF


class LoadResult(NamedTuple):
    sp_chann_num: int
    ap_minmax: float
    number_of_deselected_points: int
    station_list_ivs: List[str]
    station_list_corr: List[str]
    source_names: Tuple[str]
    obs_number: int


@pytest.fixture(scope="module")
def work_dir(tmpdir_factory):
    """
    Create and prepare a temporary working directory.
    """
    wdir = tmpdir_factory.mktemp("work_dir")
    sdir = tmpdir_factory.mktemp("pima_scr")

    for cnt_file in SAMPLE_RA_CNT.values():
        shutil.copy(cnt_file, str(wdir))

    return (str(wdir), str(sdir))


def fri_rec_equivalent(rec1, rec2):
    """
    Check two records of fri-file for equivalence.
    """
    return (
        (rec1["exper_code"] == rec2["exper_code"])
        and (rec1["session_code"] == rec2["session_code"])
        and (rec1["source"] == rec2["source"])
        and (rec1["sta1"] == rec2["sta1"])
        and (rec1["sat2"] == rec2["sta2"])
    )


class TestPima:
    load_results = {
        ("raes03eo", "l"): LoadResult(
            64,
            (1.0, 1.0),
            0,
            ["RADIO-AS", "ZELENCHK"],
            ["RA", "ZC"],
            ("0529+483", "J0533+4822", "0529+483"),
            1,
        ),
        ("raks16nq", "c"): LoadResult(
            64,
            (0.5, 0.5),
            0,
            ["EFLSBERG", "IRBENE16", "RADIO-AS", "TORUN"],
            ["EF", "IB", "RA", "TR"],
            ("0716+714", "J0721+7120", "0716+714"),
            14,
        ),
    }

    def test_environ(self):
        """
        Test environment variables.
        """
        pima_dir = os.getenv("PIMA_DIR")
        assert pima_dir is not None
        pima_exec = os.path.join(pima_dir, "bin", "pima")
        assert os.path.isfile(pima_exec)

    @pytest.mark.parametrize(("exper", "band"), [("raes03eo", "l"), ("raks16nq", "c")])
    def test_pima_init(self, work_dir, exper, band):
        """
        Test Pima.__init__
        """
        wdir, _ = work_dir
        pima = Pima(exper, band, work_dir=wdir)

        assert pima.exper == exper
        assert pima.band == band
        assert os.path.basename(pima.cnt_file_name) == os.path.basename(
            SAMPLE_RA_CNT[exper, band]
        )

    @pytest.mark.parametrize(("exper", "band"), [("raes03eo", "l"), ("raks16nq", "c")])
    def test_pima_update_cnt(self, work_dir, exper, band):
        """
        Test pima.update_cnt
        """
        wdir, sdir = work_dir
        pima = Pima(exper, band, work_dir=wdir)

        fits_file = SAMPLE_RA_FITS[exper, band]
        scf_file = SAMPLE_SCF[exper]
        pima_dir = os.getenv("PIMA_DIR")
        source_names = os.path.join(pima_dir, "share", "vtd_data", "source.names")
        sta_names = os.path.join(pima_dir, "share", "pima", "ra_station.names")
        fft_wis_file = os.path.join(
            pima_dir, "share", "pima", "pima_small_measure_1thr.wis"
        )
        vtd_config = os.path.join(pima_dir, "share", "vtd_data", "vtd_pima.cnf")

        pima.update_cnt(
            {
                "UV_FITS:": fits_file,
                "STAGING_DIR:": "NO",
                "EXPER_DIR:": str(sdir),
                "EPHEMERIDES_FILE:": scf_file,
                "SOU_NAMES:": source_names,
                "STA_NAMES:": sta_names,
                "FFT_CONFIG_FILE:": fft_wis_file,
                "VTD_CONFIG_FILE:": vtd_config,
                "MKDB.VCAT_CONFIG:": "NO",
            }
        )

        assert pima.cnt_params["UV_FITS:"] == [fits_file]
        assert pima.cnt_params["STAGING_DIR:"] == "NO"
        assert pima.cnt_params["EXPER_DIR:"] == str(sdir)
        assert pima.cnt_params["EPHEMERIDES_FILE:"] == scf_file
        assert pima.cnt_params["SOU_NAMES:"] == source_names
        assert pima.cnt_params["STA_NAMES:"] == sta_names
        assert pima.cnt_params["FFT_CONFIG_FILE:"] == fft_wis_file
        assert pima.cnt_params["VTD_CONFIG_FILE:"] == vtd_config

    @pytest.mark.parametrize(("exper", "band"), [("raes03eo", "l"), ("raks16nq", "c")])
    def test_pima_load(self, work_dir, exper, band):
        """Test pima.load."""
        wdir, _ = work_dir
        pima = Pima(exper, band, work_dir=wdir)

        pima.load()

    @pytest.mark.parametrize(("exper", "band"), [("raes03eo", "l"), ("raks16nq", "c")])
    def test_pima_exper_info(self, work_dir, exper, band):
        """Test pima.exper_info and Pima properties."""
        wdir, _ = work_dir
        pima = Pima(exper, band, work_dir=wdir)

        result = self.load_results[exper, band]

        assert pima.exper_info.exper == exper
        assert pima.exper_info.band == band
        assert pima.exper_info.sp_chann_num == result.sp_chann_num
        assert (
            pima.exper_info.deselected_points_num
            == result.number_of_deselected_points
        )

        # Check PIMA version format
        assert re.fullmatch(r"\d+\.\d+\w*", pima.exper_info.pima_version, re.ASCII)

        assert pima.ap_minmax == result.ap_minmax
        assert pima.number_of_deselected_points == result.number_of_deselected_points
        assert pima.station_list(ivs_name=False) == result.station_list_corr
        assert pima.station_list(ivs_name=True) == result.station_list_ivs
        assert pima.chan_number == result.sp_chann_num
        assert pima.obs_number == result.obs_number
        assert pima.source_list[0] == result.source_names

    @pytest.mark.parametrize(("exper", "band"), [("raes03eo", "l"), ("raks16nq", "c")])
    def test_pima_set_frq_grp(self, work_dir, exper, band):
        """Test pima.set_frq_grp."""
        wdir, _ = work_dir
        pima = Pima(exper, band, work_dir=wdir)

        pima.set_frq_grp(1)

        assert int(pima.cnt_params["FRQ_GRP:"]) == 1

        with pytest.raises(PimaError):
            pima.set_frq_grp(2)

    @pytest.mark.parametrize(("exper", "band"), [("raes03eo", "l"), ("raks16nq", "c")])
    def test_pima_coarse(self, work_dir, exper, band):
        """Test pima.coarse."""
        wdir, _ = work_dir
        pima = Pima(exper, band, work_dir=wdir)

        fri_file = pima.coarse()

        assert os.path.isfile(fri_file)

        fri = Fri(fri_file)
        ref_fri = Fri(SAMPLE_FRI[exper, band, "nobps"])

        rec1 = fri[0]
        rec2 = ref_fri[0]

        assert rec1["obs"] == rec2["obs"]
        assert rec1["exper_code"] == rec2["exper_code"]
        assert rec1["session_code"] == rec2["session_code"]
        assert rec1["polar"] == rec2["polar"]
        assert rec1["source"] == rec2["source"]
        assert rec1["sta1"] == rec2["sta1"]
        assert rec1["sta2"] == rec2["sta2"]
        assert rec1["SNR"] == pytest.approx(rec2["SNR"], rel=1e-2)
        assert rec1["ampl_lsq"] == pytest.approx(rec2["ampl_lsq"], rel=1e-2)
        assert rec1["delay"] == pytest.approx(rec2["delay"], abs=5e-9)
        assert rec1["rate"] == pytest.approx(rec2["rate"], abs=5e-14)
        assert rec1["accel"] == pytest.approx(rec2["accel"], abs=1e-17)
        assert rec1["U"] == pytest.approx(rec2["U"])
        assert rec1["V"] == pytest.approx(rec2["V"])

    @pytest.mark.parametrize(("exper", "band"), [("raes03eo", "l"), ("raks16nq", "c")])
    def test_pima_bpas(self, work_dir, exper, band):
        """Test pima.bpas."""
        wdir, _ = work_dir
        pima = Pima(exper, band, work_dir=wdir)

        mseg = pima.chan_number // 2

        bpas_params = {
            "BPS.MODE:": "ACCUM",
            "BPS.NOBS_ACCUM:": "6",
            "BPS.MSEG_ACCUM:": mseg,
            "BPS.NOBS_FINE:": "12",
            "BPS.MINOBS_FINE:": "8",
            "BPS.MSEG_FINE:": mseg,
            "BPS.SNR_MIN_ACCUM:": "5.5",
            "BPS.SNR_MIN_FINE:": "5.5",
            "BPS.AMPL_REJECT:": "0.4",
            "BPS.PHAS_REJECT:": "0.2",
            "BPS.INTRP_METHOD:": "LINEAR",
            "BPS.DEG_AMP:": "0",
            "BPS.DEG_PHS:": "1",
            "BPS.AMP_MIN:": "0.01",
            "BPS.NORML:": "IF",
            "BPS.SEFD_USE:": "NO",
        }

        pima.update_cnt(bpas_params)

        pima.bpas()

        bps_file = pima.cnt_params["BANDPASS_FILE:"]
        assert os.path.isfile(bps_file)

    @pytest.mark.parametrize(("exper", "band"), [("raes03eo", "l"), ("raks16nq", "c")])
    def test_pima_fine(self, work_dir, exper, band):
        """Test pima.fine."""
        wdir, _ = work_dir
        pima = Pima(exper, band, work_dir=wdir)

        polar = pima.cnt_params["POLAR:"]
        fri_file = f"{exper}_{band}_{polar}.fri"
        fri_file = os.path.join(wdir, fri_file)
        frr_file = f"{exper}_{band}_{polar}.frr"
        frr_file = os.path.join(wdir, frr_file)
        pima.update_cnt({"FRINGE_FILE:": fri_file, "FRIRES_FILE:": frr_file})

        fri_file = pima.fine()

        assert os.path.isfile(fri_file)

        ref_fri = Fri(SAMPLE_FRI[exper, band])
        fri = Fri(fri_file)

        rec1 = fri[0]
        rec2 = ref_fri[0]

        assert rec1["obs"] == rec2["obs"]
        assert rec1["exper_code"] == rec2["exper_code"]
        assert rec1["session_code"] == rec2["session_code"]
        assert rec1["polar"] == rec2["polar"]
        assert rec1["source"] == rec2["source"]
        assert rec1["sta1"] == rec2["sta1"]
        assert rec1["sta2"] == rec2["sta2"]
        assert rec1["SNR"] == pytest.approx(rec2["SNR"], rel=1e-2)
        assert rec1["ampl_lsq"] == pytest.approx(rec2["ampl_lsq"], rel=1e-2)
        assert rec1["delay"] == pytest.approx(rec2["delay"], abs=5e-9)
        assert rec1["rate"] == pytest.approx(rec2["rate"], abs=5e-14)
        assert rec1["accel"] == pytest.approx(rec2["accel"], abs=1e-17)
        assert rec1["U"] == pytest.approx(rec2["U"])
        assert rec1["V"] == pytest.approx(rec2["V"])
