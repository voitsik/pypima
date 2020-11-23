#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 17:14:31 2018

@author: Petr Voytsik
"""

from collections import namedtuple

import pytest

from .. import Fri
from .data import SAMPLE_FRI

FriResult = namedtuple(
    "FriResult",
    [
        "length",
        "polar_0",
        "max_snr_obs",
        "max_snr_ra_obs",
        "status_0",
    ],
)

fri_results = {
    ("raes03eo", "l"): FriResult(
        1,
        "RR",
        1,
        1,
        "y",
    ),
    ("raks16nq", "c"): FriResult(
        14,
        "LL",
        9,
        12,
        "y",
    ),
}


@pytest.fixture(
    params=[("raes03eo", "l"), ("raks16nq", "c")],
    ids=["raes03eo_l", "raks16nq_c"],
    scope="module",
)
def fri_to_test(request):
    fri = Fri(SAMPLE_FRI[request.param])
    result = fri_results[request.param]

    return fri, result


class TestFri:
    """
    Test for Fri class.
    """

    def test_init(self, fri_to_test):
        fri, result = fri_to_test

        assert len(fri) == result.length
        assert fri[0]["polar"] == result.polar_0

    def test_max_snr(self, fri_to_test):
        fri, result = fri_to_test

        assert fri.max_snr()["obs"] == result.max_snr_obs
        assert fri.max_snr("RADIO-AS")["obs"] == result.max_snr_ra_obs

    def test_update_status(self, fri_to_test):
        fri, result = fri_to_test

        assert fri[0]["status"] == "u"
        fri.update_status(64)
        assert fri[0]["status"] == result.status_0
