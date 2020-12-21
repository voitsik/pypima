#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from distutils.core import setup

import pypima

setup(
    name="pypima",
    version=pypima.__version__,
    description="Python interface to PIMA",
    author="Petr Voytsik",
    author_email="voitsik@asc.rssi.ru",
    url="https://github.com/voitsik/pypima",
    packages=["pypima"],
    scripts=[
        "bin/pima_coher.py",
        "bin/pima_plot_fringe3d.py",
        "bin/pypima_batch.py",
        "bin/pypima_single.py",
    ],
)
