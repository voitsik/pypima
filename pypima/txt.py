# -*- coding: utf-8 -*-
"""
Created on Fri Apr 18 01:05:23 MSK 2014

@author: Petr Voytsik
"""

import array


class TextTable1D(object):
    """
    Read 1D PIMA graphics in TXT format.

    """

    FORMAT_STRING = "# 1D text table.  Format version of 2012.12.30"

    def __init__(self, file_name):

        self.plot_title = ""
        self.subtitle = ""
        self.num_points = 0
        self.axis1_title = ""
        self.axis1_units = ""
        self.axis1_min = 0.0
        self.axis1_max = 0.0
        self.axis2_title = ""
        self.axis2_units = ""
        self.axis2_min = 0.0
        self.axis2_max = 0.0
        self.axis1_data = array.array("d")
        self.axis2_data = array.array("d")

        with open(file_name, "r") as fil:
            magic = fil.readline().strip()
            if magic != self.FORMAT_STRING:
                raise Exception('Bad format string in file "{}"'.format(file_name))
            for line in fil:
                if line.startswith("#"):
                    continue

                key, val = line.split(":", 1)
                key = key.strip()
                val = val.strip()
                if key == "PLOT_TITLE":
                    self.plot_title = val
                elif key == "SUBTITLE":
                    self.subtitle = val
                elif key == "NUM_POINTS":
                    self.num_points = int(val)
                elif key == "AXIS1_TITLE":
                    self.axis1_title = val
                elif key == "AXIS1_UNITS":
                    self.axis1_units = val
                elif key == "AXIS1_MIN":
                    self.axis1_min = float(val.replace("D", "e"))
                elif key == "AXIS1_MAX":
                    self.axis1_max = float(val.replace("D", "e"))
                elif key == "AXIS2_NAME":
                    self.axis2_title = val
                elif key == "AXIS2_UNITS":
                    self.axis2_units = val
                elif key == "AXIS2_MIN":
                    self.axis2_min = float(val.replace("D", "e"))
                elif key == "AXIS2_MAX":
                    self.axis2_max = float(val.replace("D", "e"))
                elif key == "POINT":
                    data = val.split()
                    self.axis1_data.append(float(data[1].replace("D", "e")))
                    self.axis2_data.append(float(data[2].replace("D", "e")))
