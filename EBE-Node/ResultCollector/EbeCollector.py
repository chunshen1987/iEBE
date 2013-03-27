#!/usr/bin/env python
"""
    This module consists of functions dealing with the collection event-by-event
    results into databases.

"""

from os import path
from sys import argv
from DBR import SqliteDB

import numpy as np


class EbeCollector:
    """
        This class contains functions that collect results from event-by-event
        calculations into databases.
    """

    def collectEccentricitiesAndRIntegrals(self, sqliteDBDatabase, folder, eventId):
        """
            This function collects initial eccentricities and r-integrals into
            the specified SqliteDB object "sqliteDBDatabase". More specifically,
            this functions fills table "ecc_id_lookup", "eccentricity", and
            "r_integrals".

            Eccentricity and r-integral files will be looked for in "folder" and
            when filling tables the specified event_id="eventId" will be used.
        """
        # sd and ed eccentricity file names have the following patterns
        sd_ecc_file_re = "ecc-init-sd-r-power-(\d*).dat"
        ed_ecc_file_re = "ecc-init-r-power-(\d*).dat"



