#!/usr/bin/env python
"""
    This module consists of functions dealing with the collection event-by-event
    results into databases.

"""

from os import path, listdir
import re
from DBR import SqliteDB

dbHolder = SqliteDB()


class EbeCollector:
    """
        This class contains functions that collect results from event-by-event
        calculations into databases.
    """

    def collectEccentricitiesAndRIntegrals(self, folder, event_id, db):
        """
            This function collects initial eccentricities and r-integrals into
            the specified SqliteDB object "db". More specifically,
            this functions fills table "ecc_id_lookup", "eccentricity", and
            "r_integrals".

            Eccentricity and r-integral files will be looked for in "folder" and
            when filling tables the specified "event_id" will be used.
        """
        # collection of file name patterns, ecc_id, and ecc_type_name
        typeCollections = (
            (
                re.compile("ecc-init-sd-r_power-(\d*).dat"), # filename pattern
                1, # ecc_id
                "sd", # ecc_type_name
            ),
            (
                re.compile("ecc-init-r_power-(\d*).dat"),
                2,
                "ed",
            )
        )
        # they have the following formats (column indices)
        ecc_real_col = 0 # real part of ecc
        ecc_imag_col = 1 # imag part of ecc
        r_inte_col = 3 # r-integral

        # first write the ecc_id_lookup table, makes sure there is only one such table
        if not db.createTableIfNotExists("ecc_id_lookup", ("ecc_id", "ecc_type_name"), ("integer", "text")):
            for pattern, ecc_id, ecc_type_name in typeCollections:
                db.insertIntoTable("ecc_id_lookup", (ecc_id, ecc_type_name))

        # next create the eccentricity and r_integrals table, if not existing
        db.createTableIfNotExists("eccentricity",
                                    ("event_id", "ecc_id", "r_power", "n", "ecc_real", "ecc_imag"),
                                    ("integer", "integer", "integer", "integer", "real", "real")
                                )
        db.createTableIfNotExists("r_integrals",
                                    ("event_id", "ecc_id", "r_power", "r_inte"),
                                    ("integer", "integer", "integer", "real", "real")
                                )

        # the big loop
        for aFile in listdir(folder): # get all file names
            for pattern, ecc_id, ecc_type_name in typeCollections: # loop over ecc types
                matchResult = pattern.match(aFile) # try to match file names
                if not matchResult: continue # not matched!
                filename = matchResult.group()
                r_power = matchResult.groups()[0] # indicated by the file name
                # read the eccentricity file and write database
                for n, aLine in enumerate(open(path.join(folder, filename))): # row index is "n"
                    data = aLine.split()
                    # insert into eccentricity table
                    db.insertIntoTable("eccentricity",
                                        (event_id, ecc_id, r_power, n, data[ecc_real_col], data[ecc_imag_col])
                                    )
                    # insert into r-integrals table but only once
                    if n==1:
                        db.insertIntoTable("r_integrals",
                                            (event_id, ecc_id, r_power, data[r_inte_col])
                                        )

        # close connection to commit changes
        db.closeConnection()


    def collectFLows_binUtilityFormat(self, folder, event_id, multiplicityFactor=1.0, db=dbHolder):
        """
            This function collects integrated and differential flows data
            and multiplicity and spectra data from "folder" into the
            database "db" using event id "event_id". The "multiplityFactor"
            will be multiplied to the number of particles read from file to
            form the multiplicity value.

            This function fills the following table: "pid_lookup",
            "inte_vn", "diff_vn", "multiplicities", "spectra".
        """
        # collection of file name patterns, pid, and particle name. The file format is determined from the "filename_format.dat" file
        typeCollections = (
            (
                re.compile("integrated_flow.dat"), # filename pattern
                0, # pid
                "charged", # particle name
            ),
        )

        # first write the pid_lookup table, makes sure there is only one such table
        if not db.createTableIfNotExists("pid_lookup", ("pid", "name"), ("integer", "text")):
            for pattern, pid, name in typeCollections:
                db.insertIntoTable("pid_lookup", (pid, name))

        # next create various tables
        db.createTableIfNotExists("inte_vn",
                                    ("event_id", "pid", "n", "vn_real", "vn_imag"),
                                    ("integer", "integer", "integer", "real", "real")
                                )
        db.createTableIfNotExists("diff_vn",
                                    ("event_id", "pid", "pT", "n", "vn_real", "vn_imag"),
                                    ("integer", "integer", "real", "integer", "real", "real")
                                )
        db.createTableIfNotExists("multiplicities",
                                    ("event_id", "pid", "N"),
                                    ("integer", "integer", "real")
                                )
        db.createTableIfNotExists("spectra",
                                    ("event_id", "pid", "pT", "N"),
                                    ("integer", "integer", "real", "real")
                                )

        # the big loop
        for aFile in listdir(folder): # get all file names
            for pattern, pid, name in typeCollections: # loop over ecc types
                !!!I'M HERE!!!
                matchResult = pattern.match(aFile) # try to match file names
                if not matchResult: continue # not matched!
                filename = matchResult.group()
                r_power = matchResult.groups()[0] # indicated by the file name
                # read the eccentricity file and write database
                for n, aLine in enumerate(open(path.join(folder, filename))): # row index is "n"
                    data = aLine.split()
                    # insert into eccentricity table
                    db.insertIntoTable("eccentricity",
                                        (event_id, ecc_id, r_power, n, data[ecc_real_col], data[ecc_imag_col])
                                    )
                    # insert into r-integrals table but only once
                    if n==1:
                        db.insertIntoTable("r_integrals",
                                            (event_id, ecc_id, r_power, data[r_inte_col])
                                        )

        # close connection to commit changes
        db.closeConnection()



if __name__ == '__main__':
    import doctest
    doctest.testfile("EbeCollector_readme.txt")