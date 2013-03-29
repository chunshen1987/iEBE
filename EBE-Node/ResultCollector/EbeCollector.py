#!/usr/bin/env python2
"""
    This module consists of functions dealing with the collection event-by-event
    results into databases.

"""

from os import path, listdir
import re
from DBR import SqliteDB
from assignmentFormat import assignmentExprStream2IndexDict

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
        if db.createTableIfNotExists("ecc_id_lookup", ("ecc_id", "ecc_type_name"), ("integer", "text")):
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
                                        (event_id, ecc_id, r_power, n, float(data[ecc_real_col]), float(data[ecc_imag_col]))
                                    )
                    # insert into r-integrals table but only once
                    if n==1:
                        db.insertIntoTable("r_integrals",
                                            (event_id, ecc_id, r_power, float(data[r_inte_col]))
                                        )

        # close connection to commit changes
        db.closeConnection()


    def collectFLowsAndMultiplicities_binUtilityFormat(self, folder, event_id, db, multiplicityFactor=1.0):
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
        pidDict = {
            "Charged"       : 0, # particle name, pid
            "Pion"          : 212,
            "Kaon"          : 321,
            "Proton"        : 2212,
        }
        filePattern = re.compile("([a-zA-z]*)_flow_([a-zA-Z+]*).dat") # filename pattern, the 2nd matched string needs to be among the pidTable above in order to be considered "matched"; the 1st matched string will either be "integrated" or "differential"
        tableChooser = { # will be used to decide which table to write to
                "integrated"    :   ("inte_vn", "multiplicities"),
                "differential"  :   ("diff_vn", "spectra"),
            }

        # next read in file format, which is assumed to be stored in the file "integrated_flow_format.dat" and "differential_flow_format.dat" (same)
        fmt = assignmentExprStream2IndexDict(open(path.join(folder, "integrated_flow_format.dat"))) # column index will automatically be 0-based
        N_col = fmt["count"] # number of particles for the given condition (diff or inte)
        pT_col = fmt["pT_mean_real"]
        vn_real_cols = {} # will have items (n, column index)
        vn_imag_cols = {}
        # probe for the largest n value
        largest_n = 1
        allFields = fmt.keys()
        while ("v_%d_mean_real" % largest_n) in allFields:
            vn_real_cols[largest_n] = fmt["v_%d_mean_real" % largest_n]
            vn_imag_cols[largest_n] = fmt["v_%d_mean_imag" % largest_n]
            largest_n += 1

        # first write the pid_lookup table, makes sure there is only one such table
        if db.createTableIfNotExists("pid_lookup", ("name", "pid"), ("text","integer")):
            db.insertIntoTable("pid_lookup", pidDict.items())

        # next create various tables
        db.createTableIfNotExists("inte_vn",
                                    ("event_id", "pid", "pT", "n", "vn_real", "vn_imag"),
                                    ("integer", "integer", "real", "integer", "real", "real")
                                )
        db.createTableIfNotExists("diff_vn",
                                    ("event_id", "pid", "pT", "n", "vn_real", "vn_imag"),
                                    ("integer", "integer", "real", "integer", "real", "real")
                                )
        db.createTableIfNotExists("multiplicities",
                                    ("event_id", "pid", "pT", "N"),
                                    ("integer", "integer", "real", "real")
                                )
        db.createTableIfNotExists("spectra",
                                    ("event_id", "pid", "pT", "N"),
                                    ("integer", "integer", "real", "real")
                                )

        # the big loop
        for aFile in listdir(folder): # get all file names
            matchResult = filePattern.match(aFile) # try to match file names
            if not matchResult: continue # not matched!
            flow_type, particle_name = matchResult.groups() # indicated by the file name
            if particle_name not in pidDict.keys(): continue # dont know about this particle
            pid = pidDict[particle_name] # get pid
            filename = matchResult.group() # get the file to be opened
            flow_table, multiplicity_table = tableChooser[flow_type] # choose tables to write to
            # read the flow file and write results
            for aLine in open(path.join(folder, filename)):
                data = aLine.split()
                # write flow table
                for n in range(1, largest_n):
                    db.insertIntoTable(flow_table,
                                        (event_id, pid, float(data[pT_col]), n, float(data[vn_real_cols[n]]), float(data[vn_imag_cols[n]]))
                                    )
                # write multiplicity table
                db.insertIntoTable(multiplicity_table,
                                        (event_id, pid, float(data[pT_col]), float(data[N_col])*multiplicityFactor)
                                    )

        # close connection to commit changes
        db.closeConnection()



if __name__ == '__main__':
    import doctest
    doctest.testfile("EbeCollector_readme.txt")
