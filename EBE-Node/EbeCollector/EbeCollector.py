#!/usr/bin/env python2
"""
    This module consists of functions dealing with the collection event-by-event
    results into databases.

"""

from os import path, listdir
import re
import numpy as np
from DBR import SqliteDB
from assignmentFormat import assignmentExprStream2IndexDict

class EbeCollector:
    """
        This class contains functions that collect results from event-by-event
        calculations into databases.
    """
    def __init__(self):
        """
            Define class-wise constants and tables.
        """
        self.pidDict = { # particle_name, pid
            "total"             :   0,
            "charged"           :   1,
            "pion"              :   6, # sum(7, 8, -7)
            "pion_p"            :   7,
            "pion_0"            :   8,
            "pion_m"            :   -7,
            "kaon"              :   11, # sum(12, 13)
            "kaon_p"            :   12,
            "kaon_0"            :   13,
            "anti_kaon"         :   -11, # sum(-12, -13)
            "kaon_m"            :   -12,
            "anti_kaon_0"       :   -13,
            "nucleon"           :   16, # sum(17, 18)
            "proton"            :   17,
            "neutron"           :   18,
            "anti_nucleon"      :   -16, # sum(-17, -18)
            "anti_proton"       :   -17,
            "anit_neutron"      :   -18,

            "charged_hydro"           :   301,
            "pion_hydro"              :   306, # sum(7, 8, -7)
            "pion_p_hydro"            :   307,
            "pion_0_hydro"            :   308,
            "pion_m_hydro"            :   -307,
            "kaon_hydro"              :   311, # sum(12, 13)
            "kaon_p_hydro"            :   312,
            "kaon_0_hydro"            :   313,
            "anti_kaon_hydro"         :   -311, # sum(-12, -13)
            "kaon_m_hydro"            :   -312,
            "anti_kaon_0_hydro"       :   -313,
            "nucleon_hydro"           :   316, # sum(17, 18)
            "proton_hydro"            :   317,
            "neutron_hydro"           :   318,
            "anti_nucleon_hydro"      :   -316, # sum(-17, -18)
            "anti_proton_hydro"       :   -317,
            "anit_neutron_hydro"      :   -318,

            "pion_thermal"              :   606, # sum(7, 8, -7)
            "pion_p_thermal"            :   607,
            "pion_0_thermal"            :   608,
            "pion_m_thermal"            :   -607,
            "kaon_thermal"              :   611, # sum(12, 13)
            "kaon_p_thermal"            :   612,
            "kaon_0_thermal"            :   613,
            "anti_kaon_thermal"         :   -611, # sum(-12, -13)
            "kaon_m_thermal"            :   -612,
            "anti_kaon_0_thermal"       :   -613,
            "nucleon_thermal"           :   616, # sum(17, 18)
            "proton_thermal"            :   617,
            "neutron_thermal"           :   618,
            "anti_nucleon_thermal"      :   -616, # sum(-17, -18)
            "anti_proton_thermal"       :   -617,
            "anit_neutron_thermal"      :   -618,
        }

    def collectEccentricitiesAndRIntegrals(self, folder, event_id, db, oldStyleStorage=False):
        """
            This function collects initial eccentricities and r-integrals into
            the specified SqliteDB object "db". More specifically,
            this functions fills table "ecc_id_lookup", "eccentricity", and
            "r_integrals".

            Eccentricity and r-integral files will be looked for in "folder" and
            when filling tables the specified "event_id" will be used.

            When "oldStyleStorage" is set to True, another subfolder
            with name "results" will be appended to "folder" which will
            be compatible to the old style storage format.
        """
        # compatibility treatment
        if oldStyleStorage: folder = path.join(folder, "results")
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
        if db.createTableIfNotExists("ecc_id_lookup", (("ecc_id","integer"), ("ecc_type_name","text"))):
            for pattern, ecc_id, ecc_type_name in typeCollections:
                db.insertIntoTable("ecc_id_lookup", (ecc_id, ecc_type_name))

        # next create the eccentricity and r_integrals table, if not existing
        db.createTableIfNotExists("eccentricity", (("event_id","integer"), ("ecc_id", "integer"), ("r_power", "integer"), ("n","integer"), ("ecc_real","real"), ("ecc_imag","real")))
        db.createTableIfNotExists("r_integrals", (("event_id","integer"), ("ecc_id","integer"), ("r_power","integer"), ("r_inte","real")))

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


    def collectFLowsAndMultiplicities_urqmdBinUtilityFormat(self, folder, event_id, db, multiplicityFactor=1.0):
        """
            This function collects integrated and differential flows data
            and multiplicity and spectra data from "folder" into the
            database "db" using event id "event_id". The "multiplityFactor"
            will be multiplied to the number of particles read from file to
            form the multiplicity value.

            This function fills the following table: "pid_lookup",
            "inte_vn", "diff_vn", "multiplicities", "spectra".

            This funtion should only be applied to a folder where flow
            files are generated by the binUtilities module specifically
            for urqmd.
        """
        # collection of file name patterns, pid, and particle name. The file format is determined from the "filename_format.dat" file
        toCollect = {
            "total"         :   "total", # string in filename, particle name
            "pion"          :   "pion",
            "kaon"          :   "kaon",
            "nucleon"       :   "nucleon",
        }
        toCollect_keys = toCollect.keys()
        filePattern = re.compile("([a-zA-z]*)_flow_([a-zA-Z+]*).dat") # filename pattern, the 2nd matched string needs to be among the toCollect.keys() above in order to be considered "matched"; the 1st matched string will either be "integrated" or "differential"
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
        if db.createTableIfNotExists("pid_lookup", (("name","text"), ("pid","integer"))):
            db.insertIntoTable("pid_lookup", list(self.pidDict.items()))

        # next create various tables
        db.createTableIfNotExists("inte_vn", (("event_id","integer"), ("pid","integer"), ("n","integer"), ("vn_real","real"), ("vn_imag","real")))
        db.createTableIfNotExists("diff_vn", (("event_id","integer"), ("pid","integer"), ("pT","real"), ("n","integer"), ("vn_real","real"), ("vn_imag","real")))
        db.createTableIfNotExists("multiplicities", (("event_id","integer"), ("pid","integer"), ("N","real")))
        db.createTableIfNotExists("spectra", (("event_id","integer"), ("pid","integer"), ("pT","real"), ("N","real")))

        # the big loop
        for aFile in listdir(folder): # get all file names
            matchResult = filePattern.match(aFile) # try to match file names
            if not matchResult: continue # not matched!
            flow_type, particle_string_infile = matchResult.groups() # indicated by the file name
            if particle_string_infile not in toCollect_keys: continue # dont know about this particle
            pid = self.pidDict[toCollect[particle_string_infile]] # get pid
            filename = matchResult.group() # get the file to be opened
            flow_table, multiplicity_table = tableChooser[flow_type] # choose tables to write to
            # read the flow file and write results
            for aLine in open(path.join(folder, filename)):
                data = aLine.split()
                if flow_type == "integrated": # for integrated flow and multiplicity; no pT info
                    # write flow table
                    for n in range(1, largest_n):
                        db.insertIntoTable(flow_table,
                                            (event_id, pid, n, float(data[vn_real_cols[n]]), float(data[vn_imag_cols[n]]))
                                        )
                    # write multiplicity table
                    db.insertIntoTable(multiplicity_table,
                                            (event_id, pid, float(data[N_col])*multiplicityFactor)
                                        )
                elif flow_type == "differential": # for differential flow and multiplicity
                    # write flow table
                    for n in range(1, largest_n):
                        db.insertIntoTable(flow_table,
                                            (event_id, pid, float(data[pT_col]), n, float(data[vn_real_cols[n]]), float(data[vn_imag_cols[n]]))
                                        )
                    # write spectra table
                    db.insertIntoTable(multiplicity_table,
                                            (event_id, pid, float(data[pT_col]), float(data[N_col])*multiplicityFactor)
                                        )

        # close connection to commit changes
        db.closeConnection()


    def collectFLowsAndMultiplicities_iSFormat(self, folder, event_id, db, useSubfolder="spectra"):
        """
            This function collects integrated and differential flows data
            and multiplicity and spectra data from "folder" into the
            database "db" using event id "event_id".

            This function fills the following table: "pid_lookup",
            "inte_vn", "diff_vn", "multiplicities", "spectra".

            This funtion should only be applied to a folder where flow
            files are generated by the iS (or iSS with calculate flow
            mode) module as in pure hydro calculations. As such, the
            subfolder name "useSubfolder" will be appended to "folder"
            automatically.
        """
        # add one more sub-directory
        folder = path.join(folder, useSubfolder)

        # collection of file name patterns, pid, and particle name. The file format is determined from the "filename_format.dat" file
        toCollect = {
            "Charged"       :   "charged_hydro", # string in filename, particle name
            "pion_p"        :   "pion_p_hydro",
            "Kaon_p"        :   "kaon_p_hydro",
            "proton"        :   "proton_hydro",
            "thermal_211"   :   "pion_p_thermal",
            "thermal_321"   :   "kaon_p_thermal",
            "thermal_2212"  :   "proton_thermal",
        }
        filename_inte = "%s_integrated_vndata.dat" # filename for integrated flow files, %s is the "string in filename" defined in toCollect
        filename_diff = "%s_vndata.dat" # filename for differential flow files

        # first write the pid_lookup table, makes sure there is only one such table
        if db.createTableIfNotExists("pid_lookup", (("name","text"), ("pid","integer"))):
            db.insertIntoTable("pid_lookup", list(self.pidDict.items()))

        # next create various tables
        db.createTableIfNotExists("inte_vn", (("event_id","integer"), ("pid","integer"), ("n","integer"), ("vn_real","real"), ("vn_imag","real")))
        db.createTableIfNotExists("diff_vn", (("event_id","integer"), ("pid","integer"), ("pT","real"), ("n","integer"), ("vn_real","real"), ("vn_imag","real")))
        db.createTableIfNotExists("multiplicities", (("event_id","integer"), ("pid","integer"), ("N","real")))
        db.createTableIfNotExists("spectra", (("event_id","integer"), ("pid","integer"), ("pT","real"), ("N","real")))

        # the big loop
        for particle_string_infile in toCollect.keys():
            pid = self.pidDict[toCollect[particle_string_infile]]

            # first, differential flow
            particle_filename = path.join(folder, filename_diff % particle_string_infile)
            if path.exists(particle_filename):
                # extract differential flow and spectra information
                diff_flow_block = np.loadtxt(particle_filename)
                largest_n = int(diff_flow_block.shape[1]/3) # should be an integer
                # write flow table
                for aRow in diff_flow_block:
                    for n in range(1, largest_n):
                        db.insertIntoTable("diff_vn",
                            (event_id, pid, aRow[0], n, aRow[3*n], aRow[3*n+1])
                        )
                    # write spectra table
                    db.insertIntoTable("spectra",
                        (event_id, pid, aRow[0], aRow[2])
                    )


            # next, integrated flow
            particle_filename = path.join(folder, filename_inte % particle_string_infile)
            if path.exists(particle_filename):
                # extract integrated flow and multiplicity information
                inte_flow_block = np.loadtxt(particle_filename)
                largest_n = inte_flow_block.shape[0]
                # write flow table
                for n in range(1, largest_n):
                    db.insertIntoTable("inte_vn",
                        (event_id, pid, n, inte_flow_block[n,3], inte_flow_block[n,4])
                    )
                # write multiplicity table
                db.insertIntoTable("multiplicities",
                    (event_id, pid, inte_flow_block[0,1])
                )

        # close connection to commit changes
        db.closeConnection()


    def createDatabaseFromEventFolders(self, folder, subfolderPattern="event-(\d*)", databaseFilename="CollectedResults.db", collectMode="fromUrQMD", multiplicityFactor=1.0):
        """
            This function collect all results (ecc+flow) from subfolders
            whose name have pattern "subfolderPattern" to a database
            with name "databaseFilename".

            The "subfolderPattern" argument can be such that when it
            matches, if its groups()[0] exists, it will be used as event
            id, otherwise the order of the subfolder in listdir will be
            used as event id. Only folders will be matched in either
            case.

            The "collectMode" argument controls how data files to be
            collected are stored internally, and whether oversampling is
            enabled. So far it can be set to either "fromUrQMD" or
            "fromPureHydro":

            -- "fromUrQMD": For eccentricity it will set
            "oldStyleStorage=False" in the
            collectEccentricitiesAndRIntegrals function; for flows
            collectFLowsAndMultiplicities_urqmdBinUtilityFormat will be
            called. In this mode "multiplicityFactor" will be passed
            along to collectFLowsAndMultiplicities_urqmdBinUtilityFormat
            function.

            -- "fromPureHydro": For eccentricity it will set
            "oldStyleStorage=True" in the
            collectEccentricitiesAndRIntegrals function; for flows
            collectFLowsAndMultiplicities_iSFormat will be called.

            -- "fromPureHydroNewStoring": For eccentricity it will set
            "oldStyleStorage=False" in the collectEccentricitiesAndRIntegrals
            function; for flows collectFLowsAndMultiplicities_iSFormat will be
            called with "useSubfolder=''"
        """
        # get list of (matched subfolders, event id)
        matchPattern = re.compile(subfolderPattern)
        matchedSubfolders = []
        for folder_index, aSubfolder in enumerate(listdir(folder)):
            fullPath = path.join(folder, aSubfolder)
            if not path.isdir(fullPath): continue # want only folders, not files
            if collectMode == "fromPureHydro":
                # no matching is needed; try all folders
                matchedSubfolders.append((fullPath, len(matchedSubfolders)+1))
            else:
                # for new-style calculations, folder name must match
                matchResult = matchPattern.match(aSubfolder)
                if matchResult: # matched!
                    if len(matchResult.groups()): # folder name contains id
                        event_id = matchResult.groups()[0]
                    else:
                        event_id = folder_index
                    matchedSubfolders.append((fullPath, event_id)) # matched!

        # the data collection loop
        db = SqliteDB(path.join(folder, databaseFilename))
        if collectMode == "fromUrQMD":
            print("-"*60)
            print("Using fromUrQMD mode")
            print("-"*60)
            for aSubfolder, event_id in matchedSubfolders:
                print("Collecting %s as with event-id: %s" % (aSubfolder, event_id))
                self.collectEccentricitiesAndRIntegrals(aSubfolder, event_id, db) # collect ecc
                self.collectFLowsAndMultiplicities_urqmdBinUtilityFormat(aSubfolder, event_id, db, multiplicityFactor) # collect flow
        elif collectMode == "fromPureHydro":
            print("-"*60)
            print("Using fromPureHydro mode")
            print("-"*60)
            for aSubfolder, event_id in matchedSubfolders:
                print("Collecting %s as with event-id: %s" % (aSubfolder, str(event_id)))
                self.collectEccentricitiesAndRIntegrals(aSubfolder, event_id, db, oldStyleStorage=True) # collect ecc
                self.collectFLowsAndMultiplicities_iSFormat(aSubfolder, event_id, db) # collect flow
        elif collectMode == "fromPureHydroNewStoring":
            print("-"*60)
            print("Using fromPureHydro mode")
            print("-"*60)
            for aSubfolder, event_id in matchedSubfolders:
                print("Collecting %s as with event-id: %s" % (aSubfolder, event_id))
                self.collectEccentricitiesAndRIntegrals(aSubfolder, event_id, db, oldStyleStorage=False) # collect ecc, no subfolders
                self.collectFLowsAndMultiplicities_iSFormat(aSubfolder, event_id, db, useSubfolder="") # collect flow


    def mergeDatabases(self, toDB, fromDB):
        """
            Meger the database "fromDB" to "toDB"; both are assumed to be
            databases created from ebe calculations, meaning that they only
            contain tables specified in EbeCollector_readme.
        """
        for aTable in fromDB.getAllTableNames():
            # first copy table structure
            firstCreation = toDB.createTableIfNotExists(aTable, fromDB.getTableInfo(aTable))
            if firstCreation:
                # just copy
                toDB.insertIntoTable(aTable, fromDB.selectFromTable(aTable))
            else: # treatment depends on table type
                if "lookup" in aTable: continue # if it's a lookup table, nothing to be done
                # not a lookup table: shift up event_id by the current existing max
                currentEventIdMax = toDB.selectFromTable(aTable, "max(event_id)")[0][0]
                def shiftEID(row):
                    newRow = list(row)
                    newRow[0] += currentEventIdMax
                    return newRow
                toDB.insertIntoTable(aTable, list(map(shiftEID, fromDB.selectFromTable(aTable))))
        toDB.closeConnection() # commit


if __name__ == '__main__':
    import doctest
    doctest.testfile("EbeCollector_readme.txt")
