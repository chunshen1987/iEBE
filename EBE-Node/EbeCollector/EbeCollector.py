#!/usr/bin/env python2
"""
    This module consists of functions dealing with the collection event-by-event
    results into databases.

"""

from os import path, listdir
import re
import numpy as np
from numpy import * # used by EbeDBReader.evaluateExpression function to support all math operations
from DBR import SqliteDB
from assignmentFormat import assignmentExprStream2IndexDict
from ListRNew import isIterable
from StringSubstitution import StringSubstitution





class EbeCollector(object):
    """
        This class contains functions that collect results from event-by-event
        calculations into databases. For the structure of the database see the documentation in the EbeCollector_readme.txt.
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
            "sigma"             :   21, # sum(22, 23, 24)
            "sigma_p"           :   22,
            "sigma_0"           :   23,
            "sigma_m"           :   24,
            "anti_sigma"        :   -21,
            "anti_simga_p"      :   -22,
            "anti_sigma_0"      :   -23,
            "anti_simga_m"      :   -24,
            "xi"                :   26, # sum(27, 28)
            "xi_0"              :   27,
            "xi_m"              :   28,
            "anti_xi"           :   -26,
            "anti_xi_0"         :   -27,
            "anti_xi_m"         :   -28,
            "lambda"            :   31,
            "anti_lambda"       :   -31,
            "omega"             :   36,
            "anti_omega"        :   -36,
            "phi"               :   41,
        }

        for aParticle in self.pidDict.keys():
            if self.pidDict[aParticle]>=0:
                self.pidDict[aParticle+"_hydro"] = self.pidDict[aParticle]+1000
            else:
                self.pidDict[aParticle+"_hydro"] = self.pidDict[aParticle]-1000
            if self.pidDict[aParticle]>=0:
                self.pidDict[aParticle+"_thermal"] = self.pidDict[aParticle]+2000
            else:
                self.pidDict[aParticle+"_thermal"] = self.pidDict[aParticle]-2000


    def collectEccentricitiesAndRIntegrals(self, folder, event_id, db, oldStyleStorage=False):
        """
            This function collects initial eccentricities and r-integrals into
            the specified SqliteDB object "db". More specifically,
            this functions fills table "ecc_id_lookup", "eccentricities", and
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
        db.createTableIfNotExists("eccentricities", (("event_id","integer"), ("ecc_id", "integer"), ("r_power", "integer"), ("n","integer"), ("ecc_real","real"), ("ecc_imag","real")))
        db.createTableIfNotExists("r_integrals", (("event_id","integer"), ("ecc_id","integer"), ("r_power","integer"), ("r_inte","real")))

        # the big loop
        for aFile in listdir(folder): # get all file names
            for pattern, ecc_id, ecc_type_name in typeCollections: # loop over ecc types
                matchResult = pattern.match(aFile) # try to match file names
                if not matchResult: continue # not matched!
                filename = matchResult.group()
                r_power = matchResult.groups()[0] # indicated by the file name
                # read the eccentricity file and write database
                for idx, aLine in enumerate(open(path.join(folder, filename))): # row index is "n"
                    n = idx+1
                    data = aLine.split()
                    # insert into eccentricity table
                    db.insertIntoTable("eccentricities",
                                        (event_id, ecc_id, r_power, n, float(data[ecc_real_col]), float(data[ecc_imag_col]))
                                    )
                    # insert into r-integrals table but only once
                    if n==1:
                        db.insertIntoTable("r_integrals",
                                            (event_id, ecc_id, r_power, float(data[r_inte_col]))
                                        )

        # close connection to commit changes
        db.closeConnection()


    def collectScalars(self, folder, event_id, db):
        """
            This function collects scalar info and into the "scalars" table.
            The supported scalars include: lifetime of the fireball.
        """
        # first write the scalar, makes sure there is only one such table
        db.createTableIfNotExists("scalars", (("event_id","integer"), ("lifetime","real")))
        # for lifetime
        maxLifetime = np.max(np.loadtxt(path.join(folder, "surface.dat"))[:,1])
        db.insertIntoTable("scalars", (event_id, maxLifetime))
        # for others (future)


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
            "sigma"         :   "sigma",
            "xi"            :   "xi",
            "lambda"        :   "lambda",
            "omega"         :   "omega",
            "phi"           :   "phi",
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
        pT_bins = np.loadtxt(path.join(folder, "pT_bins.dat"))
        dpT = pT_bins[1]-pT_bins[0]

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
                                            (event_id, pid, float(data[pT_col]), float(data[N_col])*multiplicityFactor/dpT)
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
            "Sigma_p"       :   "sigma_p_hydro",
            "Xi_m"          :   "xi_m_hydro",
            "Omega"         :   "omega_hydro",
            "Lambda"        :   "lambda_hydro",
            "Phi"           :   "phi_hydro",
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
                        (event_id, pid, aRow[0], aRow[2]*(2*np.pi)*aRow[0])
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
                self.collectScalars(aSubfolder, event_id, db)  # collect scalars
                self.collectFLowsAndMultiplicities_urqmdBinUtilityFormat(aSubfolder, event_id, db, multiplicityFactor) # collect flow
        elif collectMode == "fromPureHydro":
            print("-"*60)
            print("Using fromPureHydro mode")
            print("-"*60)
            for aSubfolder, event_id in matchedSubfolders:
                print("Collecting %s as with event-id: %s" % (aSubfolder, str(event_id)))
                self.collectEccentricitiesAndRIntegrals(aSubfolder, event_id, db, oldStyleStorage=True) # collect ecc
                self.collectScalars(path.join(aSubfolder,"results"), event_id, db)  # collect scalars
                self.collectFLowsAndMultiplicities_iSFormat(aSubfolder, event_id, db) # collect flow
        elif collectMode == "fromPureHydroNewStoring":
            print("-"*60)
            print("Using fromPureHydroNewStoring mode")
            print("-"*60)
            for aSubfolder, event_id in matchedSubfolders:
                print("Collecting %s as with event-id: %s" % (aSubfolder, event_id))
                self.collectEccentricitiesAndRIntegrals(aSubfolder, event_id, db, oldStyleStorage=False) # collect ecc, no subfolders
                self.collectScalars(aSubfolder, event_id, db)  # collect scalars
                self.collectFLowsAndMultiplicities_iSFormat(aSubfolder, event_id, db, useSubfolder="") # collect flow
        else:
            print("!"*60)
            print("Mode string not found")
            print("!"*60)


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

















class EbeDBReader(object):
    """
        This class is used to help reading database generated by the
        EbeCollector class.
        The database is assumed to have the exact structure as explained in the documentation of the EbeCollector class.
    """
    def __init__(self, database):
        """
            Register a SqliteDB database; set first-use flags.
        """
        # setup database
        if isinstance(database, str):
            if path.exists(database):
                database = SqliteDB(database)
            else:
                raise ValueError("EbeDBReader.__init__: the input argument must be an existing database file.")
        if isinstance(database, SqliteDB):
            self.db = database
        else:
            raise TypeError("EbeDBReader.__init__: the input argument must be a string or a SqliteDB database.")
        # setup lookup tables
        self.ecc_lookup = dict((item[1], item[0]) for item in self.db.selectFromTable("ecc_id_lookup"))
        self.pid_lookup = dict(self.db.selectFromTable("pid_lookup"))

        # set self.hasInitializedStringSubstitution to none for lazy initialization in evaluateExpression function
        self.hasInitializedStringSubstitution = False

    def _ecc_id(self, ecc_type_name):
        """
            Return "ecc_id" from "ecc_type_name".
        """
        return self.ecc_lookup[ecc_type_name]

    def _pid(self, name):
        """
            Return "pid" from particle "name".
        """
        return self.pid_lookup[name]

    def getEccentricities(self, eccType="ed", r_power=2, order=2, where="", orderBy="event_id"):
        """
            Return (real, imag) list for eccentricities for a given "eccType",
            "r_power", and for the given harmonic "order" from the
            eccentricities table with additional criteria given by the "where"
            string argument and "orderBy" string argument.

            -- eccType: the type of eccentricity; it is the "ecc_type_name"
                field in the "ecc_id_lookup" table.
            -- r_power: power of r in the weight function.
            -- order: the harmonic order "n".
            -- where: the "where" clause.
            -- orderBy: the "order by" clause.
        """
        whereClause = "ecc_id=%d and r_power=%d and n=%d" % (self._ecc_id(eccType), r_power, order)
        if where:
            whereClause += " and " + where
        return np.asarray(self.db.selectFromTable("eccentricities", ("ecc_real, ecc_imag"), whereClause=whereClause, orderByClause=orderBy))

    def get_Ecc_n(self, eccType="ed", r_power=2, order=2, where="", orderBy="event_id"):
        """
            Return the complex eccentricity vector from the getEccentricities
            function.
        """
        eccArray = self.getEccentricities(eccType=eccType, r_power=r_power, order=order, orderBy=orderBy)
        return eccArray[:,0] + 1j*eccArray[:,1]

    def getRIntegrals(self, eccType="ed", r_power=2, where="", orderBy="event_id"):
        """
            Return a list of the r-integrals for a given "eccType" for the given
            "r_power". Additional criteria can be provided by the "where" and
            "orderBy" string arguments.

            -- eccType: the type of weight function; it is the "ecc_type_name"
                field in the "ecc_id_lookup" table.
            -- r_power: power of r in the weight function.
            -- where: the "where" clause.
            -- orderBy: the "order by" clause.
        """
        whereClause = "ecc_id=%d and r_power=%d" % (self._ecc_id(eccType), r_power)
        if where:
            whereClause += " and " + where
        return np.asarray(self.db.selectFromTable("r_integrals", "r_inte", whereClause=whereClause, orderByClause=orderBy))

    def getLifetimes(self, orderBy="event_id"):
        """
            Return a list of lifetimes.

            -- orderBy: the "order by" clause.
        """
        return np.asarray(self.db.selectFromTable("scalars", "lifetime", orderByClause=orderBy))

    def getIntegratedFlows(self, particleName="pion", order=2, where="", orderBy="event_id"):
        """
            Return (real, imag) list for integrated flows for the species of
            particle with name "particleName", for the given harmonic "order"
            from the integrated flow table with additional criteria given by the
            "where" string argument and "orderBy" string argument.

            -- particleName: name of particle; "name" field in the "pid_lookup"
                table.
            -- order: the harmonic order "n".
            -- where: the "where" clause.
            -- orderBy: the "order by" clause.
        """
        whereClause = "pid=%d and n=%d" % (self._pid(particleName), order)
        if where:
            whereClause += " and " + where
        return np.asarray(self.db.selectFromTable("inte_vn", ("vn_real, vn_imag"), whereClause=whereClause, orderByClause=orderBy))

    def get_V_n(self, particleName="pion", order=2, where="", orderBy="event_id"):
        """
            Return the complex V_n vector from the getIntegratedFlows function.
        """
        VnArray = self.getIntegratedFlows(particleName=particleName, order=order, where=where, orderBy=orderBy)
        if VnArray.shape[0]:
            return VnArray[:,0] + 1j*VnArray[:,1]
        else:
            return VnArray

    def getMultiplicities(self, particleName="pion", where="", orderBy="event_id"):
        """
            Return the multiplicities for the particle with name "particleName".
            Additional criteria can be added by the "where" and "orderBy"
            arguments.
        """
        whereClause = "pid=%d" % self._pid(particleName)
        if where:
            whereClause += " and " + where
        tmp = np.asarray(self.db.selectFromTable("multiplicities", "N", whereClause=whereClause, orderByClause=orderBy))
        return tmp.reshape(tmp.size)

    get_dNdy = getMultiplicities

    def getDifferentialFlowDataForOneEvent(self, event_id=1, particleName="pion", order=2, pT_range=None, where="", orderBy="pT"):
        """
            Return the (p_T, real(v_n), imag(v_n)) list for the differential
            flow of order "order" for event with id "event_id", for particle
            with name "particleName". Only those data inside the given
            "pT_range" will be returned, otherwise only those satisfying
            pT_range(0)<=pT<=pT_range(1) will be returned.
        """
        whereClause = "event_id=%d and pid=%d and n=%d" % (event_id, self._pid(particleName), order)
        if pT_range:
            whereClause += " and %g<=pT and pT<=%g" % (pT_range[0], pT_range[1])
        if where:
            whereClause += " and " + where
        return np.asarray(self.db.selectFromTable("diff_vn", ("pT", "vn_real", "vn_imag"), whereClause=whereClause, orderByClause=orderBy))

    def getInterpretedComplexDifferentialFlowForOneEvent(self, event_id=1, particleName="pion", order=2, pTs=np.linspace(0,2.5,10)):
        """
            Return the interpreted complex differential flow for event with
            id="event_id" on pT points pTs, for order="order" and event
            id="event_id", and for particle name="particleName". The argument
            pTs must be iterable and it will not be checked.
        """
        diffVnData = self.getDifferentialFlowDataForOneEvent(event_id=event_id, particleName=particleName, order=order)
        return np.interp(pTs, diffVnData[:,0], diffVnData[:,1]) + 1j*np.interp(pTs, diffVnData[:,0], diffVnData[:,2])

    def getInterpretedComplexDifferentialFlowsForAllEvents(self, particleName="pion", order=2, pTs=np.linspace(0,2.5,10), where="", orderBy="event_id", verbose=False):
        """
            Return the interpreted values of complex differential flow for all
            events on pT points pTs, for order="order" and event id="event_id",
            and for particle name="particleName". The argument pTs must be
            iterable and it will be checked. Additional criteria can be added
            with "where" and "orderBy" arguments. Returned value will be a numpy
            matrix so that each row is a differential flow vector for an event.
        """
        # create a buffer in memory
        whereClause = "pid=%d and n=%d" % (self._pid(particleName), order)
        if verbose: print("""
Calculating differential flow involves interpolation.
Evaluating it at multiple pT values at the same time if possible.

For better effeciency part of the database is being copied to memory...""")
        databaseBuffer = SqliteDB(":memory:")
        databaseBuffer.createTableIfNotExists("diff_vn", self.db.getTableInfo("diff_vn"))
        databaseBuffer.insertIntoTable("diff_vn", self.db.selectFromTable("diff_vn", "*", whereClause=whereClause))
        if verbose: print("Copy completed.\n")
        # swap memory and the actual database
        self.oldDb = self.db
        self.db = databaseBuffer

        # perform actions
        #if not isIterable(pTs): pTs = [pTs]
        event_ids = self.db.selectFromTable("diff_vn", "event_id", whereClause=where, groupByClause="event_id", orderByClause=orderBy)
        collectedResults = []
        if verbose: print("Looping over {} events... (please be patient)".format(len(event_ids)))
        count = 0
        for event_id_tuple in event_ids:

            if verbose and count % 100 == 0: print("Events processed: {}".format(count))
            count += 1

            collectedResults.append(self.getInterpretedComplexDifferentialFlowForOneEvent(event_id=event_id_tuple[0], particleName=particleName, order=order, pTs=pTs))
        if verbose: print("Done. Thanks for waiting.")

        # swap back the actual database
        self.db = self.oldDb

        # return results
        return np.asarray(collectedResults)

    get_diff_V_n = getInterpretedComplexDifferentialFlowsForAllEvents

    def getSpectraDataForOneEvent(self, event_id=1, particleName="pion", pT_range=None, where="", orderBy="pT"):
        """
            Return the (p_T, dN/dy) spectra list for event with id "event_id",
            for particle with name "particleName". Only those data inside the
            given "pT_range" will be returned, otherwise only those satisfying
            pT_range(0)<=pT<=pT_range(1) will be returned.
        """
        whereClause = "event_id=%d and pid=%d" % (event_id, self._pid(particleName))
        if pT_range:
            whereClause += " and %g<=pT and pT<=%g" % (pT_range[0], pT_range[1])
        if where:
            whereClause += " and " + where
        return np.asarray(self.db.selectFromTable("spectra", ("pT", "N"), whereClause=whereClause, orderByClause=orderBy))

    def getInterpretedSpectraForOneEvent(self, event_id=1, particleName="pion", pTs=np.linspace(0,2.5,10)):
        """
            Return the interpreted spectra values for event with id="event_id"
            on pT points pTs, for event id="event_id" and for particle
            name="particleName". The argument pTs must be iterable and it will
            not be checked.
        """
        diffVnData = self.getSpectraDataForOneEvent(event_id=event_id, particleName=particleName)
        return np.interp(pTs, diffVnData[:,0], diffVnData[:,1])

    def getInterpretedSpectraForAllEvents(self, particleName="pion", pTs=np.linspace(0,2.5,10), where="", orderBy="event_id", verbose=False):
        """
            Return the interpreted spectra for all events on pT points pTs, for
            event id="event_id", and for particle name="particleName". The
            argument pTs must be iterable and it will be checked. Additional
            criteria can be added with "where" and "orderBy" arguments. Returned
            value will be a numpy matrix so that each row is a spectra vector
            for an event.
        """
        # create a buffer in memory
        whereClause = "pid=%d" % self._pid(particleName)
        if verbose: print("""
Calculating spectra involves interpolation.
Evaluating it at multiple pT values at the same time if possible.

For better effeciency part of the database is being copied to memory...""")
        databaseBuffer = SqliteDB(":memory:")
        databaseBuffer.createTableIfNotExists("spectra", self.db.getTableInfo("spectra"))
        databaseBuffer.insertIntoTable("spectra", self.db.selectFromTable("spectra", "*", whereClause=whereClause))
        if verbose: print("Copy completed.\n")
        # swap memory and the actual database
        self.oldDb = self.db
        self.db = databaseBuffer

        # processing
        #if not isIterable(pTs): pTs = [pTs]
        event_ids = self.db.selectFromTable("spectra", "event_id", whereClause=where, groupByClause="event_id", orderByClause=orderBy)
        collectedResults = []
        if verbose: print("Looping over {} events... (please be patient)".format(len(event_ids)))
        count = 0
        for event_id_tuple in event_ids:

            if verbose and count % 100 == 0: print("Events processed: {}".format(count))
            count += 1

            collectedResults.append(self.getInterpretedSpectraForOneEvent(event_id=event_id_tuple[0], particleName=particleName, pTs=pTs))

        if verbose: print("Done. Thanks for waiting.")

        # swap back the actual database
        self.db = self.oldDb

        # return results
        return np.asarray(collectedResults)


    get_dNdydpT = getInterpretedSpectraForAllEvents

    def getAttendance(self):
        """
            Return a list of tuple (name for particle, number of events) for all
            particles in the pid table in the order of increase |pid|. Harmonic
            order 2 will be used to probe all particles.
        """
        probes = []
        allParticles = self.pid_lookup.items()
        allParticles.sort(key=lambda x: abs(x[1]))
        for aParticle, pid in allParticles:
            numberOfEvents = len(self.get_V_n(particleName=aParticle))
            probes.append((aParticle, numberOfEvents))
        return probes

    def getNumberOfEvents(self):
        """
            Return total number of events by finding the difference between max
            and min of event_id.
        """
        maxEventId = self.db.selectFromTable("eccentricities", "max(event_id)")[0][0]
        minEventId = self.db.selectFromTable("eccentricities", "min(event_id)")[0][0]
        return maxEventId - minEventId + 1

    def evaluateExpression(self, expression):
        """
            Evaluate an expression by first applying substitution rules using
            the StringSubstitution.

            The passed expression will first be converted in to a normalized
            form, then the normalized expression will be replaced by
            corresponding function calls, after this the functionized expression
            will be evaluated.

            It returns the typle
            (value of the expression, string after normalization, string after functionization)

        """
        # remove spaces
        expression = expression.replace(" ", "")
        # perform lazy initialization
        if not self.hasInitializedStringSubstitution:
            
            # The groups of substitution rules it contains will loop until there
            # is no more changes, thus only the relative order between the
            # groups matter: make sure groups appear earlier contain expansions
            # that should be done before groups appear later.
            # Note that all the substitution strings contain no spaces
            self.useStringSubstitution_normalization = (
                
            # 0th priorities: standardize notations
            StringSubstitution((
                
                # common
                ("\(e\)", "(ed)"), # ed -> e
                ("\(s\)", "(sd)"), # sd -> s
                
                # add {} to subscript to enable expansion of [2] and [4]
                ("_([\d]+)", "_{{{0[0]}}}"), # add { } to subscripts
                
                # eccentricities
                ("Eccentricity_", "Ecc_"), # Eccentricity_ -> Ecc_
                ("E_", "Ecc_"), # E_ -> Ecc_
                ("eccentricity_", "ecc_"), # eccentricity_ -> ecc_
                ("e_", "ecc_"), # e_ -> ecc_
                # latex style support
                ("Epsilon_", "Ecc_"),
                ("epsilon_", "ecc_"),
                
                # eccentricity:
                # Ecc_{m,n}(ed) := {r^m e^{i n phi}}_e
                ("Ecc_{([\d]+)}", "Ecc_{{{0[0]},{0[0]}}}"), # Ecc_{n} -> Ecc_{n,n}
                
                # r-averages
                # {r^m}(ed) := int(r^m*ed)/int(ed)
                ("{R\^", "{{r^"),
            
                # r-integrals
                # [r^m](ed) := int(r^m*ed)
                ("\[R\^", "[r^"),
                
                # multiplicity:
                # dN/dy(pion) := pion multiplicity
                ("[^d]N\(", "dN/dy("),
                ("dN\(", "dN/dy("),
            
                # spectra:
                # dN/(dydpT)(pTs)(pion) := pion spectra at pTs values
                ("dN/dpT", "dN/(dydpT)"),
                ("dN/dydpT", "dN/(dydpT)"),                
            )),    
                
            # 1st priorities: expanding [2] [4]
            StringSubstitution((

                # support for xxx_{ooo}[2](oxox)
                ("([\w_]+)_{([\d,]+)}\[2\]\(([\w_]+)\)", 'sqrt(<{0[0]}_{{{0[1]}}}({0[2]})**2>)'), # without (pTs)
                ("([\w_]+)_{([\d,]+)}\[2\](\(.*?\))\(([\w_]+)\)", 'sqrt(<{0[0]}_{{{0[1]}}}{0[2]}({0[3]})**2>)'), # with (pTs)
                
                # support for xxx_{ooo}[4](oxox)
                ("([\w_]+)_{([\d,]+)}\[4\]\(([\w_]+)\)", '((2*<{0[0]}_{{{0[1]}}}({0[2]})**2>**2-<{0[0]}_{{{0[1]}}}({0[2]})**4>)**0.25)'), # without (pTs)
                ("([\w_]+)_{([\d,]+)}\[4\](\(.*?\))\(([\w_]+)\)", '((2*<{0[0]}_{{{0[1]}}}{0[2]}({0[3]})**2>**2-<{0[0]}_{{{0[1]}}}{0[2]}({0[3]})**4>)**0.25)'), # with (pTs)
            )),
            
            # 2nd priorities: expand special functions || <> $$ (related: ecc, v, Phi, Psi)
            StringSubstitution((

                # ecc = |Ecc|
                ("ecc_", "|Ecc|_"),
                # v = |V|
                ("v_", "|V|_"),

                # || = abs
                ("\|([\w_]+)\|(.*?)\(([\w_]+)\)", "|{0[0]}{0[1]}({0[2]})|"), # |ooo|xxx(oxox) -> |oooxxx(oxox)|; oxox is a word
                
                # <> = mean
                ("<([\w_]+)>(.*?)\(([\w_]+)\)", "<{0[0]}{0[1]}({0[2]})>"), # <ooo>xxx(oxox) -> <oooxxx(oxox)>; oxox is a word

                # Phi = $Ecc$
                ("Phi_", "$Ecc$_"),
                # Psi = $V$
                ("Psi_", '$V$_'),

                # $$ = get plane angles; only applies to Ecc and V
                ("\$([\w_]+)\$(.*?)\(([\w_]+)\)", "${0[0]}{0[1]}({0[2]})$"), # <ooo>xxx(oxox) -> <oooxxx(oxox)>; oxox is a word
            )),
            
            )

            # convert standardized notations to functions
            self.useStringSubstitution_functionization = StringSubstitution((

                # ||: absolute value
                ("\|(.*?)\|", 'abs({0[0]})'),

                # <>: mean value
                ("<(.*?)>", 'mean({0[0]},0)'),

                # $$: get plane angles; only applies to Ecc (angle(-Ecc_n)/n) and V (angle(V_n)/n)
                ("\$Ecc_{([\d\w+]),([\d\w+])}(.*?)\$", 'angle(Ecc_{{{0[0]},{0[1]}}}{0[2]})/{0[1]}'),
                ("\$V_{([\d\w+])}(.*?)\$", 'angle(V_{{{0[0]}}}{0[1]})/{0[0]}'),
                
                # eccentricity:
                # ecc_{m,n}(ed) := {-r^m e^{i n phi}}_e
                ("Ecc_{([\d]+),([\d]+)}\((\w\w)\)", 'self.get_Ecc_n(eccType="{0[2]}", r_power={0[0]}, order={0[1]})'), # to functions

                # r-averages
                # {r^m}(ed) := int(r^m*ed)/int(ed)
                ("{r\^([\d]+)}\((\w\w)\)", 'self.getRIntegrals(eccType="{0[1]}", r_power={0[0]}) / self.getRIntegrals(eccType="{0[1]}", r_power=0)'),

                # r-integrals
                # [r^m](ed) := int(r^m*ed)
                ("\[r\^([\d]+)\]\((\w\w)\)", 'self.getRIntegrals(eccType="{0[1]}", r_power={0[0]})'),

                # lifetimes
                ("lifetime", 'self.getLifetimes()'),

                # integrated flow:
                # V_{n}(pion) := pion complex flow vector of order n
                ("V_{([\d]+)}\(([\w_]+)\)", 'self.get_V_n(particleName="{0[1]}", order={0[0]})'),

                # multiplicity:
                # dN/dy(pion) := pion multiplicity
                ("dN/dy\(([\w_]+)\)", 'self.get_dNdy(particleName="{0[0]}")'),

                # differential flows
                # V_{n}(pTs)(pion) := complex differential flow vector of order n for pion at pTs values
                ("V_{([\d]+)}\((.*?)\)\(([\w_]+)\)", 'self.get_diff_V_n(particleName="{0[2]}", order={0[0]}, pTs={0[1]}, verbose=True)'),

                # spectra:
                # dN/(dydpT)(pTs)(pion) := pion spectra at pTs values
                ("dN/\(dydpT\)\((.*?)\)\(([\w_]+)\)", 'self.get_dNdydpT(particleName="{0[1]}", pTs={0[0]}, verbose=True)'),

            ))


        # perform normalization, should repeat until there is no more changes
        exprAfterNormalization = expression
        needMoreChanges = True
        while needMoreChanges:
            needMoreChanges = False
            for stringSubstitution in self.useStringSubstitution_normalization:
                exprAfterNormalization, numberOfScans = stringSubstitution.applyAllRules(exprAfterNormalization)
                if numberOfScans>0: needMoreChanges = True
        # perform functionization, should do only once
        exprAfterFunctionization, numberOfScans = self.useStringSubstitution_functionization.applyAllRules(exprAfterNormalization)
        # try to evaluate it
        try:
            value = eval(exprAfterFunctionization)
            return (value, exprAfterNormalization, exprAfterFunctionization)
        except:
            print("Error encounterred evaluating {}:".format(expression))
            print("-> {}\n-> {}".format(exprAfterNormalization, exprAfterFunctionization))
            raise

    def evaluateExpressionOnly(self, expression):
        """
            Wraps evaluateExpression function; returns only the result, not
            expressions for checking.
        """
        try:
            value, expr1, expr2 = self.evaluateExpression(expression)
            return value
        except:
            pass # ignore

    def getFactoryEvaluateExpressionOnly(self):
        """
            Return factory functions for evaluateExpressionOnly.
        """
        # factory function for evaluateExpressionOnly
        def evaluateExpressionOnly_factory(expression):
            return self.evaluateExpressionOnly(expression)

        return evaluateExpressionOnly_factory

if __name__ == '__main__':
    import doctest
    doctest.testfile("EbeCollector_readme.txt")
