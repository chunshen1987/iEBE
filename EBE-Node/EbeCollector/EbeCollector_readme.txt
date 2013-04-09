
========================================
    Document for EbeCollector Module
========================================

The EBER module is designed to collect various text data files generated during event-by-event hybrid calculations into a single database file, which contains several tables for different quantities of interest, for all events.

Because this module needs to extract data from various text data files whose formats cannot be learned via other means, this module contains "hard coded" info section (header) instruct the rest of the module how to read different files.

-----------------------------
1. Desgin of the database
-----------------------------

First of all the quantities of interest for each event include the following.
-- eccentricities: {r^m exp(i n theta)}. It has real and imaginary parts and it has m and n dependence. Its definition also depends on whether to use energy or entropy density as the weight function.
-- various radius integrals: integral(r^m). It has m dependence.
-- initial integrated entropy: dS/dy. It is a single quantity.
-- initial integrated energy: dE/dy. It is a single quantity.
-- integrated flows for different species of particles (total charged, thermal and total identified particles) {exp(i n phi)}. It has real and imaginary parts, it has n dependece, and it is different for different species of particles.
-- differential flows for different species of particles at different pTs. It has real and imaginary parts, it has n dependence, it has pT dependence, and it is different for different species of particles.
-- particle multiplicity. It is a single quantity, but it is different for different species of particles.
-- particle momentum spectra. It has pT dependence, and it is different for different species of particles.

Because there are many to many relation between several quantities listed above, the database is structured to have multiple tables, and relations between different quantities are chained via a common field "event_id", whose value uniquely identifies an event. This is the quantity that will be used as the primary key in any table that contains it.

The generated database contains the following tables and fields. Note that since sqlite 3 is dynamically typed, the type of a field (column) is not as strictly declared as in statically typed database system, and there are only three types that are needed: text, integer, and real.

First are a few tables storing quantities related to initial states of events.

1) Table "ecc_id_lookup". This table stores the correspondence between the type of weight function used in eccentricity calculation and a unique ecc_id.
-- ecc_id (integer)
-- ecc_type_name (text)

2) Table "eccentricity".
-- event_id (integer)
-- ecc_id (integer). This determines the weight function.
-- r_power (integer). This is the m value.
-- n (integer). This is the harmonic order.
-- ecc_real (real)
-- ecc_imag (real)

3) Table "r_integrals".
-- event_id (integer)
-- ecc_id (integer). Determines the weight function.
-- r_power (integer). This is the m value.
-- r_inte (real) (integral of r^m)

The final states of events depend on particle species, and we associate the name of each particle to a particle identification number (pid), recorded in the following table:

4) Table "pid_lookup". This table assoicate each particle to a unique number "pid", which is used to identify particle species in other tables.
-- name (text)
-- pid (integer)

The rest of the tables store information regarding the final states of the events, where the harmonic order is denoted as "n". The "pT" field in the following is either the mean pT value for a given bin (for differential quantities), or the mean pT value for all particles (for integrated quantities). Note that the structure of table "inte_vn" and "diff_vn", "multiplicities" and "spectra" are exactly the same.

5) Table "inte_vn".
-- event_id (interger)
-- pid (integer). This is the particle index value.
-- n (integer)
-- vn_real (real)
-- vn_imag (real)

6) Table "diff_vn".
-- event_id (integer)
-- pid (integer). This is the particle index value.
-- pT (real)
-- n (integer)
-- vn_real (real)
-- vn_imag (real)

7) Table "multiplicities".
-- event_id (integer)
-- pid (integer)
-- N (real)

8) Table "spectra".
-- event_id (integer)
-- pid (integer). This is the particle index value.
-- pT (real)
-- N (real)


-------------------------------
2. Structure of the package
-------------------------------

The main class is the EbeCollector class, which has several member functions, each used to collect data of certain types (not necessarily corresponding to one table). Each of the collector function accepts an argument specifying the path where the data files are stored. The filenames of the data file that will be collected, as well as the format of these datafiles, are all hardcoded to these functions. The following explains them in details.


1) collectEccentricitiesAndRIntegrals(folder, event_id, db)

This function read eccentricity files from the folder "folder", then write the eccentricity and r-integral results associated to this event to the SqliteDB database "db" using event id "event_id". The eccentricity files should have names that match either "ecc-init-sd-r_power-(\d*).dat" (entropy-weighted) or "ecc-init-r_power-(\d*).dat" (energy-weighted).

For example, assuming that the "testData" folder exists (should be included in the package), the following call collect the eccentricity data from it:
>>> import EbeCollector
>>> import DBR
>>> db = DBR.SqliteDB("tmp.db")
>>> collector = EbeCollector.EbeCollector()
>>> collector.collectEccentricitiesAndRIntegrals("testData", 1, db)

We can check the tables it contains by do the following:
>>> db.getAllTableNames()
[u'ecc_id_lookup', u'eccentricity', u'r_integrals']

We can inspect the eccentricity data via different ways, for example, assume we want all those real parts of the eccentricities whose value is bigger than 0.6:
>>> set(db.selectFromTable("eccentricity", "ecc_real", whereClause="ecc_real>0.6")) == set([(0.61667342,), (0.60969422,), (0.62655439,), (0.61452488,)])
True

Check the lookup table:
>>> db.selectFromTable("ecc_id_lookup") == [(1, u'sd'), (2, u'ed')]
True


2) collectFLowsAndMultiplicities_urqmdBinUtilityFormat(folder, event_id, db, multiplicityFactor)

This function read flow files from the folder "folder", then write the flow and multiplicity results associated to this event to the SqliteDB database "db" using event id "event_id". The multiplicity will be multiplied by the factor "multiplicityFactor" (with oversampling, the counted number of particles is not the actual multiplicity). Acceptable file names should match "([a-zA-z]*)_flow_([a-zA-Z+]*).dat" (e.g. "integrated_flow_Charged.dat"), where the first () can be either "integrated" or "differential", and the second () should be the particle type name.

For example, assuming that the "testData" folder exists (should be included in the package), the following call collect the flow and multiplicity data from it:

>>> collector.collectFLowsAndMultiplicities_urqmdBinUtilityFormat("testData", 1, db, multiplicityFactor=0.1)

We can check the tables it contains by do the following:
>>> db.getAllTableNames()
[u'ecc_id_lookup', u'eccentricity', u'r_integrals', u'pid_lookup', u'inte_vn', u'diff_vn', u'multiplicities', u'spectra']

>>> ('pion_p_hydro', 307) in db.selectFromTable("pid_lookup")
True

>>> ('charged', 1) in db.selectFromTable("pid_lookup")
True

Get averaged pT and integrated v_3:
>>> db.selectFromTable("inte_vn", ("vn_real", "vn_imag"), whereClause="n=3")
[(-0.010069782011460445, 0.02379270360546618)]

Get multiplicity for total charged particles:
>>> db.selectFromTable("multiplicities", "N", whereClause="pid=1")
[(67.7,)]

Get (pT, real(vn)) table for differential v_2:
>>> db.selectFromTable("diff_vn", ("pT", "vn_real"), whereClause="n=2")
[(0.09254586959349598, -0.0865085498282089), (0.2230454975609756, 0.05519660023016902), (0.3697138657718121, 0.08226567254348736), (0.5221993366336635, 0.17146781471001163), (0.670072, 0.006136178587808401), (0.8232357499999998, -0.004128786702717842), (0.9717935714285714, 0.09468668983488557), (1.081035, 0.32262107326181566), (1.26097, 0.35379764449038165), (1.4127966666666665, -0.5539404886975925), (1.5787820000000001, 0.03421922411074727)]


>>> db.selectFromTable("spectra", ("pT", "N"))
[(0.09254586959349598, 12.3), (0.2230454975609756, 20.5), (0.3697138657718121, 14.9), (0.5221993366336635, 10.100000000000001), (0.670072, 4.5), (0.8232357499999998, 2.4000000000000004), (0.9717935714285714, 1.4000000000000001), (1.081035, 0.2), (1.26097, 0.6000000000000001), (1.4127966666666665, 0.30000000000000004), (1.5787820000000001, 0.5)]


3) collectFLowsAndMultiplicities_iSFormat(folder, event_id, db, useSubfolder="spectra")
This function read flow files from the folder "folder/spectra", then write the flow and multiplicity results associated to this event to the SqliteDB database "db" using event id "event_id". This function looks for files that have name "(particle_string_infile)_vndata.dat" for differential flow and spectra data and files that have name "(particle_string_infile)_integrated_vndata.dat" for integrated flow and multiplicity data.

This function provides similar functionality to the collectFLowsAndMultiplicities_urqmdBinUtilityFormat function but for particle data generated from iS (or flow mode of iSS) program directly. Tests for this function is provided in the test section for createDatabaseFromEventFolders.


4) createDatabaseFromEventFolders(folder, subfolderPattern, databaseFilename, collectMode)

This is the high level "entry function" that collects results from "folder"'s subfolder that match the pattern "subfolderPattern", then write them into the datase under "folder" with the name "databaseFilename". The argument "collectMode" controls how data are collected. If data are from a hybrid calculation (hydro+urqmd), set it to "fromUrQMD". If data are from old pure hydrodynamics calculation, set it to "fromPureHydro". If data are generated from new run but for pure hydrodynamics calculation, use "fromPureHydroNewStoring". For details of what exactly this parameter affects see the docstring.

Assuming that the "testData_newStyle" folder exists (should be included in the package), the following call collect the flow and multiplicity data from its two folders and create a database:
>>> collector.createDatabaseFromEventFolders("testData_newStyle", multiplicityFactor=0.1)
------------------------------------------------------------
Using fromUrQMD mode
------------------------------------------------------------
Collecting testData_newStyle/event-2 as with event-id: 2
Collecting testData_newStyle/event-1 as with event-id: 1

The created database file "CollectedResults.db" can be examined in various ways. The following is just a simple peek:
>>> db_tmp = DBR.SqliteDB("testData_newStyle/CollectedResults.db")
>>> set(db_tmp.selectFromTable("multiplicities", ("event_id", "N"), whereClause="pid=1")) == set([(1, 286.7), (2, 67.7)])
True

Assuming that the "testData_oldStyle" folder exists (should be included in the package), the following call collect the flow and multiplicity data from its two folders and create a database:
>>> collector.createDatabaseFromEventFolders("testData_oldStyle", collectMode="fromPureHydro")
------------------------------------------------------------
Using fromPureHydro mode
------------------------------------------------------------
Collecting testData_oldStyle/5-9-VISH2p1V1.9.0.e-IINIT=2-IEOS=7-iEin=1-iLS=130-T0=0.6-Edec=0.18-vis=0.08-factor=1.0--2012-07-24@07:52:19 as with event-id: 1
Collecting testData_oldStyle/5-98-VISH2p1V1.9.0.e-IINIT=2-IEOS=7-iEin=1-iLS=130-T0=0.6-Edec=0.18-vis=0.08-factor=1.0--2012-07-24@09:59:06 as with event-id: 2

The created database file "CollectedResults.db" can be examined in various ways. The following is just a simple peek:
>>> db_tmp = DBR.SqliteDB("testData_oldStyle/CollectedResults.db")
>>> set(db_tmp.selectFromTable("multiplicities", ("event_id", "N"), whereClause="pid=301")) == set([(1, 1904.27097), (2, 1843.95785)])
True

Assuming that the "testData_PureHydroNewStyle" folder exists (should be included in the package), the following call collect the flow and multiplicity data from its two folders and create a database:
>>> collector.createDatabaseFromEventFolders("testData_PureHydroNewStyle", collectMode="fromPureHydroNewStoring")
------------------------------------------------------------
Using fromPureHydro mode
------------------------------------------------------------
Collecting testData_PureHydroNewStyle/event-2 as with event-id: 2
Collecting testData_PureHydroNewStyle/event-1 as with event-id: 1

The created database file "CollectedResults.db" can be examined in various ways. The following is just a simple peek:
>>> db_tmp = DBR.SqliteDB("testData_PureHydroNewStyle/CollectedResults.db")
>>> set(db_tmp.selectFromTable("multiplicities", ("event_id", "N"), whereClause="pid=301")) == set([(1, 1904.27097), (2, 1843.95785)])
True


5) mergeDatabases(toDatabase, fromDatabase)

This function merges the database "fromDatabase" into "toDatabase". The rule is that is a table is a lookup table (name contains "lookup"), then it is copied only if it does not exist in the target database already; otherwise the table must have a field called "event_id" and this field will be shifted up by the previous max value before merging.

For example, we first create another copy of the database using data under testData_newStyle:
>>> from shutil import copy
>>> copy("testData_newStyle/CollectedResults.db", "testData_newStyle/CollectedResults_copy.db")

Then we merge them:
>>> collector.mergeDatabases(DBR.SqliteDB("testData_newStyle/CollectedResults.db"), DBR.SqliteDB("testData_newStyle/CollectedResults_copy.db"))

Then examine it:
>>> set(DBR.SqliteDB("testData_newStyle/CollectedResults.db").selectFromTable("multiplicities", ("event_id", "N"))) == set([(2, 67.7), (1, 286.7), (4, 67.7), (3, 286.7)])
True






















-------------
Clean ups
-------------

>>> db.deleteDatabase(confirmation=True)
True
>>> DBR.SqliteDB("testData_newStyle/CollectedResults.db").deleteDatabase(confirmation=True)
True
>>> DBR.SqliteDB("testData_newStyle/CollectedResults_copy.db").deleteDatabase(confirmation=True)
True
>>> DBR.SqliteDB("testData_oldStyle/CollectedResults.db").deleteDatabase(confirmation=True)
True
>>> DBR.SqliteDB("testData_PureHydroNewStyle/CollectedResults.db").deleteDatabase(confirmation=True)
True

The End.
