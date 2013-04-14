
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

2) Table "eccentricities".
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

4) Table "pid_lookup". This table assoicate each particle to a unique number "pid", which is used to identify particle species in other tables. The particles are categorized into 3 groups: the 1st group contains final particles after UrQMD; the 2nd group contains particles after resonance decay from pure hydrodynamics calculation; and the 3rd group contains thermal particles. Note that depending on the type of calculations, not all of the categories will eventually be filled.
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

The main class is the EbeCollector class, which has several member functions, each used to collect data of certain types (not necessarily corresponding to one table). Each of the collector function accepts an argument specifying the path where the data files are stored. The filenames of the data file that will be collected, as well as the format of these datafiles, are all hardcoded to these functions. Section 3 explains them in details.

Another convenience class included in the package is the EbeDBReader class, which is used to read the database generated by the EbeCollector class. Its interface is explained in section 4.


------------------------------------------
3. Structure of the EbeCollector class
------------------------------------------

1) collectEccentricitiesAndRIntegrals(folder, event_id, db)

This function read eccentricity files from the folder "folder", then write the eccentricity and r-integral results associated to this event to the SqliteDB database "db" using event id "event_id". The eccentricity files should have names that match either "ecc-init-sd-r_power-(\d*).dat" (entropy-weighted) or "ecc-init-r_power-(\d*).dat" (energy-weighted).

For example, assuming that the "testData" folder exists (should be included in the package), the following call collect the eccentricity data from it:
>>> import EbeCollector
>>> import DBR
>>> db = DBR.SqliteDB("tmp.db")
>>> collector = EbeCollector.EbeCollector()
>>> collector.collectEccentricitiesAndRIntegrals("testData", 1, db)

We can check the tables it contains by do the following:
>>> set(db.getAllTableNames()) == set([u'ecc_id_lookup', u'eccentricities', u'r_integrals'])
True

We can inspect the eccentricity data via different ways, for example, assume we want all those real parts of the eccentricities whose value is bigger than 0.6:
>>> set(db.selectFromTable("eccentricities", "ecc_real", whereClause="ecc_real>0.6")) == set([(0.61667342,), (0.60969422,), (0.62655439,), (0.61452488,)])
True

Check the lookup table:
>>> db.selectFromTable("ecc_id_lookup") == [(1, u'sd'), (2, u'ed')]
True


2) collectFLowsAndMultiplicities_urqmdBinUtilityFormat(folder, event_id, db, multiplicityFactor)

This function read flow files from the folder "folder", then write the flow and multiplicity results associated to this event to the SqliteDB database "db" using event id "event_id". The multiplicity will be multiplied by the factor "multiplicityFactor" (with oversampling, the counted number of particles is not the actual multiplicity). Acceptable file names should match "([a-zA-z]*)_flow_([a-zA-Z+]*).dat" (e.g. "integrated_flow_Charged.dat"), where the first () can be either "integrated" or "differential", and the second () should be the particle type name.

For example, assuming that the "testData" folder exists (should be included in the package), the following call collect the flow and multiplicity data from it:

>>> collector.collectFLowsAndMultiplicities_urqmdBinUtilityFormat("testData", 1, db, multiplicityFactor=0.1)

We can check the tables it contains by do the following:
>>> set(db.getAllTableNames()) == set([u'ecc_id_lookup', u'eccentricities', u'r_integrals', u'pid_lookup', u'inte_vn', u'diff_vn', u'multiplicities', u'spectra'])
True

>>> ('pion_p_hydro', 1007) in db.selectFromTable("pid_lookup")
True

>>> ('charged', 1) in db.selectFromTable("pid_lookup")
True

Get averaged pT and integrated v_3:
>>> db.selectFromTable("inte_vn", ("vn_real", "vn_imag"), whereClause="n=3")
[(-0.010069782011460445, 0.02379270360546618)]

Get multiplicity for total charged particles:
>>> db.selectFromTable("multiplicities", "N", whereClause="pid=0")
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
>>> set(db_tmp.selectFromTable("multiplicities", ("event_id", "N"), whereClause="pid=0")) == set([(1, 286.7), (2, 67.7)])
True

Assuming that the "testData_oldStyle" folder exists (should be included in the package), the following call collect the flow and multiplicity data from its two folders and create a database:
>>> collector.createDatabaseFromEventFolders("testData_oldStyle", collectMode="fromPureHydro") # doctest: +ELLIPSIS
------------------------------------------------------------
Using fromPureHydro mode
------------------------------------------------------------
Collecting testData_oldStyle/5-9... as with event-id: ...
Collecting testData_oldStyle/5-9... as with event-id: ...

The created database file "CollectedResults.db" can be examined in various ways. The following is just a simple peek:
>>> db_tmp = DBR.SqliteDB("testData_oldStyle/CollectedResults.db")
>>> set(db_tmp.selectFromTable("multiplicities", "N", whereClause="pid=1001", orderByClause="N"))
set([(1843.95785,), (1904.27097,)])


Assuming that the "testData_PureHydroNewStyle" folder exists (should be included in the package), the following call collect the flow and multiplicity data from its two folders and create a database:
>>> collector.createDatabaseFromEventFolders("testData_PureHydroNewStyle", collectMode="fromPureHydroNewStoring")
------------------------------------------------------------
Using fromPureHydroNewStoring mode
------------------------------------------------------------
Collecting testData_PureHydroNewStyle/event-2 as with event-id: 2
Collecting testData_PureHydroNewStyle/event-1 as with event-id: 1

The created database file "CollectedResults.db" can be examined in various ways. The following is just a simple peek:
>>> db_tmp = DBR.SqliteDB("testData_PureHydroNewStyle/CollectedResults.db")
>>> set(db_tmp.selectFromTable("multiplicities", ("event_id", "N"), whereClause="pid=1001")) == set([(1, 1904.27097), (2, 1843.95785)])
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


------------------------------------------
4. Structure of the EbeDBReader class
------------------------------------------


The constructor takes either a SqliteDB database or a string for a SQLite database filename, and store it internally. All other queries are with repect to this database. Assuming that the "collected.db" file exists under "testDB", the following example establishes such a link:
>>> reader = EbeCollector.EbeDBReader("testDB/collected.db")


1) Eccentricities.

the getEccentricities(eccType="ed", r_power=2, order=2, where="", orderBy="event_id") function returns (real(ecc_n), imag(ecc_n)) tuple from the "eccentricities" table for all events, and by default order them by "event_id". Additional requirements can be added by modifying the "where" argument and the "orderBy" argument. The following are a few examples.

To take all energy density defined 2nd order eccentricity matrix:
>>> reader.getEccentricities()
array([[-0.00016182, -0.00449052],
       [-0.01112089, -0.03613683],
       [-0.13409214,  0.09590283],
       [ 0.00110106,  0.00146287]])

This function has a getEccn variant that returns the complex eccentricity factor instead of the the (real, imag) matrix. For example:
>>> reader.get_ecc_n()
array([-0.00016182-0.00449052j, -0.01112089-0.03613683j,
       -0.13409214+0.09590283j,  0.00110106+0.00146287j])


2) R-integrals.

The getRIntegrals(eccType="ed", r_power=2, where="", orderBy="event_id") function returns a list of r-integral value from the "r_integrals" table for all events, and by default order them by "event_id". Additional requirements can be added by modifying the "where" and "orderBy" arguments. For example:
>>> reader.getRIntegrals()
array([[   51.429047],
       [ 1628.8231  ],
       [  675.07625 ],
       [   45.966488]])


3) Integrated flows.

The getIntegratedFlows(particleName="pion_p", order=2, where="", orderBy="event_id") function returns (real(V_n), imag(V_n)) tuple from the "inte_vn" table for all events, and by default order them by "event_id". Additional requirements can be added by modifying the "where" argument and the "orderBy" argument. The following are a few example.

To take all 2nd order total Pion integrated flow matrix:
>>> reader.getIntegratedFlows()
array([[ 0.01927809,  0.05239437],
       [-0.00390843,  0.05045172],
       [ 0.07518625,  0.02788557],
       [-0.03141065, -0.01320664]])

Similarly this function has a getVn variant returning the complex vector:
>>> reader.get_V_n()
array([ 0.01927809+0.05239437j, -0.00390843+0.05045172j,
        0.07518625+0.02788557j, -0.03141065-0.01320664j])


4) Multiplicities.

The getMultiplicities(particleName="pion", where="", orderBy="event_id") function returns a list of multiplicities from the "multiplicities" table for all events, and by default order by "event_id". Additional requirements can be added by modifying the "where" and "orderBy" arguments. For example:
>>> reader.getMultiplicities()
array([[  22.8],
       [ 290.4],
       [  84. ],
       [  26.6]])

This function is also aliased as the get_dNdy function.


5) Differential flows.

The getDifferentialFlowDataForOneEvent(event_id=1, particleName="pion", order=2, pT_range=None, where="", orderBy="pT") function returns the raw (p_T, real(v_n), imag(v_n)) list for the given event with id "event_id", and for all p_T that lies in the range [pT_range[0], pT_range[1]]. Additional requirements can be added by modifying the "where" and "orderBy" arguments. For example:
>>> reader.getDifferentialFlowDataForOneEvent(pT_range=(1,2))
array([[ 1.07688333,  0.02276404,  0.33534406],
       [ 1.22453   , -0.95678683, -0.29079023],
       [ 1.83739   ,  0.14930991,  0.98879045]])

The getInterpretedComplexDifferentialFlowForOneEvent(event_id=1, particleName="pion", order=2, pTs=np.linspace(0,2.5,10)) function returns the interpreted complex differential flows for pT points in pTs. For example:
>>> reader.getInterpretedComplexDifferentialFlowForOneEvent(pTs=[0.5, 1.5])
array([ 0.11968908-0.06324814j, -0.45961542+0.28435922j])

Use the getInterpretedComplexDifferentialFlowsForAllEvents(self, particleName="pion", order=2, pTs=np.linspace(0,2.5,10), where="", orderBy="event_id") to return the interpreted complex differential flow from all events.
>>> reader.getInterpretedComplexDifferentialFlowsForAllEvents(pTs=[0.5, 1.5])
array([[ 0.11968908-0.06324814j, -0.45961542+0.28435922j],
       [-0.04217610+0.07601679j,  0.05340598+0.10338167j],
       [ 0.23127808+0.05106764j, -0.23767418+0.57741666j],
       [-0.21871027-0.25311593j,  0.08848211+0.16774454j]])

The getInterpretedComplexDifferentialFlowsForAllEvents function is also aliased as get_diff_V_n function, and this alias is the recommended function for accessing differential flows.


6) Spectra.

The getSpectraDataForOneEvent(event_id=1, particleName="pion", pT_range=None, where="", orderBy="pT") function returns the raw (p_T, dN/(dy p_T dp_T dphi)) list for the given event with id "event_id", and for all p_T that lies in the range [pT_range[0], pT_range[1]]. Additional requirements can be added by modifying the "where" and "orderBy" arguments. For example:
>>> reader.getSpectraDataForOneEvent(pT_range=(1,2))
array([[ 1.07688333,  0.3       ],
       [ 1.22453   ,  0.1       ],
       [ 1.83739   ,  0.1       ]])

Similarly to differential flow there is a getInterpretedSpectraForOneEvent(event_id=1, particleName="pion", pTs=np.linspace(0,2.5,10)) function that perform interpolation and another getInterpretedSpectraForAllEvents(particleName="pion", pTs=np.linspace(0,2.5,10), where="", orderBy="event_id") that does so for all events. The later is aliased as get_dNdypTdpT function, which is recommended way for accessing spectra. For example:
>>> reader.get_dNdypTdpT(pTs=[0.5, 1.5])
array([[  3.37855688,   0.1       ],
       [ 48.18746303,   1.91488976],
       [ 11.62703839,   0.31319786],
       [  3.81290396,   0.24769716]])












































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
