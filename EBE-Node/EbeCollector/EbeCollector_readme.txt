
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
-- N (real): dN/dy

8) Table "spectra".
-- event_id (integer)
-- pid (integer). This is the particle index value.
-- pT (real)
-- N (real): dN/(dy dpT)

There is an additional table storing scalar quantity for each event.

9) Table "scalars"
-- event_id (integer)
-- lifetime (real)

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

1*) collectScalars(folder, event_id, db)

This function collect scalar info from the folder "folder", then write it into the "scalars" table. So far the lifetime is included.

>>> collector.collectScalars("testData", 1, db)
>>> set(db.selectFromTable("scalars", "lifetime")) == set([(3.34351882,)])
True

2) collectFLowsAndMultiplicities_urqmdBinUtilityFormat(folder, event_id, db, multiplicityFactor)

This function read flow files from the folder "folder", then write the flow and multiplicity results associated to this event to the SqliteDB database "db" using event id "event_id". The multiplicity will be multiplied by the factor "multiplicityFactor" (with oversampling, the counted number of particles is not the actual multiplicity). Acceptable file names should match "([a-zA-z]*)_flow_([a-zA-Z+]*).dat" (e.g. "integrated_flow_Charged.dat"), where the first () can be either "integrated" or "differential", and the second () should be the particle type name.

For example, assuming that the "testData" folder exists (should be included in the package), the following call collect the flow and multiplicity data from it:

>>> collector.collectFLowsAndMultiplicities_urqmdBinUtilityFormat("testData", 1, db, multiplicityFactor=0.1)

We can check the tables it contains by do the following:
>>> set(db.getAllTableNames()) == set([u'ecc_id_lookup', u'eccentricities', u'r_integrals', u'pid_lookup', u'inte_vn', u'diff_vn', u'multiplicities', u'spectra', u'scalars'])
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
[(0.09254586959349598, 82.00000000000001), (0.2230454975609756, 136.66666666666669), (0.3697138657718121, 99.33333333333334), (0.5221993366336635, 67.33333333333334), (0.670072, 30.0), (0.8232357499999998, 16.000000000000004), (0.9717935714285714, 9.333333333333334), (1.081035, 1.3333333333333335), (1.26097, 4.000000000000001), (1.4127966666666665, 2.0000000000000004), (1.5787820000000001, 3.3333333333333335)]


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
>>> reader.get_Ecc_n()
array([-0.00016182-0.00449052j, -0.01112089-0.03613683j,
       -0.13409214+0.09590283j,  0.00110106+0.00146287j])


2) R-integrals.

The getRIntegrals(eccType="ed", r_power=2, where="", orderBy="event_id") function returns a list of r-integral value from the "r_integrals" table for all events, and by default order them by "event_id". Additional requirements can be added by modifying the "where" and "orderBy" arguments. For example:
>>> reader.getRIntegrals()
array([[   51.429047],
       [ 1628.8231  ],
       [  675.07625 ],
       [   45.966488]])

2*) Lifetimes.

The getLifetime(orderBy="event_id") function returns a list of lifetimes for all events, ordered by their "event_id".
>>> reader.getLifetimes()


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
array([  22.8,  290.4,   84. ,   26.6])

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

Similarly to differential flow there is a getInterpretedSpectraForOneEvent(event_id=1, particleName="pion", pTs=np.linspace(0,2.5,10)) function that perform interpolation and another getInterpretedSpectraForAllEvents(particleName="pion", pTs=np.linspace(0,2.5,10), where="", orderBy="event_id") that does so for all events. The later is aliased as get_dNdydpT function, which is recommended way for accessing spectra. For example:
>>> reader.get_dNdydpT(pTs=[0.5, 1.5])
array([[  3.37855688,   0.1       ],
       [ 48.18746303,   1.91488976],
       [ 11.62703839,   0.31319786],
       [  3.81290396,   0.24769716]])


7) Ultimate evaluation function.

The evaluateExpression(expression) function evaluates an expression. The expression will first be converted into a normalized form, then the stardard expression will be converted into function calls, after this the expression will be evaluated. It returns a tuple:
(value, expression after normalization, expression after functionization)

The expressions that it can recognize is listed in the following. Note that all spaces will be elliminated first so any spaces in the original expression will not influence the matching process.

<1> Eccentricity.

The standard form of the complex eccentricity vector is of the form: "Ecc_{m,n}(ed)", which is defined to be { (-1) r^m e^{i n phi} }_e (e is the weight function), and the one using entropy density has the form "Ecc_{m,n}(sd)".

For loose matching the following are some rules:
Eccentricity_, E_ -> Ecc_
(e) -> (ed)
(s) -> (sd)
Ecc_n, Ecc_{n} -> Ecc_{n,n}

A few examples.
For 3rd order complex eccentricity vector using energy density as the weight function:
>>> res = reader.evaluateExpression("Ecc_3(ed)")
>>> "{} -> {}".format(res[1], res[2])
'Ecc_{3,3}(ed) -> self.get_Ecc_n(eccType="ed", r_power=3, order=3)'
>>> res[0]
array([-0.01453327 +5.60409700e-03j,  0.03313195 -5.51618030e-02j,
       -0.15194745 -3.41929410e-01j,  0.00647459 +9.97194020e-05j])

Conventional complex eccentricity vector using entropy density as the weight function:
>>> res = reader.evaluateExpression("E_2 (s)")
>>> "{} -> {}".format(res[1], res[2])
'Ecc_{2,2}(sd) -> self.get_Ecc_n(eccType="sd", r_power=2, order=2)'
>>> res[0]
array([ 0.00034019-0.00471038j, -0.01500890-0.0461648j ,
       -0.11382886+0.07742261j,  0.00026303+0.00115929j])

The r^3 weighted 1st order eccentricity vector:
>>> res = reader.evaluateExpression("E_{3, 1} (e)")
>>> "{} -> {}".format(res[1], res[2])
'Ecc_{3,1}(ed) -> self.get_Ecc_n(eccType="ed", r_power=3, order=1)'
>>> res[0]
array([-0.20594250+0.03413583j, -0.17964724+0.30506617j,
        0.64069616+0.43195141j, -0.02206065-0.0115441j ])

The "ecc" symbol stands for the magnitude of the "Ecc" vector, and it has the same loose matching rules as stated above. In addition, "|Ecc|" is also recoginzable.

A few examples:
Same as the first example above, but for magnitudes:
>>> res = reader.evaluateExpression("ecc_3(ed)")
>>> "{} -> {}".format(res[1], res[2])
'|Ecc_{3,3}(ed)| -> abs(self.get_Ecc_n(eccType="ed", r_power=3, order=3))'
>>> res[0]
array([ 0.01557632,  0.06434711,  0.37417075,  0.00647536])

Showing the "|Ecc_n|=|E_n|=|E|_n" syntax:
>>> res = reader.evaluateExpression("|E|_{7, 5} (e)")
>>> "{} -> {}".format(res[1], res[2])
'|Ecc_{7,5}(ed)| -> abs(self.get_Ecc_n(eccType="ed", r_power=7, order=5))'
>>> res[0]
array([ 0.02074981,  0.08562445,  0.17040317,  0.03101743])

You can also use "Epsilon" and "epsilon":
>>> res = reader.evaluateExpression("epsilon_{2} (e)")
>>> "{} -> {}".format(res[1], res[2])
'|Ecc_{2,2}(ed)| -> abs(self.get_Ecc_n(eccType="ed", r_power=2, order=2))'
>>> res[0]
array([ 0.00449343,  0.03780932,  0.16485768,  0.00183093])

Replacing "Ecc_{m,n}" by "Phi_{m,n}" gives the corresponding event plane angle lies between +-pi/n (using short axis of the deformation as + direction):
>>> res = reader.evaluateExpression("Phi_2 (e)")
>>> "{} -> {}".format(res[1], res[2])
'$Ecc_{2,2}(ed)$ -> angle(self.get_Ecc_n(eccType="ed", r_power=2, order=2))/2'
>>> res[0]
array([-0.80340868, -0.93467121,  1.26036869,  0.4627947 ])

<2> R-averages.

The standard form of the r-integral is of the form: "{r^m}(ed)", which is defined as int(r^m*ed)/int(ed) (e is the weight function), and the one using entropy density as weight function is "{r^m}(sd)".

For loose matching the following are some rules:
{R^ -> {r^

A few examples:
>>> res = reader.evaluateExpression("{ R^2 }(e)")
>>> "{} -> {}".format(res[1], res[2])
'{r^2}(ed) -> self.getRIntegrals(eccType="ed", r_power=2) / self.getRIntegrals(eccType="ed", r_power=0)'
>>> res[0]
array([[ 1.12722936],
       [ 2.5233512 ],
       [ 4.60676773],
       [ 0.98940458]])

Replacing "{}" by "[]" will give the r-integrals instead of r-averages:
>>> res = reader.evaluateExpression("[r^2](e)")
>>> "{} -> {}".format(res[1], res[2])
'[r^2](ed) -> self.getRIntegrals(eccType="ed", r_power=2)'
>>> res[0]
array([[   51.429047],
       [ 1628.8231  ],
       [  675.07625 ],
       [   45.966488]])

<2*> Lifetimes.

The lifetimes of all the fireballs can be returned simply by using the keyword "lifetime":
>>> reader.evaluateExpression("lifetime")[0]


<3> Integrated flows.

The standard form of the complex integrated flow vector is of the form "V_n(pion)", which is the n-th order integrated complex flow vector for "pion".

A few examples:
>>> res = reader.evaluateExpression("V_2(pion)")
>>> "{} -> {}".format(res[1], res[2])
'V_{2}(pion) -> self.get_V_n(particleName="pion", order=2)'
>>> res[0]
array([ 0.01927809+0.05239437j, -0.00390843+0.05045172j,
        0.07518625+0.02788557j, -0.03141065-0.01320664j])

The particle with name "pion_p" is not in the test database, so the following example gives an empty result:
>>> res = reader.evaluateExpression("V_2(pion_p)")
>>> "{} -> {}".format(res[1], res[2])
'V_{2}(pion_p) -> self.get_V_n(particleName="pion_p", order=2)'
>>> res[0]
array([], dtype=float64)

For magnitudes of the flow use "v_n=|V|_n=|V_n|", for example:
>>> res = reader.evaluateExpression("v_2(pion)")
>>> "{} -> {}".format(res[1], res[2])
'|V_{2}(pion)| -> abs(self.get_V_n(particleName="pion", order=2))'
>>> res[0]
array([ 0.05582844,  0.05060289,  0.08019088,  0.0340741 ])

Replacing "V_n" by "Psi_n" gives the n-th order event plane angle between +- pi/n:
>>> res = reader.evaluateExpression("Psi_2(pion)")
>>> "{} -> {}".format(res[1], res[2])
'$V_{2}(pion)$ -> angle(self.get_V_n(particleName="pion", order=2))/2'
>>> res[0]
array([ 0.6091139 ,  0.82405529,  0.17757978, -1.37179069])

<4> Multiplicities.

The standard form of the multiplicity is of the form "dN/dy(pion)", which is the multiplcity of pion.

For loose matching the following are some rules:
N( -> dN/dy(
dN( -> dN/dy(

A few examples:
>>> res = reader.evaluateExpression("dN/dy(pion)")
>>> "{} -> {}".format(res[1], res[2])
'dN/dy(pion) -> self.get_dNdy(particleName="pion")'
>>> res[0]
array([  22.8,  290.4,   84. ,   26.6])

>>> res = reader.evaluateExpression("dN(total)")
>>> "{} -> {}".format(res[1], res[2])
'dN/dy(total) -> self.get_dNdy(particleName="total")'
>>> res[0]
array([  31.7,  395. ,  116.6,   37.7])

<5> Differential flows.

The standard form of the complex differential flow vector is of the form "V_n(pTs)(pion)", which is the differential flow of pion at given pT value(s).

A few examples.
Differential flow of pion at pT=0.5 GeV:
>>> res = reader.evaluateExpression("V_2(0.5)(pion)")
<BLANKLINE>
Calculating differential flow involves interpolation.
Evaluating it at multiple pT values at the same time if possible.
<BLANKLINE>
For better effeciency part of the database is being copied to memory...
Copy completed.
<BLANKLINE>
Looping over 4 events... (please be patient)
Events processed: 0
Done. Thanks for waiting.
>>> "{} -> {}".format(res[1], res[2])
'V_{2}(0.5)(pion) -> self.get_diff_V_n(particleName="pion", order=2, pTs=0.5, verbose=True)'
>>> res[0]
array([ 0.11968908-0.06324814j, -0.04217610+0.07601679j,
        0.23127808+0.05106764j, -0.21871027-0.25311593j])

Differential flow of total particles at pT=0 and pT=1 GeV:
>>> res = reader.evaluateExpression("V_2(linspace(0,1,2))(total)")
<BLANKLINE>
Calculating differential flow involves interpolation.
Evaluating it at multiple pT values at the same time if possible.
<BLANKLINE>
For better effeciency part of the database is being copied to memory...
Copy completed.
<BLANKLINE>
Looping over 4 events... (please be patient)
Events processed: 0
Done. Thanks for waiting.
>>> "{} -> {}".format(res[1], res[2])
'V_{2}(linspace(0,1,2))(total) -> self.get_diff_V_n(particleName="total", order=2, pTs=linspace(0,1,2), verbose=True)'
>>> res[0]
array([[-0.08786127+0.04011986j,  0.02558784-0.02479659j],
       [ 0.01758154+0.01620859j, -0.14460473+0.17226352j],
       [ 0.05387421-0.04663111j,  0.19836007+0.31369321j],
       [ 0.05789534+0.191375j  , -0.09130198+0.03386241j]])

For magnitudes of the flow use "v_n=|V|_n=|V_n|", for example:
>>> res = reader.evaluateExpression("v_7(0.15)(kaon)")
<BLANKLINE>
Calculating differential flow involves interpolation.
Evaluating it at multiple pT values at the same time if possible.
<BLANKLINE>
For better effeciency part of the database is being copied to memory...
Copy completed.
<BLANKLINE>
Looping over 4 events... (please be patient)
Events processed: 0
Done. Thanks for waiting.
>>> "{} -> {}".format(res[1], res[2])
'|V_{7}(0.15)(kaon)| -> abs(self.get_diff_V_n(particleName="kaon", order=7, pTs=0.15, verbose=True))'
>>> res[0]
array([ 0.63721643,  0.13419388,  0.07407628,  0.0893462 ])

Replacing "V_n" by "Psi_n" gives the n-th order event plane angle between +- pi/n:
>>> res = reader.evaluateExpression("Psi_3(0.15)(pion)")
<BLANKLINE>
Calculating differential flow involves interpolation.
Evaluating it at multiple pT values at the same time if possible.
<BLANKLINE>
For better effeciency part of the database is being copied to memory...
Copy completed.
<BLANKLINE>
Looping over 4 events... (please be patient)
Events processed: 0
Done. Thanks for waiting.
>>> "{} -> {}".format(res[1], res[2])
'$V_{3}(0.15)(pion)$ -> angle(self.get_diff_V_n(particleName="pion", order=3, pTs=0.15, verbose=True))/3'
>>> res[0]
array([-0.34120512,  0.63215072,  0.44368578, -0.32626036])

<6> Spectra.

The standard form of the spectra is of the form "dN/(dydpT)(pTs)(pion)".

For loose matching the following are some rules:
dN/dpT, dN/dydpT, dN/(dy dpT) -> dN/(dydpT)

A few examples.
Pion spectra at pT=0.5 GeV:
>>> res = reader.evaluateExpression("dN/dydpT(0.5)(pion)")
<BLANKLINE>
Calculating spectra involves interpolation.
Evaluating it at multiple pT values at the same time if possible.
<BLANKLINE>
For better effeciency part of the database is being copied to memory...
Copy completed.
<BLANKLINE>
Looping over 4 events... (please be patient)
Events processed: 0
Done. Thanks for waiting.
>>> "{} -> {}".format(res[1], res[2])
'dN/(dydpT)(0.5)(pion) -> self.get_dNdydpT(particleName="pion", pTs=0.5, verbose=True)'
>>> res[0]
array([  3.37855688,  48.18746303,  11.62703839,   3.81290396])

Total spectra at pT=0 and 1 GeV:
>>> res = reader.evaluateExpression("dN/dydpT([0,1])(total)")
<BLANKLINE>
Calculating spectra involves interpolation.
Evaluating it at multiple pT values at the same time if possible.
<BLANKLINE>
For better effeciency part of the database is being copied to memory...
Copy completed.
<BLANKLINE>
Looping over 4 events... (please be patient)
Events processed: 0
Done. Thanks for waiting.
>>> "{} -> {}".format(res[1], res[2])
'dN/(dydpT)([0,1])(total) -> self.get_dNdydpT(particleName="total", pTs=[0,1], verbose=True)'
>>> res[0]
array([[  4.7       ,   1.44748483],
       [ 49.6       ,  18.99128464],
       [ 17.1       ,   2.82785324],
       [  5.8       ,   1.30400747]])

<7> Advanced evaluation features.

All the expression above can be mixed algebraically. Any math function supported by numpy package is supported. In addition, the "|Q|" is the same as "abs(Q)", and "<Q>" is the same as "mean(Q)".

A few examples.

Calculate v_2[2](pion):
>>> res = reader.evaluateExpression("sqrt(< v_2(pion)**2 >)")
>>> "{} -> {}".format(res[1], res[2])
'sqrt(<|V_{2}(pion)|**2>) -> sqrt(mean(abs(self.get_V_n(particleName="pion", order=2))**2,0))'
>>> res[0]
0.05759576300663468

In fact the "[2]" syntex is supported:
>>> res = reader.evaluateExpression("v_2[2](pion)")
>>> "{} -> {}".format(res[1], res[2])
'sqrt(<|V_{2}(pion)|**2>) -> sqrt(mean(abs(self.get_V_n(particleName="pion", order=2))**2,0))'
>>> res[0]
0.05759576300663468

It works for differential flow too:
>>> res = reader.evaluateExpression("v_2[2](0.5)(pion)")
<BLANKLINE>
Calculating differential flow involves interpolation.
Evaluating it at multiple pT values at the same time if possible.
<BLANKLINE>
For better effeciency part of the database is being copied to memory...
Copy completed.
<BLANKLINE>
Looping over 4 events... (please be patient)
Events processed: 0
Done. Thanks for waiting.
>>> "{} -> {}".format(res[1], res[2])
'sqrt(<|V_{2}(0.5)(pion)|**2>) -> sqrt(mean(abs(self.get_diff_V_n(particleName="pion", order=2, pTs=0.5, verbose=True))**2,0))'
>>> res[0]
0.22016044857620587

And also for eccentricity:
>>> res = reader.evaluateExpression("e_2[2](ed)")
>>> "{} -> {}".format(res[1], res[2])
'sqrt(<|Ecc_{2,2}(ed)|**2>) -> sqrt(mean(abs(self.get_Ecc_n(eccType="ed", r_power=2, order=2))**2,0))'
>>> res[0]
0.08460369752054607

Calculate v_2[4](pion):
>>> res = reader.evaluateExpression(" ( 2*<v_2(pion)**2>**2 - <v_2(pion)**4> )**(1.0/4) ")
>>> "{} -> {}".format(res[1], res[2])
'(2*<|V_{2}(pion)|**2>**2-<|V_{2}(pion)|**4>)**(1.0/4) -> (2*mean(abs(self.get_V_n(particleName="pion", order=2))**2,0)**2-mean(abs(self.get_V_n(particleName="pion", order=2))**4,0))**(1.0/4)'
>>> res[0]
0.051918047896339255

Note how they are compared to <v_2(pion)>:
>>> res = reader.evaluateExpression(" <v_2(pion)> ")
>>> "{} -> {}".format(res[1], res[2])
'<|V_{2}(pion)|> -> mean(abs(self.get_V_n(particleName="pion", order=2)),0)'
>>> res[0]
0.055174074938951469

Difference between the event plane and participant plane angles:
>>> res = reader.evaluateExpression(" < | Phi_2(ed) - Psi_2(total) | > ")
>>> "{} -> {}".format(res[1], res[2])
'<|$Ecc_{2,2}(ed)$-$V_{2}(total)$|> -> mean(abs(angle(self.get_Ecc_n(eccType="ed", r_power=2, order=2))/2-angle(self.get_V_n(particleName="total", order=2))/2),0)'
>>> res[0]
1.3872205573450411

>>> res = reader.evaluateExpression(" < | Phi_3(ed) - Psi_3(total) | > ")
>>> "{} -> {}".format(res[1], res[2])
'<|$Ecc_{3,3}(ed)$-$V_{3}(total)$|> -> mean(abs(angle(self.get_Ecc_n(eccType="ed", r_power=3, order=3))/3-angle(self.get_V_n(particleName="total", order=3))/3),0)'
>>> res[0]
0.74089569544563649

<8> A few convenience functions.

The evaluateExpressionOnly function wraps the evaluateExpression function and returns only its o-th component, that is, the result:
>>> reader.evaluateExpressionOnly("dN/dydpT(0.5)(pion)")
<BLANKLINE>
Calculating spectra involves interpolation.
Evaluating it at multiple pT values at the same time if possible.
<BLANKLINE>
For better effeciency part of the database is being copied to memory...
Copy completed.
<BLANKLINE>
Looping over 4 events... (please be patient)
Events processed: 0
Done. Thanks for waiting.
array([  3.37855688,  48.18746303,  11.62703839,   3.81290396])

The getFactoryEvaluateExpressionOnly function return the factory function for evaluateExpressionOnly, which can then be used without referring to the objects. For example:
>>> e = reader.getFactoryEvaluateExpressionOnly()
>>> e("<e_2(e)>")
0.052247841125730748

The getAttendance function returns a list of tuple (name of particle, number of events) for all particles given in the pid table:
>>> reader.getAttendance() # doctest: +ELLIPSIS
[(u'total', 4), (u'charged', 0), (u'pion', 4),...]

The getNumberOfEvents function returns the number of events:
>>> reader.getNumberOfEvents()
4

























































































===========
The END
===========




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
