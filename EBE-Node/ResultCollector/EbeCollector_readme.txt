
================================
    Document for EBER Module
================================

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
-- pid (integer)
-- name (text)

The rest of the tables store information regarding the final states of the events, where the harmonic order is denoted as "n".

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

The main class is the EbeCollector class, which has several member functions, each used to collect data of certain types (not necessarily corresponding to one table). Each of the collector function accepts an argument specifying the path where the data files are stored. The filenames of the data file that will be collected, as well as the format of these datafiles, are all internal to these functions. The following explains them in details.

1) collectEccentricitiesAndRIntegrals(folder, event_id, db)
This function read eccentricity files from the folder "folder", then write the eccentricity and r-integral results associated to this event to the SqliteDB database "db" using event id "event_id".

For example, assuming that the "testData" folder exists (should be included in the package), the following call collect the eccentricity data from it:
>>> import EbeCollector
>>> import DBR
>>> db = DBR.SqliteDB("tmp.db")
>>> collector = EbeCollector.EbeCollector()
>>> collector.collectEccentricitiesAndRIntegrals("testData", 1, db)

We can check the tables it contains by do the following:
>>> db.getAllTableNames()
[u'ecc_id_lookup', u'eccentricity', u'r_integrals']

We can inspect the eccentricity data via different ways, for example, assume we want all those real parts of the eccentricities whose value is bigger than 0.4:
>>> db.selectFromTable("eccentricity", "ecc_real", whereClause="ecc_real>0.4")
[(0.4266273,), (0.41968203,), (0.42908949,), (0.40137673,)]

>> db.deleteDatabase(confirmation=True)
True

The End.
