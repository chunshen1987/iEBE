
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
-- event_id (integer primary key)
-- ecc_id (integer). This determines the weight function.
-- r_power (integer). This is the m value.
-- n (integer). This is the harmonic order.
-- ecc_real (real)
-- ecc_imag (real)

3) Table "r_integrals".
-- event_id (integer primary key)
-- r_power (integer). This is the m value.
-- r_inte (real) (integral of r^m)

4) Table "singles". This table is used to store those quantities whose only dependence is "event_id".
-- event_id (integer primary key)
-- dSdy (real)
-- dEdy (real)

The final states of events contains dependence of particle species and various pT values. Since sometimes pT needs to be compared it is better to use interger values so each pT will be associated to a corresponding "pT_id" value, and similarly each particle is associated to a particle identification number (pid), and these two correspondence are recorded in the next two tables.

5) Table "pT_id_lookup". This table records the correspondence between "pT_id" the real "pT" value, so in other tables the presense of "pT" will be replaced by "pT_id".
-- pT_id (integer)
-- pT (real)

6) Table "pid_lookup". This table assoicate each particle to a unique number "pid", which is used to identify particle species in other tables.
-- pid (integer)
-- name (text)

The rest of the tables store information regarding the final states of the events, where the harmonic order is denoted as "n".

7) Table "inte_vn".
-- event_id (interger primary key)
-- pid (integer). This is the particle index value.
-- n (integer)
-- vn_real (real)
-- vn_imag (real)

8) Table "diff_vn".
-- event_id (integer primary key)
-- pid (integer). This is the particle index value.
-- pT_id (integer). This is the pT index value.
-- n (integer)
-- vn_real (real)
-- vn_imag (real)

9) Table "multiplicities".
-- event_id (integer primary key)
-- pid (integer)
-- N (real)

10) Table "spectra".
-- event_id (integer primary key)
-- pid (integer). This is the particle index value.
-- pT_id (integer). This is the pT index value.
-- N (real)


-------------------------------
2. Structure of the package
-------------------------------

The main class is the EbeCollector class, which has several member functions, each used to collect data of certain types (not necessarily corresponding to one table). Each of the collector function accepts an argument specifying the path where the data files are stored. The filenames of the data file that will be collected, as well as the format of these datafiles, are all internal to these functions. The following explains them in details.

1) collectEccentricitiesAndRIntegrals
This function has




The End.
