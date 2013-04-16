
===================
    iEBE Readme
===================

The iEBE package is a convenience package to automate event-by-event hybrid calculations. It divide calculatons into "jobs", where each job consists multiple "ebe-calculations". Each "ebe-calculation" is a complete hybrid calculation that in execution order performs: heavy-ion event generation (superMC), hydrodynamics simulation (VISHNew), particle emission sampling (iSS), hadron rescattering simulation (osc2u and urqmd), flow calculation (binUtilities), and finally collect important results to databases (EbeCollector). Each "job" runs the given number of "ebe-calculations" sequentially, and "jobs" are run in parallel. The package has utility scripts that can combine the generated SQLite database files from different jobs into one that can be analyzed later.

The main programs are contained in the subfolder "EBE-Node", which is used to perform one job, and which will be duplicated when multiple jobs are desired. The package needs two locations to perform multi-job calculations: one folder is used to store duplications of "EBE-Node" and intermediate results generated during the calculation (refer to as "working folder" in the following), and another folder is used to store final results (refer to as "result folder" in the following). By default the working folder is named as "PlayGround" and the result folder is named as "RESULTS", both in the root directory of the package.

------------------------------------------------------------------
<<1>> How to use the package to perform multi-job calculations
------------------------------------------------------------------

This section explains how to use the highest-level scripts provided by the package to perform event-by-event hybrid calculations. This section should be the only section requires reading for user who are not interested in modifying the package.

*VERY IMPORTANT* Make sure you have Python 2.7+ (or Python 3) installed before proceeding.

In the following all path are relative to the root directory of the package.

Step 1) Generate jobs using the ./generateJobs.py script.

To generate jobs use the generateJobs.py script in the root directory. Most of the runnable scripts in this package provide the feature that if you run it without additional arguments it will print the usage echo, for example:
$ ./generateJobs.py
And you should see the output:
Usage: generateJobs.py number_of_jobs number_of_events_per_job [working_folder="./PlayGround"] [results_folder="./RESULTS"] [walltime="03:00:00" (per event)] [compress_results_folder="yes"]

The echo says that the 1st argument for the script should specify the number of jobs you want to generate; the 2nd argument specifying the number of ebe-calculations for each job; the 3rd argument points to the result folder; the 4th argument specifies the "wall time" (used in torque system, explain later); the 5th argument points to the working folder; and the 6th argument is for whether to compress final results. Except for the first two, all other arguments have default values. The simplest way to generate jobs is just to accept the default values, and as an example to generate 2 jobs, eaching performing 5 ebe-calculations, simply do the following:
$ ./generateJobs.py 2 5

This script will first check required libraries, then compile all the programs if not existing already, the generate the actual folders for jobs.

After you see the echo "Jobs generated.", you should see the working folder "PlayGround" and the result folder "RESULTS" in the root directory if not existing previously.

Step 2) Submit jobs.

The way to submit jobs depends on the system. For a cluster that has "torque" scheduling system (therefore the "qsub" command is availible), submit jobs use the submitJobs_qsub.py script; for a local computation use the submitJobs_local.py script. The difference is that the local computation is only parallelled for the local CPUs and calculation on cluster, via the torque system, will be distributed to multiple nodes. To submit a local calculation simply do (The script knows how to get the location of the working folder automatically):
$ ./submitJobs_local.py

You should see some feedbacks listing the jobs that have been submitted.

Step 3) Checking progress.

Progresses for each job can be checked by the progressReport.py script in the root directory (The script knows how to get the location of the working folder automatically):
$ ./progressReport.py

It will list the current progress for all jobs.

Step 4) Combining databases.

Once all calculations are finished, the generated database files from all events will be combined automatically, and a single file "collected.db" will be generated in the results folder.

---------------------------------------
<<2>> How to analyze generated data
---------------------------------------

The "collected.db" generated from previous steps is the SQLite database file that can be analyzed by any desired means. The recommended way is to use the uhg.py script in the utilities folder. This script not only reads the database, but also performs additional analysis like interpolation along pT, calculation of mean, and etc. It can either be run from command line to evaluate a single expression, or interactively from a shell. To evaluate a single expression run the uhg.py script in the utilities folder:
$ ./uhg.py database_filename "expression to evaluate"

A more convenient way to evaluate multiple expressions as well as perform additional analysis, is to run the uhg.py script interactively. For example:
$ python -ic "from uhg import *"

The interactive mode will also print out a simple help showing recognizable symbols that can be included in the expression.

Another way is to use the databaseQuery.py script to evaluate a single piece of SQL query from command line:
# ./databaseQuery.py database_filename "SQL_query"

A third way is to use the ./unpackDatabase.py script located under /EBE-Node/EbeCollector/ to dump the whole database into separated space-separated text files, each for individual table. For example, running the following under /RESULTS/:
$ ../EBE-Node/EbeCollector/unpackDatabase.py ./collected.db .
will generate several ".dat" files, each containing data for the corresponding type. Each file have a one-line header to indicate what data each column records, and the rest are data separated by spaces.

For more details about the structure of the database and the uhg.py script see /EBE-Node/EbeCollector/EbeCollector_readme.txt.

--------------------------------
<<3>> How to tune parameters
--------------------------------

After you familiarized yourself about how to perform multi-job hybrid calculations, finally the time comes to the question about how to tune parameters for the simulations. The most commonly tuned parameters are in the ParameterDict.py file in the directory, which should be the only file use to direct the simluations. This file will be copied to the result folder for record-keeping purpose when generating jobs.


Enjoy!
