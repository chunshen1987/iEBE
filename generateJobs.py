#! /usr/bin/env python
"""
    This script duplicates the EBE-Node folder and generate a collection of pbs
    files to be batch-submitted. For efficiency all codes inside EBE-Node should
    be compiled.
"""


from sys import argv, exit
from os import makedirs, path
from shutil import copytree, copy, rmtree

# check argv
try:
    # set parameters
    numberOfJobs = int(argv[1])
    numberOfEventsPerJob = int(argv[2])
except:
    print('Usage: generateJobs.py number_of_jobs number_of_events_per_job [results_folder="./RESULTS"] [walltime="03:00:00"] [working_folder="./PlayGround"] [compress_results_folder="yes"]')
    exit()

# set optional parameters
if len(argv)>=4: # folder to store results
    resultsFolder = path.abspath(argv[3])
else:
    resultsFolder = path.abspath("./RESULTS")

if len(argv)>=5: # set wall time
    walltime = argv[4]
else:
    walltime = "30:00:00"

if len(argv)>=6: # set working folder
    workingFolder = path.abspath(argv[5])
else:
    workingFolder = path.abspath("./PlayGround")

if len(argv)>=7: # whether to compress final results folder
    compressResultsFolderAnswer = argv[6]
else:
    compressResultsFolderAnswer = "yes"

# prepare directories
if not path.exists(resultsFolder): makedirs(resultsFolder)
if path.exists(workingFolder): rmtree(workingFolder)
makedirs(workingFolder)

# copy parameter file into the entrance folder
copy("./ParameterDict.py", "./EBE-Node/entrance")

# backup parameter files to the result folder
copy("./EBE-Node/entrance/SequentialEventDriver.py", resultsFolder)
copy("./EBE-Node/entrance/ParameterDict.py", resultsFolder)

# duplicate EBE-Node folder to working directory, write .pbs file
for i in range(1, numberOfJobs+1):
    targetWorkingFolder = path.join(workingFolder, "job-%d" % i)
    # copy folder
    copytree("./EBE-Node", targetWorkingFolder)
    open(path.join(targetWorkingFolder, "job-%d.pbs" % i), "w").write(
"""
#!/usr/bin/env bash
#PBS -N iEBE-%d
#PBS -l walltime=%s
#PBS -j oe
#PBS -S /bin/bash
(cd entrance
    python ./SequentialEventDriver_shell.py %d | tee RunRecord.txt
    mv RunRecord.txt ../finalResults/
)
mv ./finalResults %s/job-%d
""" % (i, walltime, numberOfEventsPerJob, resultsFolder, i)
    )
    if compressResultsFolderAnswer == "yes":
        open(path.join(targetWorkingFolder, "job-%d.pbs" % i), "a").write(
"""
(cd %s
zip -r -m job-%d.zip job-%d
)
""" % (resultsFolder, i, i)
        )

print("Jobs generated. Submit them using submitJobs scripts.")
