#! /usr/bin/env python
"""
    This utility checks the existence of all zip files in the given directory,
    then calls the combineEbeDatabasesFromZippedResults when all the expected
    files are present in the directory.
"""

from sys import argv, exit, stdout
from time import sleep
from os import path, listdir
from datetime import datetime
from subprocess import call
import re

# get options
try:
    parentFolder = path.abspath(argv[1])
    numberOfZipFiles = int(argv[2])
    # get optional parameters
    if len(argv)>=4:
        subfolderPattern = argv[3]
    else:
        subfolderPattern = "job-(\d*).zip"
    if len(argv)>=5:
        timeInterval = int(argv[4])
    else:
        timeInterval = 60
    if len(argv)>=6:
        databaseFilename = argv[5]
    else:
        databaseFilename = "collected.db"
    if len(argv)>=7:
        databaseFilename_particles = argv[6]
    else:
        databaseFilename_particles = "particles.db"
except:
    print("Usage: autoZippedResultsCombiner.py parent_folder expected_number_of_zip_files [subfolder_pattern] [watch_time_interval (seconds)] [database_filename]")
    exit()

# watch loop
matchFilename = re.compile(subfolderPattern)
keepWatching = True
while keepWatching:
    count = 0
    if path.exists(parentFolder):
        fileList = listdir(parentFolder)
        for aFile in fileList:
            if matchFilename.match(aFile): count += 1
        print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        print("Number of zipped job files found: %d" % count)
    # keep watching?...
    if count < numberOfZipFiles:
        print("Keep watching...")
        stdout.flush()
        sleep(timeInterval)
        continue
    # enough!
    if count > numberOfZipFiles:
        print("Warning: found %d zipped job files, more than the expected number: %d" % (count, numberOfZipFiles))
    keepWatching = False

stdout.flush()
call("python ./combineEbeDatabasesFromZippedResults.py %s %s" % (parentFolder, databaseFilename), shell=True)
call("python ./combineEbeDatabasesFromZippedResults.py %s %s" % (parentFolder, databaseFilename_particles), shell=True)
