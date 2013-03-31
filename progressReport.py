#! /usr/bin/env python
"""
    Search and report progress of all current running jobs by looking at the the
    RunRecord.txt file inside direct subdirectories under the given path.
"""

from sys import argv
from os import path, listdir

# get optional arguments
if len(argv)>=2:
    targetWorkingDirectory = argv[1]
else:
    targetWorkingDirectory = "PlayGround"

from subprocess import call

for aFolder in listdir(targetWorkingDirectory):
    subFolder = path.join(targetWorkingDirectory, aFolder, "crank")
    recordFile = path.join(subFolder, "RunRecord.txt")
    if not path.exists(recordFile): continue # looking at the wrong subdirectory
    print("For %s:" % aFolder)
    commandString = 'grep "events out of" RunRecord.txt | tail -n 1'
    call(commandString, shell=True, cwd=subFolder)
