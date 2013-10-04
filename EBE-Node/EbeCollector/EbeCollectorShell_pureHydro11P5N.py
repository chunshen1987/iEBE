#!/usr/bin/env python
"""
    This is one of the shells to the EbeCollector class. This one
    creates a database using data from subfolders containing multiple
    pure hydro events from 11P5N.
"""

from sys import argv, exit
from os import path

try:
    from_folder = path.abspath(argv[1])
except:
    print("Usage: shell from_folder [database_filename]")
    exit()

# get optional parameters
if len(argv)>=3:
    database_filename = argv[2]
else:
    database_filename = "collected.db"

# call EbeCollector
from EbeCollector import EbeCollector
EbeCollector().createDatabaseFromEventFolders(from_folder, databaseFilename=database_filename, collectMode="fromPureHydro11P5N")
