#!/usr/bin/env python
"""
    This is one of the shells to the EbeCollector class. This one
    creates a database using data from subfolders containing multiple
    hybrid (hydro+UrQMD) events.
"""

from sys import argv, exit
from os import path

try:
    from_folder = path.abspath(argv[1])
    multiplicity_factor = float(argv[2])
except:
    print("Usage: shell from_folder multiplicity_factor [sub_folder_pattern] [database_filename]")
    exit()

# get optional parameters
if len(argv)>=4:
    subfolder_pattern = argv[3]
else:
    subfolder_pattern = "event-(\d*)"
if len(argv)>=5:
    database_filename = argv[4]
else:
    database_filename = "collected.db"

# call EbeCollector
from EbeCollector import EbeCollector
EbeCollector().collectParticleinfo(from_folder, subfolder_pattern)
EbeCollector().createDatabaseFromEventFolders(from_folder, subfolder_pattern, database_filename, collectMode="fromUrQMD", multiplicityFactor=multiplicity_factor)
