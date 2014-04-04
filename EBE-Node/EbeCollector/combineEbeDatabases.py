#!/usr/bin/env python
"""
    This script looks for Ebe database files under direct subdirectories of a
    given directory and combined them into a large one in the given directory.
    This script uses the mergeDatabases function from EbeCollector module. Only
    databases with the same name will be merged together, and the resulting
    database has the same name, only different locations.
"""

from sys import argv, exit
from os import path, listdir

try:
    parentFolder = path.abspath(argv[1])
except:
    print("Usage: combineEbeDatabases.py parent_folder")
    exit()

from DBR import SqliteDB
from EbeCollector import EbeCollector
collector = EbeCollector()
# loop over subdirectories
for aSubfolder in listdir(parentFolder):
    subfolder = path.join(parentFolder, aSubfolder)
    if not path.isdir(subfolder): continue # not a directory
    for aFile in listdir(subfolder):
        if path.splitext(aFile)[1] == ".db":
            print("Merging %s from %s..." % (aFile, aSubfolder))
            if path.splitext(aFile)[0] == "particles":
                collector.mergeparticleDatabases(SqliteDB(path.join(parentFolder, aFile)), SqliteDB(path.join(subfolder, aFile))) # merge a database to a database in parent folder with the same name.
            elif path.splitext(aFile)[0] == "minbiasEcc":
                collector.mergeMinbiasDatabases(SqliteDB(path.join(parentFolder, aFile)), SqliteDB(path.join(subfolder, aFile))) # merge a database to a database in parent folder with the same name.
            else:
                collector.mergeDatabases(SqliteDB(path.join(parentFolder, aFile)), SqliteDB(path.join(subfolder, aFile))) # merge a database to a database in parent folder with the same name.

print("Done.")
