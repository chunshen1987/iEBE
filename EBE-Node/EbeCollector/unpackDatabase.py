#!/usr/bin/env python
"""
    Use SqliteDB to read a database file then unpack it into separated data
    files using the unpackDatabase function.
"""

# first check argv
from sys import argv
if len(argv)<2:
    print("Usage: unpackDatabase dataFilename [unpackToFolder] [useHeaderWithCommentSymbol]")
    exit(-1)

# next call the function from SqliteDB
from os import path
from DBR import SqliteDB

# get write-to folder
if len(argv)>=3:
    unpackToFolder = path.abspath(argv[2])
else:
    unpackToFolder = "."

# determine whether to use header
if len(argv)>=4:
    useHeader = (True, argv[3])
else:
    useHeader = (True, "# ")

SqliteDB(argv[1]).unpackDatabase(sep="\t", writeToFolder=unpackToFolder, writeHeader=useHeader)

print("Thanks for using. Zhi Qiu 03/2013")

