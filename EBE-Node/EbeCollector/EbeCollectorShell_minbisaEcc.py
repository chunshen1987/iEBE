#!/usr/bin/env python
"""
    This is one of the shells to the EbeCollector class. This one
    creates a database using data from superMC output files for 
    initial condition statistics.
"""

from sys import argv, exit
from os import path

try:
    from_folder = path.abspath(argv[1])
    multiplicityFactor = float(argv[2])
except:
    print("Usage: %s from_folder multiplicityFactor [database_filename]" % argv[0])
    exit()

# get optional parameters
if len(argv)>=4:
    database_filename = argv[3]
else:
    database_filename = "minbiasEcc.db"

# call EbeCollector
from EbeCollector import EbeCollector
EbeCollector().collectMinbiasEcc(from_folder, database_filename, multiplicityFactor)
