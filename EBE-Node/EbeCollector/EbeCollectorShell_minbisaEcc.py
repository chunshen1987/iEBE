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
except:
    print("Usage: %s from_folder [multiplicityFactor deformedNuclei(yes or no) database_filename]" % argv[0])
    exit()

# get optional parameters
if len(argv) >= 5:
    database_filename = argv[4]
else:
    database_filename = "minbiasEcc.db"
if len(argv) >= 4:
    deformedNuclei = argv[3]
else:
    deformedNuclei = "no"
if len(argv) >= 3:
    multiplicityFactor = float(argv[2])
else:
    multiplicityFactor = 1.0

if deformedNuclei == "no":
    deformedFlag = False
else:
    deformedFlag = True

# call EbeCollector
from EbeCollector import EbeCollector
EbeCollector().collectMinbiasEcc(from_folder, database_filename, multiplicityFactor, deformedFlag)
