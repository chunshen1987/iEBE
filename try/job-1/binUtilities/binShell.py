#! /usr/bin/env python

from sys import argv

if len(argv) < 3:
    print("Usage: shell.py sample_file sample_format_file sample_block_division_file")
    exit(-1)

import binISS
from dirR import expandPath
#import cProfile

binISS.use_bin_processes["calculate integrated flow"] = False
binISS.use_bin_processes["count particles in pT range"] = False

#if len(argv)==4:
#    cProfile.run(binISS.binISSDataFile(expandPath(argv[1]), expandPath(argv[2]), expandPath(argv[3])))
#else:
#    cProfile.run('binISS.binISSDataFile(expandPath(argv[1]), expandPath(argv[2]))', expandPath('~/Downloads/report.prof'))

print("Executing: " + "binISS.binISSDataFile"+"('"+"','".join([expandPath(var) for var in argv[1:]])+"')")
exec("binISS.binISSDataFile"+"('"+"','".join([expandPath(var) for var in argv[1:]])+"')")