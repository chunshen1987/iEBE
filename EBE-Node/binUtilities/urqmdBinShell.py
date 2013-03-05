#! /usr/bin/env python

from sys import argv
from os import unlink
import binISS
import UrqmdOutputFormatter

tmpFile = "./results/tmp.dat"

if len(argv) < 2:
    print("Usage: shell.py urqmd_output_file")
    exit(-1)

# format urqmd output file
UrqmdOutputFormatter.formatUrqmdOutputFile(argv[1], tmpFile)

# binning
binISS.use_bin_processes["calculate integrated flow"] = True
binISS.use_bin_processes["calculate differential flow"] = False
binISS.use_bin_processes["count particles in pT range"] = False
binISS.binISSDataFile(tmpFile, "./formats/binUrqmd_default_format.dat") # note that the format needs to agree with those dumped by UrqmdOutputFormatter.formatUrqmdOutputFile

# delete temp file
unlink(tmpFile)
