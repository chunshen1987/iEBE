#! /usr/bin/env python2

from sys import argv
from os import unlink, rename
import binISS
import UrqmdOutputFormatter

tmpFile = "./results/tmp.dat"

if len(argv) < 2:
    print("Usage: shell.py urqmd_output_file")
    exit(-1)

# format urqmd output file; so far only produce results for charged particle
UrqmdOutputFormatter.formatUrqmdOutputFile(argv[1], tmpFile)

# binning
binISS.use_bin_processes["calculate integrated flow"] = True
binISS.use_bin_processes["calculate differential flow"] = True
binISS.use_bin_processes["count particles in pT range"] = False
binISS.binISSDataFile(tmpFile, "./formats/binUrqmd_default_format.dat") # note that the format needs to agree with those dumped by UrqmdOutputFormatter.formatUrqmdOutputFile

rename("results/integrated_flow.dat", "results/integrated_flow_Charged.dat")
rename("results/differential_flow.dat", "results/differential_flow_Charged.dat")

# delete temp file
unlink(tmpFile)
