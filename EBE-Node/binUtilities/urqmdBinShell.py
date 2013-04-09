#! /usr/bin/env python2

from sys import argv
from os import unlink, rename
import binISS
import UrqmdOutputFormatter

tmpFile = "./results/tmp.dat"

if len(argv) < 2:
    print("Usage: shell.py urqmd_output_file")
    exit(-1)

# format urqmd output file
UrqmdOutputFormatter.formatUrqmdOutputFile(argv[1], tmpFile)

# binning
useBinProcesses = []
useBinProcesses.append(binISS.differentialFlowCharged)
useBinProcesses.append(binISS.integratedFlowCharged)
useBinProcesses.extend(binISS.generateFlowActionsForPids(
        ( ("pion", 101), ("kaon", 106), ("nucleon", 1) )
    ))
binISS.binISSDataFileSimple(tmpFile, "./formats/binUrqmd_default_format.dat", useBinProcesses) # note that the format needs to agree with those dumped by UrqmdOutputFormatter.formatUrqmdOutputFile

# delete temp file
unlink(tmpFile)
