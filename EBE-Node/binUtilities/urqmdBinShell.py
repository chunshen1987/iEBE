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
if UrqmdOutputFormatter.formatUrqmdOutputFile(argv[1], tmpFile):

    # binning
    useBinProcesses = []
    useBinProcesses.append(binISS.differentialFlowCharged)
    useBinProcesses.append(binISS.integratedFlowCharged)
    useBinProcesses.extend(binISS.generateFlowActionsForPids(
            (
                ("pion", 101),
                ("kaon", 106), ("anti_kaon", -106),
                ("nucleon", 1), ("anti_nucleon", -1),
                ("sigma", 40), ("anti_sigam", -40),
                ("xi", 49), ("anti_xi", -49),
                ("lambda", 27), ("anti_lambda", -27),
                ("omega", 55), ("anti_omega", -55),
                ("phi", 109),
            )
        ))
    binISS.binISSDataFileSimple(tmpFile, "./formats/binUrqmd_default_format.dat", useBinProcesses) # note that the format needs to agree with those dumped by UrqmdOutputFormatter.formatUrqmdOutputFile

    # delete temp file
    unlink(tmpFile)
