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
                #("pion_p", "101_2"), ("pion_m", "101_-2"), ("pion_0", "101_0"),
                #("kaon_p", "106_1"), ("kaon_0", "106_-1"),
                #("kaon_m", "-106_-1"), ("anti_kaon_0", "-106_1"),
                #("proton", "1_1"), ("neutron", "1_-1"),
                #("anti_proton", "-1_-1"), ("anti_neutron",-1_1"),
                #("sigma_p", "40_2"), ("sigma_m", "40_-2"), ("sigma_0", "40_0"),
                #("anti_sigma_p", "-40_-2"), ("anti_sigma_m", "-40_2"), ("anti_sigma_0", "-40_0"),
                #("xi_0", "49_1"), ("xi_m", "49_-1"),
                #("anti_xi_0", "-49_-1"), ("anti_xi_m", "-49_1"),
                #("lambda", "27_0"), ("anti_lambda", "-27_0"),
                #("omega", "55_0"), ("anti_omega", "-55_0"),
                #("phi", "109_0"),

                # new pid: isospin*2000 + pid
                ("pion_p", 2101), ("pion_m", -1899), ("pion_0", 101),
                ("kaon_p", 1106), ("kaon_0", -894),
                ("kaon_m", -1106), ("anti_kaon_0", 894),
                ("proton", 1001), ("neutron", -999),
                ("anti_proton", -1001), ("anti_neutron", 999),
                ("sigma_p", 2040), ("sigma_m", -1960), ("sigma_0", 40),
                ("anti_sigma_p", -2040), ("anti_sigma_m", 1960), ("anti_sigma_0", -40),
                ("xi_0", 1049), ("xi_m", -951),
                ("anti_xi_0", -1049), ("anti_xi_m", 951),
                ("lambda", 27), ("anti_lambda", -27),
                ("omega", 55), ("anti_omega", -55),
                ("phi", 109),
            )
        ))
    binISS.binISSDataFileSimple(tmpFile, "./formats/binUrqmd_default_format.dat", useBinProcesses) # note that the format needs to agree with those dumped by UrqmdOutputFormatter.formatUrqmdOutputFile

    # delete temp file
    unlink(tmpFile)
