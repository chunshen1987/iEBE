#!/usr/bin/env python

from sys import argv

if len(argv) < 5:
    print("Usage: freeStreaming2DAlongZShell.py input_file starting_z velocity_z evolution_time")
    exit(-1)

import assignmentFormat
import freeStreaming

from parameters_freeStreaming import time_interval

input_filename = argv[1]
starting_z = float(argv[2])
velocity_z = float(argv[3])
evolution_time = float(argv[4])

output_filename = "results/streaming_particle3d.dat"
output_dict_filenamt = "results/streaming_particle3d_format.dat"

streaming_format_dict = freeStreaming.freeStream2dAlongZ(
    file(input_filename), file(output_filename,"w"),
    starting_z, velocity_z, evolution_time,
    time_interval=time_interval)

streaming_format_file = file(output_dict_filenamt,"w")
for aTerm in assignmentFormat.dict2AssignmentExprList(streaming_format_dict):
    streaming_format_file.write(aTerm + "\n")
streaming_format_file.close()
