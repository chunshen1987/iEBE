#!/usr/bin/env python

from sys import argv

if len(argv) < 3:
    print("Usage: freeStreamingShell.py sample_file sample_format_file")
    exit(-1)

import assignmentFormat
import freeStreaming

sample_filename = argv[1]
sample_format_filename = argv[2]
#sample_filename = "/home/qiu/Downloads/fancy_3d_movie/hydro_movie/en-movie.dat"
#sample_format_filename = "formats/hydro_movie_format.dat"

streaming_filename = "results/cells.dat"
streaming_format_filename = "results/hydro_streaming_format.dat"

from parameters_freeStreaming import time_interval, starting_time, hydro_endding_time, ed_dec, ed_smallness, hydro_extended_time

sample_format_dict = assignmentFormat.assignmentExprStream2IndexDict(file(sample_format_filename))

streaming_format_dict = freeStreaming.streamingHydroZ(
        file(sample_filename), sample_format_dict,
        file(streaming_filename, "w"),
        starting_time=starting_time,
        endding_time=hydro_endding_time,
        extended_time=hydro_extended_time,
        time_interval=time_interval,
        ed_dec=ed_dec, ed_smallness=ed_smallness,
        )

streaming_format_file = file(streaming_format_filename,"w")
for aTerm in assignmentFormat.dict2AssignmentExprList(streaming_format_dict):
    streaming_format_file.write(aTerm + "\n")
streaming_format_file.close()
