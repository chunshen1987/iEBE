#!/usr/bin/env python
"""
    Read in a urqmd output file and output a space separated version.
"""

from sys import argv
from math import sqrt, atan2
import os


def formatUrqmdOutputFile(urqmdOutputFilePath, formattedFilePath):
    """
        Scan the output of urqmd and reformat it into space separated data.
        The current version outputs E, pT, phi, and particle id.
        Return False when failed.
    """

    # check input file
    if not os.access(urqmdOutputFilePath, os.F_OK):
        print("Cannot open urqmd file: "+urqmdOutputFilePath)
        return False

    # check output file
    formattedFileHandler = file(formattedFilePath, 'w')
    if not formattedFileHandler:
        print("Cannot create file "+formattedFileHandler+" for output.")
        return False

    # perform the conversion
    read_mode = "header_first_part"
    header_count = 1 # the first read line is already part of the header line
    data_row_count = 0
    for aLine in file(urqmdOutputFilePath):
        if read_mode=="header_first_part":
            if header_count <= 14: # skip first 14 lines
                header_count += 1
                continue
            # now at 15th line
            assert header_count==15, "No no no... Stop here."
            try:
                data_row_count = int(aLine.split()[0])
            except ValueError as e:
                print("The file "+urqmdOutputFilePath+" does not have a valid urqmd output file header!")
                print(e)
                return False
            read_mode = "header_second_part"
        elif read_mode=="header_second_part":
            # skip current line by switching to data reading mode
            read_mode = "data_part"
        elif read_mode=="data_part":
            if data_row_count>0:
                # still have data to read
                try:
                    p0, px, py, pz = map(lambda x: float(x.replace("D","E")), aLine[98:193].split())
                    pid = int(aLine[216:222])
                    pT = sqrt(px*px + py*py)
                    phi = atan2(py, px)
                    E = p0
                    formattedFileHandler.write("%d  %g  %g  %g\n" % (pid, E, pT, phi))
                except ValueError as e:
                    print("The file "+urqmdOutputFilePath+" does not have valid urqmd data!")
                    print(e)
                    return False
                data_row_count -= 1
            if data_row_count == 1: # note: not 0, but 1
                # switch back to header mode
                data_row_count = 0
                header_count = 0 # not pointing at the header line yet
                read_mode = "header_first_part"



if __name__ == '__main__':
    if len(argv)<3:
        print("Usage: UrqmdOutputFormatter urqmdFile formattedFile")
    else:
        formatUrqmdOutputFile(argv[1], argv[2])
        #print(argv)
