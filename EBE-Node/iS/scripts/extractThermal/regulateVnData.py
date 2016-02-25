#!/usr/bin/env python

# This file uses the new format:
# order numerator_real numerator_imag flow_real flow_imag flow_magnitude

from os import path
from sys import argv

from dirR import listSubDirectories
from fileR import readData
from fileR import writeData

number_of_moments = 9; # how many moments are calculated

# The format table for reading v_n data (within each piece), it has the form [[moment_index (starting from 1), line_to_read (starting from 0)], ...]
format_vn = map(lambda ii: [ii+1,ii+1], range(number_of_moments+1)); #+2 means to skip the line for particle species, and the line for the zero-th moment
data_indices_vn = [1,2,3,4,5]; # indices of data to be added from each line

def regulateVnData(dir_path, piece=2, reg_path="thermal_211_integrated_vn.dat", vn_path="v2data-inte.dat"):
    """
        Regulate the file specified by vn_path by extracting only the n-th particle data
        and write it to a separated file specified by reg_path. The variable n is specified
        by "piece".
    """

    piece = int(piece);
    data = map(lambda x:[], range(number_of_moments+1)); # build an empty list for results. data[i] holds results for moment i (thus index 0 is not used)

    if (not path.exists(path.join(dir_path, vn_path))):
        print("Missing data files in directory "+dir_path);
        return;

    # declare data buffer
    tmp_data = [];

    # read and add v_n file:
    vn_data = readData(path.join(dir_path, vn_path), "!");

    # the data in the n-th piece starts from line (format_vn[-1][1]+1)*(n-1)
    vn_data_piece = vn_data[(format_vn[-1][1]+1)*(piece-1):(format_vn[-1][1]+1)*piece];

    order = 0;
    for idx in format_vn: # loop through indices
        tmp_line = map(lambda ii: vn_data_piece[idx[1]][ii], data_indices_vn); # v_n data
        tmp_line.insert(0, order);
        # tmp_line.append(vn_data_piece[format_vn[0][1]-1][1]); # add number of particles for this species
        tmp_data.append(tmp_line); # add desired v_n data
        order += 1;
    # write to file
    reg_full = reg_path;
    writeData(path.join(dir_path, reg_full), tmp_data);

    print("Finished.")

if __name__ == "__main__":
    if len(argv)==1 or len(argv)>4:
        print("Usage: regulateVnData.py directory_path [piece_index] [generated_datafile_relative_path] [v_n_datafile_relative_path] ")
    elif len(argv)==2:
        regulateVnData(argv[1]);
    elif len(argv)==3:
        regulateVnData(argv[1],argv[2]);
    elif len(argv)==4:
        regulateVnData(argv[1],argv[2],argv[3]);
    elif len(argv)==5:
        regulateVnData(argv[1],argv[2],argv[3],arg[4]);
