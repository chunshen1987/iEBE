#!/usr/bin/env python

from os import path
from sys import argv

from dirR import listSubDirectories
from fileR import readData
from fileR import writeData

number_of_pts = 15; # how many p_T

# The format table for reading v_n data (within each piece), it has the form [[moment_index (starting from 1), line_to_read (starting from 0)], ...]
format_vn = map(lambda ii: [ii+1,ii+2], range(number_of_pts)); #+2 means to skip the line for particle species, and the line for the zero-th moment
data_indices_vn = [1]; data_indices_vn.extend(range(3,31)); # indices of data to be added from each line

lines_to_read = lambda idx: range((number_of_pts)*(idx-1), (number_of_pts)*(idx-1)+number_of_pts); # which line to read for given index
# cols_to_read = range(0,30);

def regulateDiffVnData(dir_path, piece=1, reg_path="spectra/thermal_211_vn.dat", diff_vn_path="spectra/v2data.dat"):
	"""
		Regulate the file specified by diff_vn_path by extracting only the n-th particle data
		and write it to a separated file specified by reg_path. The variable n is specified
		by "piece".
	"""

	piece = int(piece);
	data = range(number_of_pts); # resulting list place holder
	#data = map(lambda x:[], range(number_of_moments+1)); # build an empty list for results. data[i] holds results for moment i (thus index 0 is not used)

	for sub_dir in listSubDirectories(dir_path):
		if (not path.exists(path.join(sub_dir, diff_vn_path))):
			print("Missing data files in subdirectory "+sub_dir);
			continue;

		# declare data buffer
		tmp_data = [];

		# read and add v_n file:
		diff_vn_data = readData(path.join(sub_dir, diff_vn_path));

		# get the n-th piece of data
		diff_vn_data_piece = map(lambda ii: diff_vn_data[ii], lines_to_read(piece));

		for aLine in diff_vn_data_piece: # loop through indices
			tmp_data.append(aLine[:]); # add desired v_n data
			# tmp_data.append(map(lambda ii: aLine[ii], cols_to_read)); # add desired v_n data

		# write to file
		reg_full = reg_path;
		writeData(path.join(sub_dir, reg_full), tmp_data);

	print("Finished.")

if __name__ == "__main__":
	if len(argv)==1 or len(argv)>4:
		print("Usage: regulateDiffVnData.py directory_path [piece_index] [generated_datafile_relative_path] [diff_v_n_datafile_relative_path]")
	elif len(argv)==2:
		regulateDiffVnData(argv[1]);
	elif len(argv)==3:
		regulateDiffVnData(argv[1],argv[2]);
	elif len(argv)==4:
		regulateDiffVnData(argv[1],argv[2],argv[3]);
	elif len(argv)==5:
		regulateDiffVnData(argv[1],argv[2],argv[3],arg[4]);
