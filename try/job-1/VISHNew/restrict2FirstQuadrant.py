#!/usr/bin/env python

from sys import argv
from os import path
from os import path
from os import getcwd

# my library
from fileR import readData
from fileR import writeData

surf_from = "surface.dat"
dec_from = "decdat2.dat"

surf_to = "surface-1st.dat"
dec_to = "decdat-1st.dat"

def restrictToFirstQuadrant(file_path=getcwd()):
	""" Generate surface and decoupling data files in the first quadrant
	based on the data files on the whole transverse plane. """

	surf_from_full = path.join(file_path, surf_from);
	dec_from_full = path.join(file_path, dec_from);
	print("Read surface data file from "+surf_from_full);
	print("Read decoupling data file from "+dec_from_full);

	if not path.exists(surf_from_full) or not path.exists(dec_from_full):
		print("Surface or decoupling data files not found at "+file_path);
		return
	else:
		surf_to_full = path.join(file_path, surf_to);
		dec_to_full = path.join(file_path, dec_to);
		print("Write surface data file to "+surf_to_full);
		print("Write decoupling data file to "+dec_to_full);
		surf_buffer_read = readData(surf_from_full);
		dec_buffer_read = readData(dec_from_full);
		surf_buffer_write = [];
		dec_buffer_write = [];
		if len(2*surf_buffer_read) != len(dec_buffer_read):
			print("Number of data unmatch:"+"surface data file has "+str(len(surf_buffer_read))+"entries; "+"decoupling data file has "+str(len(dec_buffer_read))+"entries.");
			return;
		else:
			record = 0;
			for i in range(len(surf_buffer_read)):
				if surf_buffer_read[i][2]>=0 and surf_buffer_read[i][3]>=0:
					surf_buffer_write.append(surf_buffer_read[i]);
					dec_buffer_write.append(dec_buffer_read[2*i]);
					dec_buffer_write.append(dec_buffer_read[2*i+1]);
					record = record + 1
			writeData(surf_to_full, surf_buffer_write);
			writeData(dec_to_full, dec_buffer_write);
			print("Number of entries written: "+str(record));

if __name__ == "__main__":
	print("Usage: restrict2FirstQuadrant path_to_data_files")
	if (len(argv)>=2):
		restrictToFirstQuadrant(argv[1]);
	else:
		restrictToFirstQuadrant();
