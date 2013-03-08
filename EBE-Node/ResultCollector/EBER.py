#!/usr/bin/env python
# Version 0.5

"""
    This module consists of functions dealing with collection and processing of
    data files.

    Major types of files include: 11P5N, 2NP1Ecc.

    ***** 11P5N *****:
    The 11P5N file has 11+5*N columns. The first 4 columns are:
        real(s_n), imag(s_n), real(s'_n), imag(s'_n)
    which are the eccentricities defined using entropy profile, and 5-8 columns
    gives:
        real(e_n), imag(e_n), real(e'_n), imag(e'_n)
    which are the eccentricities defined using energy profile. Here the
    quantities with prime (') are defined using r^n weight, otherwise with r^2
    weight. Columns 9-10 are the real and imaginary part of v_n, and column 11
    is the particle number N. Starting from column 12, the data are grouped
    every 5 columns to store p_T differential quantities. The total number of
    columns depending on the actual number of p_T points used in calculations.
    For each group, the 5 columns give:
        p_T, m_T, dN/(2pi p_T dp_T), real(v_n(p_T)), imag(v_n(p_T)).
    The order n and particle species info are indicated by the file names.

    ***** 2NP1Ecc (r_power) *****:
    The 2NP1Ecc file has 2*N+1 columns. The first 2*N columns records the real
    and imagniary parts of the eccentricity <r^m cos(n*phi)> / <r^m> where
    n=1..N and m is indicated by the filename. The last column records
    int(r^m*weight) on the transverse plane. The row index is the n index. The
    type of profile used for weighting (sd or ed) and order n are indiciated by
    the file name. Here int(r^n) are the integral of r^n*profile. It can be used
    to calculate <r^n> etc.

    ***** 6Ecc *****:
    The 6Ecc file has 6 columns. The 1-3 columns are the real part, imaginary
    part, and magnitude of the eccentricity calculated using r^2 weight and the
    4-6 columns are similar but for r^n weight. The row index is the order index
    that starts at 1. The choice of weighting function is indicated by the file
    name.

    ***** 4Ecc *****:
    The 4Ecc file has 4 columns. The 1-3 columns are the real part, imaginary
    part, and magnitude of the eccentricity calculated using r^m weight, and m
    is indicated by the file name. The weight function is also indicated by the
    file name. The 4th column is int(r^m*weight). For m=0 the file records
    int(weight).

    ***** 6InteFlow *****:
    The 6InteFlow file has 6 columns. The 1st column gives the flow order index
    starting from 0. The 2nd and 3rd columns are the real and imaginary part of
    the numerator of the fraction that defines the integrated anisotropy flows.
    The 4-6 columns are the real part, imaginary part, and the magnitude of the
    integrated flow. For which particle are such flows calculated for and
    weather it is for thermal particles or after resonance decay are all
    indicated by the file name.

    ***** 3P3NDiffFlow *****
    The 3P3NDiffFlow file has 3+3*N columns. The 1-3 columns are for p_T, m_T,
    and dN/(2*pi p_T dp_T). The rest of the columns are grouped every 3 columns.
    For each group, the 1st, 2nd, and 3rd columns are for the real part,
    imaginary part, and the magnitude of the differential flow of order n, where
    n=1..N. Different rows correspond to different p_T values, and the total
    number of p_T points are determined by the actual integration method used in
    the program.

"""

from os import path
from sys import argv, exit

import numpy as np

from dirR import _toList, listDirectSubDirectoriesHavingFiles, listDirectSubDirectoriesHavingFilesSpanningMostIndices

amount_of_feedback = 99 # number of printed lines on the terminal

# declare default parameters
# all the particles that can be collected
default_particle_names=[
    "thermal_211", "thermal_321", "thermal_2212",
    'pion_p', 'Kaon_p', 'proton', 'Lambda', 'Sigma_p', 'Sigma_m', 'Omega', 'Xi_m',
    'Charged', 'Charged_ptcut02', 'Charged_ptcut03', 'Charged_ptcut05', # for dN/dy
    'R_0.5_Charged_eta', 'R_2.5_Charged_eta', 'R_2.0_Charged_eta', 'R_0.5_2.5_Charged_eta', 'R_0.5_4.8_Charged_eta',
    'R_0.5_Charged_ptcut02_eta', 'R_2.5_Charged_ptcut02_eta', 'R_2.0_Charged_ptcut02_eta', 'R_0.5_2.5_Charged_ptcut02_eta', 'R_0.5_4.8_Charged_ptcut02_eta',
    'R_0.5_Charged_ptcut05_eta', 'R_2.5_Charged_ptcut05_eta', 'R_2.0_Charged_ptcut05_eta', 'R_0.5_2.5_Charged_ptcut05_eta', 'R_0.5_4.8_Charged_ptcut05_eta',
    'R_0.5_Charged_ptcut03_eta', 'R_2.5_Charged_ptcut03_eta', 'R_2.0_Charged_ptcut03_eta', 'R_0.5_2.5_Charged_ptcut03_eta', 'R_0.5_4.8_Charged_ptcut03_eta',
    'R_0.5_Charged_ptcut03_3_eta', 'R_2.5_Charged_ptcut03_3_eta', 'R_2.0_Charged_ptcut03_3_eta', 'R_0.5_2.5_Charged_ptcut03_3_eta', 'R_0.5_4.8_Charged_ptcut03_3_eta',
    ]

# what name to save as (11P5N)
default_save_as_abr="11P5N_moments_order_%%d_%s.dat"

# flow and differential flow files (6InteFlow, 3P3NDiffFlow)
# use new format: both thermal and RD particles share the same pattern. Use extractThermal3 script from iS with version >= 1.2.1.9 to ensure this.
default_vn_file_path_abr="spectra/%s_integrated_vndata.dat"
default_diff_vn_file_path_abr="spectra/%s_vndata.dat"
# for backward compatibility
default_historical_thermal_vn_file_path_abr="spectra/%s_integrated_vn.dat"
default_historical_thermal_diff_vn_file_path_abr="spectra/%s_vn.dat"

# "conventional" eccentricity files (6Ecc)
default_sn_file_path="results/ecc-init-sd.dat"
default_en_file_path="results/ecc-init.dat"

# r-power eccentricity files (3Ecc, 3Ecc, 2NP1Ecc, 2NP1Ecc)
# these exotic eccentriciy files are only generated by VISHNew with version>=1.10.0
default_rp_sn_file_path='results/ecc-init-sd-r_power-%d.dat' # for r-power exotic ecc files from ebe folder; using sd
default_rp_en_file_path='results/ecc-init-r_power-%d.dat' # using ed
default_save_rp_sn_as='sd_ecc_r_power_%d.dat' # using what names to save exotic sd ecc files
default_save_rp_en_as='ed_ecc_r_power_%d.dat' # using what names to save exotic ed ecc files
default_r_power_start_with=0 # for backward compatibility use 1


#------------------------------------------------------------------------------
def particleNameIsThermal(particle_name):
    """
        Return True if the given particle_name is for a thermal particle.
    """
    if "thermal" in particle_name.lower():
        return True
    else:
        return False


#------------------------------------------------------------------------------
def loadSeriesOfFiles(filename_pattern, from_index=None, to_index=None):
    """ Load a series of files with filenames given by the pattern
    string filename_pattern, which contains a "%d" string specifying the
    order.
    -- filename_pattern: this string will be used with a index to
        generate the actual filename to be loaded. It should contain one
        "%d" to be filled by an acctual index.
    -- from_index, to_index: files with index between (include) two
        these two indices will be loaded. They will be auto-dectected by
        a simple checking-for-existence method starting from index 0 and
        increment 1.
    -- return: a dictionary whose keys are the indices and values are
        the loaded data. It also contains special keys "from_index",
        "to_index", and "filename_pattern" for control info.
    """
    # check indices boundary
    if from_index==None: from_index=0
    if to_index==None: # auto detect
        to_index=from_index
        while path.exists(filename_pattern % to_index):
            to_index += 1
        else:
            to_index -= 1

    # triviality check
    if to_index < from_index: return {}

    # load data
    data_collection = {"from_index": from_index, "to_index": to_index, "filename_pattern":filename_pattern};
    for ii in range(from_index, to_index+1): data_collection[ii] = np.loadtxt(filename_pattern % ii)
    return data_collection



#------------------------------------------------------------------------------
# All the data collection dictionary have two keys: The "header"
# identify the type of data stored in the dictionary and the "data"
# indentify the actual data read using the loadSeriesOfFiles function.
# The possible values for headers are listed below. Dictionary with
# key="header" that matches one of the following given values are
# assumed to store the corresponding data using the key="data".

# header are used to identify the dictionary that contains a type of data to all orders
collection_11P5N_header = "11P5NCollection" # 11P5N
collection_2NP1Ecc_header = "2NP1Ecc" # 2NP1Ecc

# To read 11P5N files to all orders
load11P5NFilesAllOrder = lambda foldername="/home/qiu/Downloads/Pb_w_resonance_moms-08-08-2012/Pb_0_5_Glb", filename_pattern="11P5N_moments_order_%d_Charged.dat", from_index=1, to_index=None: {"header": collection_11P5N_header, "data": loadSeriesOfFiles(path.join(foldername, filename_pattern), from_index, to_index)}

# To read eccentricity files to all r^m powers
load2NP1EccFilesAllOrder = lambda foldername="/home/qiu/Downloads/Pb_w_resonance_moms-08-08-2012/Pb_0_5_Glb", filename_pattern="ed_ecc_r_power_%d.dat", from_index=0, to_index=None: {"header": collection_2NP1Ecc_header, "data": loadSeriesOfFiles(path.join(foldername, filename_pattern), from_index, to_index)}




#------------------------------------------------------------------------------
def getCplxMomFromSingle11P5NBlock(block, moment_name, pT_idx=0):
    """ Return a complex moment from the 11P5N data block.
    -- block: the 11P5N data block.
    -- moment_name: it can be one of the following:
        -- "sn": the eccentricity calculated using entropy density and
            r^2 weight (0,1 columns).
        -- "spn": the eccentricity calculated using entropy density and
            r^n weight (2,3 columns).
        -- "en": the eccentricity calculated using energy density and
            r^2 weight (4,5 columns).
        -- "epn": the eccentricity calculated using energy density and
            r^n weight (6,7 columns).
        -- "vn": the flow (8,9 columns).
        -- "N": total particle yield (10 column).
        -- "diffVn": the differential flow v_n(p_T). Only when choosing
            this string will the pT_idx parameter be used.
    -- pT_idx: for which p_T index the info will be extracted. It starts
        from 1.
    -- return: The returned value is a np.matrix object containing
        the required info for all the events. When moment_name is not
        "diffVn", the returned matrix has only 1 column containing the
        complex values of the moment (or real in the case of "N"); when
        moment_name is "diffVn", the required value is a list of 4
        1-column matrices:
            p_T, m_T, dN/(2*pi p_Tdp_T), v_n(p_T)
    """
    if moment_name=="sn":
        return block[:,0] + 1j*block[:,1]
    elif moment_name=="spn":
        return block[:,2] + 1j*block[:,3]
    elif moment_name=="en":
        return block[:,4] + 1j*block[:,5]
    elif moment_name=="epn":
        return block[:,6] + 1j*block[:,7]
    elif moment_name=="vn":
        return block[:,8] + 1j*block[:,9]
    elif moment_name=="N":
        return block[:,10]
    elif moment_name=="diffVn":
        if pT_idx<1:
            print("getCplxMomFromSingle11P5NBlock error: pT_idx cannot be less than 1.")
            exit(-1)
        column_shift = 11 + (pT_idx-1)*5
        return [ block[:,column_shift], block[:,column_shift+1], block[:,column_shift+2], block[:,column_shift+3]+1j*block[:,column_shift+4] ]


#------------------------------------------------------------------------------
def taking_conjugate(mom, order):
    """
        Return mom if order>=0, otherwise np.conj(mom).
    """
    if order >= 0:
        return mom
    else:
        return np.conj(mom)


#------------------------------------------------------------------------------
def getCplxMomFrom11P5NCollection(collection, moment_name, order=2, pT_idx=1):
    """
        Return the complex moment from collection of 11P5N data blocks.
        The block is assumed to be read using the load11P5NFilesAllOrder
        function.

        -- collection: the returned dict by load11P5NFilesAllOrder.

        -- moment_name: one of the following:
            "sn", "spn", "en", "epn", "vn", "N", "diffVn"
             See getCplxMomFromSingle11P5NBlock.

        -- order: the n value of, e.g. v_n. Can be a list, in which case the
            returned value is an np.array whose 1st index is the order index
            and 2nd index is the event index. If it is an empty list, then all
            the specified moments of all order will be returned. When
            setting to negative values, the conjugate of of the moment
            is returned.

        -- pT_idx: used for differential data, see
            getCplxMomFromSingle11P5NBlock.

        -- return: the returned value is the same as by the function
            getCplxMomFromSingle11P5NBlock. This function only addes the
            possibility to choose order.
    """
    order = _toList(order)
    try:
        # check header
        if collection["header"] != collection_11P5N_header: raise KeyError
        # check special grammar of "order"
        if len(order)==0: order = range(collection["data"]["from_index"], collection["data"]["to_index"]+1)
        # check index
        if min(np.abs(order))<collection["data"]["from_index"] or max(np.abs(order))>collection["data"]["to_index"]:
            print("getCplxMomFrom11P5NCollection error: The requested order (index) does not lie inside the index bound for the read data.")
            exit(-1)
        return np.transpose(np.array([ taking_conjugate(getCplxMomFromSingle11P5NBlock(collection["data"][abs(anOrder)], moment_name, pT_idx), anOrder) for anOrder in order]))

    except TypeError and KeyError:
        print("getCplxMomFrom11P5NCollection error: The passed-in data collection has the wrong structure or empty or wrong identification key value.")
        exit(-1)


#------------------------------------------------------------------------------
# To return the integrated v_n all at once:
getAllVnFrom11P5NCollection = lambda collection: getCplxMomFrom11P5NCollection(collection, 'vn', order=[])


#------------------------------------------------------------------------------
def getCplxEccFrom2NP1EccCollection(collection, r_power=2, angle_order=2):
    """
        Return the complex eccentricity from collection of 2NP1Ecc data
        blocks. The block is assumed to be read using the
        load2NP1EccFilesAllOrder function.
        -- collection: the returned dict by load2NP1EccFilesAllOrder
            function.

        -- r_power: returned r^r_power weighted eccentricity

        -- angle_order: this is the order of anisotropy. If this value
            is negative, then the conjugate is returned.

        -- return: the return value is a 1 column matrix that stores the
            eccentricity <r^m cos(n phi)> where m=r_power and n=angle_order.
    """
    try:
        # check header
        if collection["header"] != collection_2NP1Ecc_header: raise KeyError
        # check r power
        if r_power<collection["data"]["from_index"] or r_power>collection["data"]["to_index"]:
            print("getCplxEccFrom2NP1EccCollection error: The requested r power (index) does not lie inside the index bound for the read data.")
            exit(-1)
        if angle_order==0:
            # return <r^m>
            if collection["data"]["from_index"] > 0:
                print("getCplxEccFrom2NP1EccCollection error: The calculation of <r^m> requires the present of int(r^0> data (i.e. r_power=0), which is missing.")
                exit(-1)
            else:
                return collection["data"][r_power][:,-1] / collection["data"][0][:,-1] # assume eqaul length and 1-1 correspondence
        else:
            # return eccentricity
            number_of_columns = collection["data"][r_power].shape[1]
            if abs(angle_order) > (number_of_columns-1)/2:
                print("getCplxEccFrom2NP1EccCollection error: The requested angular order is out-of-bound.")
                exit(-1)
            else:
                return taking_conjugate(collection["data"][r_power][:,(abs(angle_order)-1)*2] + 1j*collection["data"][r_power][:,(abs(angle_order)-1)*2+1], angle_order)

    except TypeError and KeyError:
        print("getCplxEccFrom2NP1EccCollection error: The passed-in data collection has the wrong structure or wrong identification key value.")
        exit(-1)


#------------------------------------------------------------------------------
def getAllEccFrom2NP1EccCollection(collection, type="rn_weighted"):
    """
        Return the specified eccentricity of all orders at once.

        This function calls the getCplxEccFrom2NP1EccCollection function.

        -- collection: The returned dict by load2NP1EccFilesAllOrder
            function.
        -- type: It can be one of the following.
            "rn_weighted": Return r^3-weighted e_1 and r^n-weighted e_n.
            "r2_weighted": Return r^2-weighted e_n.
        -- return: It is a np.array whose 1st index is the order index and
            2nd index is the event index.
    """
    try:
        # check header
        if collection["header"] != collection_2NP1Ecc_header: raise KeyError

        min_r_power = collection["data"]["from_index"]
        max_order = (collection["data"][min_r_power].shape[1] - 1) / 2

        if type=="r2_weighted":
            return np.array([ getCplxEccFrom2NP1EccCollection(collection, 2, order) for order in range(1,max_order+1) ])
        elif type=="rn_weighted":
            to_return = np.array([ getCplxEccFrom2NP1EccCollection(collection, order, order) for order in range(1,max_order+1) ])
            to_return[0, :] = getCplxEccFrom2NP1EccCollection(collection, 3, 1)
            return np.transpose(to_return)

    except TypeError and KeyError:
        print("getAllEccFrom2NP1EccCollection error: The passed-in data collection has the wrong structure or empty.")
        exit(-1)


#------------------------------------------------------------------------------
def collect11P5NColMomDataFromFolder(folder_path,
        particle_names=["pion_p"],
        save_as_abr=default_save_as_abr,
        vn_file_path_abr=default_vn_file_path_abr,
        diff_vn_file_path_abr=default_diff_vn_file_path_abr,
        sn_file_path=default_sn_file_path,
        en_file_path=default_en_file_path,
        up_to_order=None,
        echo_mode=amount_of_feedback):
    """
        Collect all moments from all subfolder of a given folder. The
        collected moments are:

        en(Sd), en'(Sd), en(Ed), en'(Ed), v_n (inte & diff)

        The first 4 quantities are collected from 2 6ColEcc files. The
        integrated v_n's are collected from a 6InteFlow file and the
        differential v_n's are collected from a 3P3NDiffFlow file. The
        collected data from all subfolders are written into a 11P5N file.

        All direct subfolders will be scanned, but only those containing all
        the required files will be collected. Required files are those
        specified by sn_file_path, en_file_path, and those by vn_file_path_abr
        and diff_vn_file_path_abr with their "%s" string replaced by those
        strings in particle_names.

        -- folder_path: Subfolders of this folder which contain all the
            required files will be collected.

        -- paticle_names: A list of strings, each of which will be used to
            substitue the "%s" in save_as_abr, vn_file_path_abr and
            diff_vn_file_path_abr to form complete file names.

        -- save_as_abr: The collected file will be saved to files with this
            name. It has a "%s" part that will be replaced by strings in the
            particle_names list; it also has a "%%d" part that will be
            replaced by the order of anisotropy. The "%s" holder will be
            evaluated first and the "%%d" holder will be evaluated second.

        -- vn_file_path_abr, diff_vn_path, sn_file_path, en_file_path: Names and
            relative path for the required data files. The v_n related files
            have "%s" part that will be filled with strings from the
            particle_names list.

        -- up_to_order: Up to which order will the anisotropies (ecc and flow)
            will be collected. If it is set to None, then it is determined
            from the file vn_file_path_abr % particle_names[0]. The starting
            value is always 1. (See 6InteFlow.)

        -- echo_mode: The larger this value is, the more output to the
            screen.
    """
    #
    # the beginning
    #
    if echo_mode>10: print("\n----- collect11P5NColMomDataFromFolder -----\n")

    particle_names = _toList(particle_names)
    number_of_particles = len(particle_names)

    # expand save_as_abr, vn_file_path_abr, and diff_vn_file_path_abr
    vn_file_path = [vn_file_path_abr % aStr for aStr in particle_names]
    diff_vn_file_path = [diff_vn_file_path_abr % aStr for aStr in particle_names]

    # get the list of sub-folders
    if echo_mode>10: print("Generating list for sub-directories...")
    # check for the following files
    files_to_check = [sn_file_path, en_file_path]
    files_to_check.extend(vn_file_path)
    files_to_check.extend(diff_vn_file_path)
    # now get the list of valid sub-folders with complete path
    if echo_mode>13:
        feedback_msg = "collect11P5NColMomDataFromFolder warning: subdirectory %s does not contain file %s; skipped..."
    else:
        feedback_msg = None
    valid_folders = map(lambda x: path.join(folder_path, x), listDirectSubDirectoriesHavingFiles(folder_path, files_to_check, feedback_msg)) # list of sub-directories with complete path
    # triviality check
    if len(valid_folders)==0:
        if echo_mode>10: print("No valid sub-directories to process.")
        return
    else:
        number_of_events = len(valid_folders)

    # take one particular sample to get the max anisotropy and number of p_T points info
    if up_to_order==None:
        vn_file = np.loadtxt(path.join(valid_folders[0], vn_file_path[0]))
        up_to_order = vn_file.shape[0]-1 # number of rows - 1 = max anisotropy order (assume the eccentricities use the same)

    diff_vn_file = np.loadtxt(path.join(valid_folders[0], diff_vn_file_path[0]))
    number_of_pTs = diff_vn_file.shape[0] # number of rows = number of p_T's

    # echo parameters
    if echo_mode>15:
        print("Total number of events: %d" % number_of_events)
        print("Use maximum anisotropy order: %d" % up_to_order)
        print("Use number of pT points: %d" % number_of_pTs)

    # pre-allocating space
    momsTable = [np.zeros([up_to_order, number_of_events, 11+5*number_of_pTs]) for particle_id in range(number_of_particles)] # 11P5N blocks


    #
    # the big loop; start to collect events
    #
    if echo_mode>10: print("Start collecting data from events...")
    event_id = 0 # always point to the next to-be-processed event
    for event_id in xrange(number_of_events):
        current_path = valid_folders[event_id]
        if echo_mode>=12:
            print('\nCurent event index: %d' % event_id)
            print('\nProcessing:\n %s\n' % current_path)

        # first, sd ecctricities
        block = np.loadtxt(path.join(current_path, sn_file_path))
        for particle_id in xrange(number_of_particles):
            for order_id in xrange(up_to_order):
                momsTable[particle_id][order_id, event_id, 0:2] = block[order_id, 0:2] # r^2 weighted
                momsTable[particle_id][order_id, event_id, 2:4] = block[order_id, 3:5] # r^n weighted

        # next, ed eccntricities from en file
        block = np.loadtxt(path.join(current_path, en_file_path))
        for particle_id in xrange(number_of_particles):
            for order_id in xrange(up_to_order):
                momsTable[particle_id][order_id, event_id, 4:6] = block[order_id, 0:2] # r^2 weighted
                momsTable[particle_id][order_id, event_id, 6:8] = block[order_id, 3:5] # r^n weighted

        # semi-finally, v_n data and hadrons yield (of this species)
        for particle_id in xrange(number_of_particles):
            block = np.loadtxt(path.join(current_path, vn_file_path[particle_id]))
            for order_id in xrange(up_to_order):
                momsTable[particle_id][order_id, event_id, 8:10] = block[order_id+1, 3:5] # order_id+1: the 1st row is for order=0
                momsTable[particle_id][order_id, event_id, 10] = block[0, 1] # 1st row, 2nd element is the particle yield

        # finale, differential v_n data
        for particle_id in xrange(number_of_particles):
            block = np.loadtxt(path.join(current_path, diff_vn_file_path[particle_id]))
            for order_id in xrange(up_to_order):
                momsTable[particle_id][order_id, event_id, 11::5] = block[:,0] # p_T
                momsTable[particle_id][order_id, event_id, 12::5] = block[:,1] # m_T
                momsTable[particle_id][order_id, event_id, 13::5] = block[:,2] # dN
                momsTable[particle_id][order_id, event_id, 14::5] = block[:,3+order_id*3] # v_n, real part
                momsTable[particle_id][order_id, event_id, 15::5] = block[:,4+order_id*3] # v_n, imaginary part

    # save collected data files
    if echo_mode>10: print("Writing data files...")
    for particle_id in xrange(number_of_particles):
        for order_id in xrange(up_to_order):
            to_save_filename = path.join(folder_path, save_as_abr % particle_names[particle_id] % (order_id+1)) # for human, order starts with 1 instead of 0
            np.savetxt(to_save_filename, momsTable[particle_id][order_id, :,:], fmt="%15.8e")

    if echo_mode>10: print("Done")



#------------------------------------------------------------------------------
def collectExoticEccFromFolder(folder_path,
        rp_eccn_file_path=default_rp_en_file_path,
        save_rp_eccn_as=default_save_rp_en_as,
        r_power_start_with=default_r_power_start_with,
        echo_mode=amount_of_feedback):
    """
        Collect all the exotic eccentricities from subdirectories of a folder.
        The exotic eccentricities are read from rp_sn_file_path and
        rp_en_file_path (3Ecc files). The collected files are 2NP1Ecc files and
        they are saved to save_rp_sn_as and save_rp_en_as files.

        -- folder_path: Data from subdirectories of this folder will be
            collected.

        -- rp_eccn_file_path: Exotic eccentricity files to read.
            They should have a "%d" part to be filled by the r-power.

        -- save_rp_eccn_as: Exotic eccentricity files to save. They
            should have a "%d" part to be filled by the r-power.

        -- r_power_start_with: Starting from which r^m power are the
            eccentricities calculated.

        -- echo_mode: The larger this value is, the more output to the
            screen.
    """

    #
    # the beginning
    #
    if echo_mode>10: print("\n----- collectExoticEccFromFolder -----\n")

    # get the list of sub-folders
    if echo_mode>10: print("Generating list for sub-directories...")
    if echo_mode>13:
        feedback_msg = "collectExoticEccFromFolder warning: subdirectory %s does not contain only %d exotic eccentricity files; skipped..."
    else:
        feedback_msg = None
    valid_folders, up_to_r_power = listDirectSubDirectoriesHavingFilesSpanningMostIndices(folder_path, rp_eccn_file_path, start_from_index=r_power_start_with, error_msg=feedback_msg) # list of sub-directories with complete path
    # triviality check
    if len(valid_folders)==0 or up_to_r_power<r_power_start_with:
        if echo_mode>10: print("No valid sub-directories to process.")
        return
    else:
        number_of_events = len(valid_folders)
    # make complete path
    valid_folders = map(lambda x: path.join(folder_path, x), valid_folders)
    number_of_events = len(valid_folders)


    # take one particular sample to get the max order for angular anisotropy
    ecc_file = np.loadtxt(path.join(valid_folders[0], rp_eccn_file_path % r_power_start_with))
    up_to_order = ecc_file.shape[0] # number of rows = order of angular anisotropy

    # echo parameters
    if echo_mode>15:
        print("Total number of events: %d" % number_of_events)
        print("Use maximum anisotropy order: %d" % up_to_order)
        print("Use maximum r-power: %d" % up_to_r_power)

    # pre-allocating space
    eccsTable = np.zeros([up_to_r_power-r_power_start_with+1, number_of_events, 2*up_to_order+1]) # 3NP1 blocks

    #
    # the big loop; start to collect events
    #
    if echo_mode>10: print("Start collecting data from events...")
    event_id = 0 # always point to the next to-be-processed event
    for event_id in xrange(number_of_events):
        current_path = valid_folders[event_id]
        if echo_mode>=12:
            print('\nCurent event index: %d' % event_id)
            print('\nProcessing:\n %s\n' % current_path)

        # get eccentricities
        for rpower_id in xrange(up_to_r_power-r_power_start_with+1):
            rpower = rpower_id + r_power_start_with
            block = np.loadtxt(path.join(current_path, rp_eccn_file_path % rpower))
            eccsTable[rpower_id, event_id, 0:-1:2] = block[:, 0] # real part of eccentricity
            eccsTable[rpower_id, event_id, 1:-1:2] = block[:, 1] # imaginary part of eccentricity
            eccsTable[rpower_id, event_id, -1] = block[0, -1] # int(r^rpower*weight)

    # save collected data files
    if echo_mode>10: print("Writing data files...")
    for rpower_id in xrange(up_to_r_power-r_power_start_with+1):
        rpower = rpower_id + r_power_start_with
        to_save_filename = path.join(folder_path, save_rp_eccn_as % rpower)
        np.savetxt(to_save_filename, eccsTable[rpower_id, :,:], fmt="%15.8e")

    if echo_mode>10: print("Done.")
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
if __name__ == "__main__":
    print("----------------------------------------------------")
    print("Welcome! -- Zhi Qiu")
    print("----------------------------------------------------")
    if len(argv) == 1:
        print("Use one of the following:\n")
        print("See the source file.")
    else:
        print("Executing: "+argv[1]+"('"+"','".join(argv[2:])+"')\n")
        exec(argv[1]+"('"+"','".join(argv[2:])+"')")
