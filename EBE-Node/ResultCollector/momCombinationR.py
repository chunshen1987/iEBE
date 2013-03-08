#!/usr/bin/env python
# Version 0.1

"""
    This module consists of functions about forming moment combinations among
    flows or exotic eccentricities. The required data files are in 11P5N or
    2NP1Ecc formats, formed by functions from EBER. Here a moment is either a
    flow v_n or an eccentricity e_n, or anything complex.

    A combination is a product of several moments, for example:
        C235 = v_2 * v_3 / v_5
    The average of real part of the normalized combination is the "angle
    correlation" studied by several groups. For example:
        < real(C235) > = < cos( 2*Psi_2 + 3*Psi_3 - 5*Psi_5 ) >

    A combination is internally represented by the powers of the moments,
    including the skipped ones. For example, the combination
        C235 = v_1^0 * v_2^1 * v_3^1 * v_4^0 * v_5^(-1)
    is represented as [0, 1, 1, 0, -1].

    A (combination) relation list is a list of internally stored combinations.
    For example, a valid relation list can be:
        [ [0,  2,  0, -1], [0,  4,  0, -2] ]

    ***** momCombination_relation *****
    The relation file records the combinations used when the momCombination file
    is formed. It should be placed in the same folder as the momCombination
    file. It is a 1-line csv file, with each element a string indicating a
    combination. For example, if the first element is "[0 2 0 -1]", it gives the
    internal combination [0, 2, 0, -1].

"""

import numpy
import csv
from os import path
from sys import argv

from EBER import load11P5NFilesAllOrder, load2NP1EccFilesAllOrder, getAllVnFrom11P5NCollection, getAllEccFrom2NP1EccCollection
from fileRVer2 import writeCplxData

# default name for combination file
default_momCombination_filename = "flow_combinations.dat"
# for more sophisticated use
default_save_sn_combination_as = "sd_ecc_combinations.dat"
default_save_en_combination_as = "ed_ecc_combinations.dat"
default_save_flow_combination_as = "flow_combination_%s.dat"

# default name for combination relation file; the file contains info for how combinations are formed
default_momCombination_relation_filename = "moms_combination_relation.csv"

# default combination relation list; should only be used when no momCombination_relation file is present
default_momCombination_relation = [
    # fig 6
    [0,  2,  0, -1], # v2^2/v4
    [0,  4,  0, -2], # v2^4/v4^2
    [0,  6,  0, -3], # v2^6/v4^3
    [0,  3, -2], # v2^3/v3^2
    [0,  3,  0,  0,  0, -1], # v2^3/v6
    [0,  0,  2,  0,  0, -1], # v3^2/v6
    [0,  0,  4, -3], # v3^4/v4^3
    [0,  5,  0,  0, -2], # v2^5/v5^2
    # fig 7
    [0,  1,  1,  0, -1], # v2*v3/v5
    [0, -4,  1,  0,  1], # v2^(-4)*v3*v5
    [0,  1,  0,  1,  0, -1], # v2*v4/v6
    [0, -5,  0,  1,  0,  1], # v2^(-5)*v4*v6
    [0,  1, -2,  1], # v2*v3^(-2)*v4
    [0, -5,  2,  1], # v2^(-5)*v3^2*v4
]

amount_of_feedback = 99 # number of printed lines on the terminal


#------------------------------------------------------------------------------
def loadCombinationRelationFile(filename):
    """
        Load a momCombination_relation file and return the relation list it
        represents.
        -- filename: The complete path of the file to be read.
        -- return: The relation list from the file.
    """
    return [ numpy.matrix(element).tolist()[0] for element in [row for row in csv.reader(open(filename))][0] ]


#------------------------------------------------------------------------------
def saveCombinationRelationFile(filename, relation_list):
    """
        Save a relation list to a momCombination_relation file.
        -- filename: The complete path of the file to be save to.
        -- relation_list: The relation list to be saved.
    """
    csv.writer(file(filename, 'w')).writerow([str(numpy.array(element)) for element in relation_list])


#------------------------------------------------------------------------------
def formCombinationFromRelationList(relation_list, moms_cplx):
    """
        Combine moments' into combinations using relation list relation_list.

        -- relation_list: The relation list. Each element of relation_list
            corresponds to one relation between moments; it corresponds to a
            product of all moments, with each moment raised to the power given
            by the corresponding element, starting from moment_1. For example:
            [0 1 1 0 -1 0 ] -> v1^0*v2^1*v3^1*v4^0*v5^(-1)*v6^0 = v2*v3/v5
            Here "vn" means the complex vn vector: v_n*exp(i*n*Psi_n).

        -- moms_cplx: A numpy.array or numpy.matrix of all the complex moments.
            Its 2nd index is the order index. Thus moms_cplx[:,1] is the vector
            for all the complex moment_2's, etc.

        -- return: It has the same number of rows (1st index) as relation_cell,
            i.e. the number of relations. It has the same number of columns (2nd
            index) as moms_cplx, i.e. the number of events. Its (i,j) element is
            the combination calculated using i-th relation and from the j-th
            event.
    """

    # initialization
    number_of_relations = len(relation_list)
    number_of_events = moms_cplx.shape[0]
    results = numpy.zeros([number_of_relations, number_of_events]) + 1j*numpy.zeros([number_of_relations, number_of_events])

    # the big loop
    for relation_id in xrange(number_of_relations):
        power = numpy.array(relation_list[relation_id]) # a particular relation
        moms_cplx_restricted = moms_cplx[:, :len(power)] # get relevant mom's
        #power = diag(relation_this)*ones(number_of_used_vns, number_of_events); % the power matrix for all vn's for all events
        moms_cplx_powered = moms_cplx_restricted**power[None, :] # take power
        moms_cplx_product = moms_cplx_powered.prod(1) # multiplying power of vn's
        results[relation_id, :] = moms_cplx_product # record results

    return results


#------------------------------------------------------------------------------
def produceVnCombinationFrom11P5NColMomFiles(folderName, momsFilename,
    combination_save_as = default_momCombination_filename,
    relation_save_as = default_momCombination_relation_filename,
    momCombination_relation = default_momCombination_relation,
    echo_mode = amount_of_feedback):
    """
        Generate files containing combination for flows using the 11P5N moms
        files. It is assumed that the mom files are under folderName, and they
        are loaded using EBER.load11P5NFilesAllOrder function. The combination
        are formed using the relation momCombination_relation.

        -- folderName: The folder containing the 11P5N files.

        -- momsFilename: 11P5N file that contains the required data. Its
            filename should contain a "%d" symbol allowing the program to choose
            different orders.

        -- combination_save_as: Filename used to save the resulting combination
            files.

        -- relation_save_as: Filename used to save the corresponding combination
            relation file.

        -- momCombination_relation: This relation list will be used to form
            combinations.

        -- echo_mode: The larger this value is, the more output to the screen.
    """
    if echo_mode>10: print("\n----- produceVnCombinationFrom11P5NColMomFiles -----\n")

    # read in vn
    moms_all = load11P5NFilesAllOrder(folderName, momsFilename)
    vn_cplx_all = getAllVnFrom11P5NCollection(moms_all)

    # process relation
    if (echo_mode>10): print("Forming combinations...")
    momCombination = formCombinationFromRelationList(momCombination_relation, vn_cplx_all);

    # write to file
    if (echo_mode>10): print("Writing to files...")
    writeCplxData(path.join(folderName, combination_save_as), momCombination.transpose())

    # save combination relation file too
    saveCombinationRelationFile(path.join(folderName, relation_save_as), momCombination_relation)

    if (echo_mode>10): print("Done.")


#------------------------------------------------------------------------------
def produceEccCombinationFromExoticEccFiles(folderName, momsFilename,
    combination_save_as = default_momCombination_filename,
    relation_save_as = default_momCombination_relation_filename,
    momCombination_relation = default_momCombination_relation,
    ecc_type = "rn_weighted",
    echo_mode = amount_of_feedback):
    """
        Generate files containing combination for eccentricities using the
        2NP1Ecc exotic eccentricity files. It is assumed that the mom files are
        under folderName, and they are loaded using
        EBER.load2NP1EccFilesAllOrder function. The combination are formed using
        the relation momCombination_relation.

        -- folderName: The folder containing the exotic eccentricity files.

        -- momsFilename: 2NP1Ecc file that contains the required data. Its
            filename should contain a "%d" symbol allowing the program to choose
            different r powers.

        -- combination_save_as: Filename used to save the resulting combination
            files.

        -- relation_save_as: Filename used to save the corresponding combination
            relation file.

        -- momCombination_relation: This relation list will be used to form
            combinations.

        -- echo_mode: The larger this value is, the more output to the screen.
    """
    if echo_mode>10: print("\n----- produceEccCombinationFromExoticEccFiles -----\n")

    # read in ecc
    moms_all = load2NP1EccFilesAllOrder(folderName, momsFilename)
    ecc_cplx_all = getAllEccFrom2NP1EccCollection(moms_all, ecc_type)

    # process relation
    if (echo_mode>10): print("Forming combinations...")
    momCombination = formCombinationFromRelationList(momCombination_relation, ecc_cplx_all);

    # write to file
    if (echo_mode>10): print("Writing to files...")
    writeCplxData(path.join(folderName, combination_save_as), momCombination.transpose())

    # save combination relation file too
    saveCombinationRelationFile(path.join(folderName, relation_save_as), momCombination_relation)

    if (echo_mode>10): print("Done.")


#------------------------------------------------------------------------------
def getStrFromRelationList(relation_list, mode='actual', momName='vn', hideESd=False):
    """
        Convert a 1d relation list giving combinations among moments to a list
        of strings. Each element of relation_list corresponds to one relation
        between moments; it corresponds to a product of all moments, with each
        moment raised to the power given by the corresponding element, starting
        from v1. For example, when moment is taken to be vn, the relation:
            [ 0 1 1 0 -1 0 ] -> v1^0*v2^1*v3^1*v4^0*v5^(-1)*v6^0 = v2*v3/v5
        gives corresponding string is "v2^(1)*v3^(1)*v5^(-1)".

        -- relation_list: The list containing all the relations.

        -- mode: Can be 'actual', 'cos_angle', 'cos_angle_latex'. Different
           string are used for different modes.

        -- momName: Can be 'vn', 'en', 'sn', or just 'v', 'e', 's'.

        -- hideESd: The '(ed)' or '(sd)' string will be hided if set to True.
    """

    # this is the list of relation strings that will be returned
    relation_string_list = [];

    if mode == 'actual':
        # used for actual book-keeping
        if momName == 'vn' or momName == 'v':
            atom = 'v%d^(%d)'
        elif momName == 'en' or momName == 'e':
            atom = 'e%d^(%d)'
        elif momName == 'sn' or momName == 's':
            atom = 's%d^(%d)'
        # build up the string using robust loop
        for aRelation in relation_list:
            string_this = '';
            for order_id in range(len(aRelation)):
                element = aRelation[order_id]
                if element==0: continue
                string_this = string_this + atom % (order_id+1,element)
            # add this relation
            relation_string_list.append(string_this)


    elif mode == 'cos_angle':
        # used for generating readable human strings
        for aRelation in relation_list:
            string_this = 'cos(';
            for order_id in range(len(aRelation)):
                element = aRelation[order_id]
                if element==0: continue
                if element > 0:
                    string_this = string_this + '+%dPsi%d' % (element*(order_id+1), order_id+1)
                else:
                    string_this = string_this + '%dPsi%d' % (element*(order_id+1), order_id+1)
            # add a distinguishing tail
            if momName == 'vn' or momName == 'v':
                string_this = string_this + ') (EP)'
            elif momName == 'en' or momName == 'e':
                string_this = string_this + ') (PP-ed)'
            elif momName == 'sn' or momName == 's':
                string_this = string_this + ') (PP-sd)'
            # add this relation
            relation_string_list.append(string_this)


    elif mode == 'cos_angle_latex':
        # used for generating latex string for human
        if momName == 'vn' or momName == 'v':
            atom = '%d\\Psi_{%d}^\\mathrm{EP}'
        elif momName == 'en' or momName == 'e' or momName == 'sn' or momName == 's':
            atom = '%d\\Psi_{%d}^\\mathrm{PP}'
        # loop over body
        for aRelation in relation_list:
            string_this = '$\cos(';
            being_first = True # the 1st atom does not carry the "+" sign
            for order_id in range(len(aRelation)):
                element = aRelation[order_id]
                if element==0: continue
                if element > 0:
                    if not being_first: string_this = string_this + '+'
                else: # element<0
                    string_this = string_this + '{-}'
                string_this = string_this + atom % (element*(order_id+1), order_id+1)
                if being_first: being_first=False
            # add a distinguishing tail
            if momName == 'vn' or momName == 'v':
                string_this = string_this + ')$'
            elif momName == 'en' or momName == 'e':
                if hideESd:
                    string_this = string_this + ')$'
                else:
                    string_this = string_this + ')$ (ed)'
            elif momName == 'sn' or momName == 's':
                if hideESd:
                    string_this = string_this + ')$'
                else:
                    string_this = string_this + ')$ (sd)'
            # add this relation
            relation_string_list.append(string_this)

    return relation_string_list




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
