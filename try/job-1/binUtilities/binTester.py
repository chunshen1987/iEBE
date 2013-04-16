#! /usr/bin/env python

#import binUtilities

# result = [x for x in filterFuncs.str2listFilter(["1 2", "3, 4"])]

# [x for x in filterFuncs.str2listFilter(filterFuncs.divideDataFilter(["1 1","2 2","3 3"], [1,1]))]

#for aBlock in filterFuncs.divideDataFilter(["1 1","2 2","3 3"], [1,1]):
#    aNumBlock = [var for var in filterFuncs.strStream2BlockFilter(aBlock)]
#    print(aNumBlock)

#dataStream = ["1 2 3", "4 5 6"]
#formatDict = {"pT=1", "phi=2"}

#dataStream = file("/home/qiu/Downloads/fluctuation_study/results_Edec_0.1/samples_111_short.dat")
#controlStream = file("/home/qiu/Downloads/fluctuation_study/results_Edec_0.1/samples_control_111_short.dat")
#
#f = lambda x: sum([y[0]*y[1] for y in x])
#
#test_object = binUtilities.binUtilities_iSS()
#
#result = test_object.binDataStream(dataStream, controlStream)
#
##result = test_object.binData([f], dataStream)
#
#print(result)
#binObject = binUtilities.singleVarBin([2.5], "A")
#blockBinObject = binUtilities.blockBin([2, 3, 5, 3, 2])
#actionObject = binUtilities.singleVarValue("A")
#binProc = binUtilities.binProcess(binObject, actionObject)
#blockBinProc = binUtilities.binProcess(blockBinObject, actionObject)
#binProc.saveTo = "test.dat"
#blockBinProc.saveTo = "blockTest.dat"
#binUtilities.binDataStream(range(15),  {"A":0},  [binProc, blockBinProc])

from os import path
from dirR import expandPath
import binISS

target_folder = expandPath("~/Downloads/iSS_V2.1.1.0/results")

binISS.binISSDataFile(path.join(target_folder, "samples_211.dat"),
                      path.join(target_folder, "samples_format.dat"),
                      path.join(target_folder, "samples_control_211.dat"))
