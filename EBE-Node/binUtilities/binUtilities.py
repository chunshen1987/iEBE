#! /usr/bin/env python

"""
    Philosophy:
        A bin process (BinProcess) is uniquely given by 2 object:
        1) The one determining the choice of bin (BinObject),
        2) The one taks action (ActionObject).
        The average process is another class (binner) that is inernal in BinProcess.
        By combining different BinObject and ActionObject different BinProcesses
        can be produced.
        The actual binning function binDataStream read the file and then just
        blindly calls all the given BinProcesses.
        The "format" of a data file is defined as a {"varaible_name": column_index}
        dictionary, see formatter module for more.
"""

from os import path
import numpy

from fileRVer2 import writeData, writeCplxData
from assignmentFormat import dict2AssignmentExprList

use_separator = ' ' # data read in will be separated by this character
use_zero = 1e-15

# the following conventions are used for the format files
default_count_string = "count"
default_avg_string = "_mean"
default_std_string = "_std"
default_real_part_string = "_real"
default_imag_part_string = "_imag"

######################################################################################
class DataBinner():
    """ This class separate and store data into bins and return the bin average and counts.
    """
    # initialize average and count buffer
    def __init__(self, numberOfBins):
        self.reinitialize(numberOfBins)
    def reinitialize(self, numberOfBins):
        self.binAvgBuffer = [None]*numberOfBins # buffer for averages
        self.binStdBuffer = [None]*numberOfBins # buffer for standard deviations
        self.countBuffer = [0]*numberOfBins
        self.numberOfBins = numberOfBins

    # perform average and increase count
    def pushSample(self, bin_idx, sample):
        """
            Average the current sample with the recorded averages in
            self.binAvgBuffer then increase the count. The input sample must be
            a numpy.array object. The bin_idx is the index of the bin (starting
            with 0) the sample belongs to.
        """
        if bin_idx<0: return # a bin_idx<0 can be used to skip samples
        if self.binAvgBuffer[bin_idx] == None:
            self.binAvgBuffer[bin_idx] = numpy.copy(sample)
            self.binStdBuffer[bin_idx] = numpy.copy(sample*sample)
            self.countBuffer[bin_idx] = 1
        else:
            self.binAvgBuffer[bin_idx] += sample
            self.binStdBuffer[bin_idx] += sample*sample
            self.countBuffer[bin_idx] += 1

    # Return average and count
    def getAvgAndCount(self):
        """ The returned binAvgBuffer is a numpy.matrix object for easier post process.
        """
        for ii in range(self.numberOfBins): # look for empty entries
            if self.countBuffer[ii] != 0: break # found one
        else:
            return [ [0], 0 ]*self.numberOfBins # all are empty
        blank = self.binAvgBuffer[ii]*0 # empty entries will be replaced by this entry
        for ii in range(self.numberOfBins): # replace empty entries by blank
            if self.countBuffer[ii] == 0:
                self.binAvgBuffer[ii] = numpy.copy(blank)
            else:
                self.binAvgBuffer[ii] /= float(self.countBuffer[ii])
                self.binStdBuffer[ii] /= float(self.countBuffer[ii])
                self.binStdBuffer[ii] -= self.binAvgBuffer[ii]**2
                self.binStdBuffer[ii] = numpy.sqrt(self.binStdBuffer[ii])
        return [self.binAvgBuffer, self.binStdBuffer, self.countBuffer]



#####################################################################################
class BinObject(object):
    """ This class defines a bin. Its interface is given below.
    """
    def decideBin(self, sample, sampleFormat):
        """
            Return a bin index. This is the (only) required interface.
        """
        raise TypeError("BinObject.decideBin is not defined.")


####################################################################################
class SingleVarBin(BinObject):
    """ This class defines a bin using a single variable.
    """
    def __init__(self, bins, variableName):
        """ Copy and sort the bin boundary array. The actual number of bins is the
            number of elements in the boundary plus 1. variableName gives the
            one used for binning, the absolute index is read using sampleFormat,
            see "farmatter".
        """
        self.binBoundaries = bins
        self.binBoundaries.sort()
        self.binMin = min(self.binBoundaries)
        self.binMax = max(self.binBoundaries)
        self.binLength = len(self.binBoundaries)
        self.numberOfBins = len(self.binBoundaries)-1
        self.variableName = variableName
        # for optimation; see if the bins are uniformly separated
        if sum(abs(numpy.linspace(self.binMin, self.binMax, self.binLength)-self.binBoundaries))<use_zero:
            self.isBinUniform = True
            self.binStep = (self.binMax-self.binMin) / float(self.numberOfBins)
            def wrap_decideBin(self,sample,sampleFormat):
                elem = sample[sampleFormat[variableName]]
                if elem < self.binMin or elem > self.binMax:
                    return -1 # skip such samples
                else:
                    return int((elem-self.binMin)/self.binStep)
            self.decideBin_kernal = wrap_decideBin
        else:
            self.isBinUniform = False
            self.decideBin_kernal = lambda self,sample,sampleFormat: numpy.searchsorted(self.binBoundaries, sample[sampleFormat[variableName]])

    def decideBin(self, sample, sampleFormat):
        """ Inherit the interface.
        """
        return self.decideBin_kernal(self, sample, sampleFormat)


#####################################################################################
class SingleVarBinCheckingField(SingleVarBin):
    """
        This class extend the SingleVarBin class by first performing a check on
        a field, then only return a positive bin number when the passed-in data
        agrees with the registered one for the desired field.
    """
    def __init__(self, bins, variableName, desiredValue, fieldName):
        """
            Register fieldName and its desiredValue; call superclass's
            constructor.
        """
        super(SingleVarBinCheckingField, self).__init__(bins, variableName)
        self.fieldName = fieldName
        self.desiredValue = desiredValue

    def decideBin(self, sample, sampleFormat):
        """
            First check the value for the registered field then call the
            superclass' decideBin method. The check is:
            sample[sampleFormat[fieldName]] == desiredValue.
        """
        if sample[sampleFormat[self.fieldName]] == self.desiredValue:
            return self.decideBin_kernal(self, sample, sampleFormat)
        else:
            return -1


#####################################################################################
class BlockBin(BinObject):
    """ This class uses data block index as the bin index. It uses elements read
        froma blockSizeStream as the size of each block and it returns the block
        index as bin index.
        Do NOT use one object from this class for multiple processes since this
        class replies on number of calls to determine bin index; use multiple
        objects instead.
    """
    def __init__(self, blockSizeStream):
        """ Get and store elements in blockSizeStream. Its length is used as the
            bin size and the number of times the decideBin is called is used to
            determine which block (=bin index) the current sample is in.
        """
        self.blockSizes = [var for var in blockSizeStream]
        self.numberOfBins = len(self.blockSizes)
        self.blockCounter = 0 # which block the current sample is in
        self.numberOfRecentCalls = 0 # number of calls in the current block

    def decideBin(self, sample, sampleFormat):
        """ Return a bin index = block index. The convention is that the bin index
            starts with 0 (used by DataBinner). The sample and sampleFormat
            argument is still carried to mimic the same interface as BinObject.
        """
        self.numberOfRecentCalls += 1
        if self.numberOfRecentCalls <= self.blockSizes[self.blockCounter]:
            # still in the same block
            return self.blockCounter
        else:
            # in the next block
            self.blockCounter += 1
            self.numberOfRecentCalls = 1 # not 0
            return self.blockCounter



#####################################################################################
class ActionObject():
    """ This is the class that defines the interface for the "action" class.
        All "action" class should have the same interface.
    """
    def action(self, sample, sampleFormat):
        """
            This is a required interface. It should return the quantitiy that
            will be analyzed (mean + std).
        """
        raise TypeError("ActionObject.action is not defined.")
    def getDataFormatStrings(self):
        """
            This is also a required interface. It should return a list of
            strings describing the functionalities of each quantity (in order)
            returned by action().
        """
        raise TypeError("ActionObject.getDataFormatStrings is not defined.")



#####################################################################################
class SingleVarValue(ActionObject):
    """ It returns the value of the given variable with name variableName.
    """
    def __init__(self, variableName):
        self.variableName = variableName

    def action(self, sample, sampleFormat):
        """ Takes action on sample and the returned values will be used for binning.
            The default behavior is just to return the variable given by variableName.
        """
        return sample[sampleFormat[self.variableName]]




####################################################################################
class CountInRange(ActionObject):
    """ Its result 1 if the specified variable variableName lies in the designated interval, specified
        by self.lowerBound and self.upperBound; otherwise 0.
    """
    def __init__(self, variableName, lowerBound,  upperBound):
        self.variableName = variableName
        self.lowerBound = lowerBound
        self.upperBound = upperBound

    def action(self, sample, sampleFormat):
        """ Return 1 if variable given by variableName lies in the range given by
            self.lowerBound and self.upperBound; otherwise 0.
        """
        return int(self.lowerBound<=sample[sampleFormat[self.variableName]]<=self.upperBound)



####################################################################################
class BinProcess():
    """ This class wrap up a BinObject and a DataBinner together to
        form a complete object used for binning. Its interface is
        defined in the following. Each time the pushSample function
        is called, the bin the pushed sample is determined then quantities
        that need to be averaged will be calculated (using action), then
        results are averaged and stored using DataBinner class.
        Note that this class assumes that the BinObject obeys the
        interface defined above.
    """
    def __init__(self, BinObject, ActionObject):
        self.binObj = BinObject
        self.actionObj = ActionObject
        self.binner = DataBinner(self.binObj.numberOfBins)
        self.saveTo = None # avg and std will be saved to this file
        self.saveFormatTo = None # corresponding format file will be saved to this file, if set to non-None
        self.useCplx = False # set to True if data are complex

    def pushSample(self, sample, sampleFormat):
        """ First decide which bin the sample falls into. The element
            that is used in binning is the one located in column given by
            sampleFormat[self.binVar] (see formatter). Then use action
            function to generate results, then pass results to DataBinner.
        """
        bin_idx = self.binObj.decideBin(sample, sampleFormat)
        results = self.actionObj.action(sample, sampleFormat)
        self.binner.pushSample(bin_idx, results)

    def returnAvgAndCount(self):
        """ Return the results. Note that the returned averages are
            numpy.matrix converted iterated lists.
        """
        avgs, count = self.binner.getAvgAndCount()
        return [[var.tolist() for var in avgs], count]

    def saveAvgAndCount(self):
        """ Save the results to the file specified by self.saveTo.
            The count will be added as an additional column to the existing
            averages. It should be overloaded when required.
        """
        # deal with numerics first
        avgs, stds, count = self.binner.getAvgAndCount()
        avgs1d = []
        stds1d = []
        nonzero_count = []
        for idx in range(len(count)):
            if count[idx]!=0:
                # non empty bin
                avgs1d.append(avgs[idx].reshape([1,avgs[idx].size]).tolist()[0]) # make avgs 1d lists
                stds1d.append(stds[idx].reshape([1,stds[idx].size]).tolist()[0]) # make stds 1d lists
                nonzero_count.append(count[idx])
        to_write = [[] for var in range(len(nonzero_count))]
        for idx in range(len(nonzero_count)):
            to_write[idx].append(nonzero_count[idx])
            to_write[idx].extend(avgs1d[idx]) # insert averages
            to_write[idx].extend(stds1d[idx]) # insert standard deviations
        if not self.useCplx:
            writeData(self.saveTo, to_write)
        else:
            writeCplxData(self.saveTo, to_write)
        # now the format strings
        if self.saveFormatTo:
            format_file = open(self.saveFormatTo,"w")
            col_idx = 1
            format_file.write(default_count_string + " = %d\n" % col_idx)
            strings_to_be_appended = self.actionObj.getDataFormatStrings()
            for aString in strings_to_be_appended:
                col_idx += 1
                format_file.write(aString + default_avg_string + default_real_part_string + " = %d\n" % col_idx)
            for aString in strings_to_be_appended:
                col_idx += 1
                format_file.write(aString + default_std_string + default_real_part_string + " = %d\n" % col_idx)
            if self.useCplx:
                col_idx += 1
                format_file.write(default_count_string + "_imag = %d\n" % col_idx)
                for aString in strings_to_be_appended:
                    col_idx += 1
                    format_file.write(aString + default_avg_string + default_imag_part_string + " = %d\n" % col_idx)
                for aString in strings_to_be_appended:
                    col_idx += 1
                    format_file.write(aString + default_std_string + default_imag_part_string + " = %d\n" % col_idx)


###################################################################################
def binDataStream(dataStream, dataFormat, BinProcesses, level_of_output=5):
    """
        This function goes over the data stream and generate bin results.
        It has a BinProcess lists whose elements are BinProcess objects that
        will be called for each line of read data.
        The data are recorded using dataFormat (see formatter).
    """
    # process data
    line_index = 1;
    for aLine in dataStream:
        aLine = numpy.fromstring(aLine, sep=use_separator)
        for aProcess in BinProcesses: aProcess.pushSample(aLine, dataFormat)
        if level_of_output >= 3:
            line_index +=1
            if line_index % 100000 == 0: print("Number of lines processed: %d" % line_index)

    # save results
    for aProcess in BinProcesses: aProcess.saveAvgAndCount()



###################################################################################
def splitDataStream(dataStream, dataFormat, binObject, tmpDirectory="split"):
    """
        This function goes over the data stream and split it into separated files.
        The index returned by binObject will be used as the filename under the
        folder whose name is given by tmpDirectory.
    """
    tmp_filename = path.join(tmpDirectory, "%d.tmp")
    already_opened = {}
    for aLine in dataStream:
        dataLine = numpy.fromstring(aLine, sep=use_separator)
        bin_idx = binObject.decideBin(dataLine, dataFormat)
        if bin_idx in already_opened.keys():
            already_opened[bin_idx].write(aLine)
        else:
            already_opened[bin_idx] = file(tmp_filename % bin_idx, 'w')
            already_opened[bin_idx].write(aLine)

####################################################################################
# 03-28-2013:
# Bug fix: In BinProcesses.saveAvgAndCount, when using complex file format, a
# column corresponding to the imaginary part of "count " (=0) should also have a
# corresponding format label in the corresponding format file; previous this was
# overlooked.
