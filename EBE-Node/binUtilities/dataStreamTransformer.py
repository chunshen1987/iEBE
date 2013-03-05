#! /usr/bin/env python

""" All functions in this module are intended to be "filters".
    They should be used to transform one type of data to another.
    All of they all use iterators for efficiency.
"""

from listR import isFloat, separateStr

# Version 1

def dataStreamToBlockStream(dataStream, separationList=[]):
    """
        Read a data stream and for each time it is called, yield a
        block. The blocks contain number of elements specified by
        the sperationList, which when left empty will causes this function
        to do nothing, reducing it to the identity function.
        -- dataStream: anything that supports the iterator interface.
        -- separationList: anything that supports iterator interface.
            Each element of this list is used to control how many number of
            lines should be read from dataStream. The element returned by
            this list will be converted to floats. If no such list is
            provided or it is exausted, then each line of the data stream
            is yielded in turn. In the former case it is possible that the
            stream is exausted before this control list is.
    """
    # initialization of iterator
    control_iter = iter(separationList)

    # get number of lines to be read, initialize data block to be yield
    try:
        number_of_lines_to_read = float(control_iter.next())
    except StopIteration:
        number_of_lines_to_read = 1
    data = []

    # now loop over the stream
    for aLine in dataStream:
        # check if data needs to be yielded
        if number_of_lines_to_read == 0:
            # first yield previous block
            yield data
            # need to know how many lines to read next, and initialize data block
            try:
                number_of_lines_to_read = float(control_iter.next())
            except StopIteration:
                number_of_lines_to_read = 1
            data = []
        # deal with the read line, append it to the block
        data.append(aLine)
        number_of_lines_to_read -= 1
    if data: yield data # only when non-empty


def strStream2BlockStream(dataStream, numericalOnly=False, separators=[",", " ", ";"], commentSymbol="#", skipEmpty=True):
    """
        Read a string data stream and for each time it is called, yield a list
        corresponding to data converted from the corresponding string.
        -- dataStream: anything that when using the iterator interface
            gives a string.
        -- numericalOnly: when set to true, only numerical data will be
            included in the return block.
        -- separators: used to break string into substrings.
        -- commentSymbol: anything after this symbol in the string read
            from the dataStream will be discarded.
        -- skipEmpty: whether empty lines are skipped (not yielded). A line
            is empty if the final result is an empty list. Note that a non-empty
            string may give an empty list because of the numericalOnly parameter.
    """
    # loop over the stream
    for aLine in dataStream:
        # deal with the read line, append it to the block
        commentSymbolLocation = aLine.find(commentSymbol)
        if commentSymbolLocation != -1: aLine = aLine[:commentSymbolLocation] # contains the comment symbol
        lineData = [] # store the converted results
        splited = separateStr(aLine, separators); # split the line
        for piece in splited:
            if isFloat(piece): # is numerical
                lineData.append(float(piece))
            else: # is a string
                if not numericalOnly: lineData.append(piece)
        if lineData: # check if it's empty
            yield lineData


def permuteDataStream(self, dataStream, actualFormatDict, desiredFormatDict, indexStartsWith1=True):
    """ Read dataStream then reorder and yield them the data.
        The read-in data have format given in actualFormatDict.
        The output data have format given in desiredFormatDict.
        See formatter module for the meaning of "format".
        If indexStartsWith1 is True then the column indices are
        assumed to start with 1.
    """
    # deal with default parameters
    if actualFormatDict=={}: actualFormatDict = self.expectedFormatDict

    # get the function needed for reordering
    reorderFunc = formatter.getReorderingFunction(actualFormatDict, desiredFormatDict, indexStartsWith1)

    # deal with data
    counter = 1;
    for aLine in dataStream:
        yield reorderFunc(aLine)
        if counter % 500 == 0: print("Blocks analyzed: "+str(counter))
        counter += 1


