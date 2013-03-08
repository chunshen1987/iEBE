###################    Last edited on Oct. 26, 2010   #################################
#                                                              version 23 --- Zhi Qiu
# level 2
"""
    Provide functions related to file operations or data file manipulations,
    and functions that are closely related to them.
"""

import dirR
import listR

from os import path, getcwd, remove, rename, mkdir, rmdir, listdir;
from sys import platform, exit
from shutil import copy;
from subprocess import Popen;

__q_debug__ = False


# --- Used in all functions that call external programs ---
_xterm_path_directory = {"darwin":"/usr/X11/bin/xterm", "linux2":"/usr/bin/xterm"}
if platform in _xterm_path_directory.keys():
    _default_xterm = _xterm_path_directory[platform]
else:
    _default_xterm = _xterm_path_directory["linux2"]

_konsole_path_directory = {"darwin":"/usr/X11/bin/konsole", "linux2":"/usr/bin/konsole"}
if platform in _konsole_path_directory.keys():
    _default_konsole = _konsole_path_directory[platform]
else:
    _default_konsole = _konsole_path_directory["linux2"]

_gnome_terminal_path_directory = {"darwin":"/usr/bin/gnome-terminal", "linux2":"/usr/bin/gnome-terminal"}
if platform in _gnome_terminal_path_directory.keys():
    _default_gnome_terminal = _gnome_terminal_path_directory[platform]
else:
    _default_gnome_terminal = _gnome_terminal_path_directory["linux2"]

_default_terminal = _default_xterm

def run(command, cwd=getcwd()):
    """ Invoke a command and wait for it to stop. """
    proc = Popen(command, shell=True, cwd=cwd)
    while proc.communicate() != (None,None):
        pass

def runInTerm(command, cwd=getcwd(), terminal=_default_terminal):
    """ Invoke a command and wait for it to stop. """
    proc = Popen(terminal+" -e "+command, shell=True, cwd=cwd)
    while proc.wait() != 0:
        pass

def runInTermUnlimited(command, cwd=getcwd(), sleepTime=0, terminal=_default_terminal):
    """ Invode a command and wait for it to stop.
    ulimit is set to unlimited before the execution of
    the program.
    """
    tmpsh = open(path.join(cwd,"QZTEMP.sh"), "w")
    tmpsh.write("ulimit -s unlimited\n")
    tmpsh.write(command+"\n")
    tmpsh.write("sleep "+str(sleepTime))
    if __q_debug__ == True:
        tmpsh.write("sleep 10")
    tmpsh.close()
    proc = Popen(terminal+" -e bash QZTEMP.sh", shell=True, cwd=cwd)
    while proc.wait() != 0:
        pass
    #remove(path.join(cwd,"QZTEMP.sh"))


def sortByColumn(data_file, column_to_sort=1):
    """
        Sort the data file by the specified column. The specified
        column must be numerical.
    """
    cmp = lambda qvar: float(qvar.split()[column_to_sort-1]) # used in "sorted" method
    in_file = open(data_file, "r")
    out_file = open(data_file+".TEMP", "w")
    for a_line in sorted(in_file.readlines(),key=cmp):
        out_file.write(a_line)
    in_file.close()
    out_file.close()
    copy(data_file+".TEMP", data_file)
    remove(data_file+".TEMP")


def switchColumn(data_file, column1, column2):
    """
        Switch two columns specified by column1 and column2 in data_file.
    """
    data = []
    for dataLine in readData(data_file):
        tmp = dataLine[column1-1]
        dataLine[column1-1] = dataLine[column2-1]
        dataLine[column2-1] = tmp
        data.append(dataLine)
    writeData(data_file, data)

def copyFile(filename, sourceDir, targetDir, renameTo=None, silent=True):
    """ Copy file from sourceDir to targetDir.
    """
    if renameTo == None: renameTo = filename
    fullname_source = path.join(sourceDir, filename)
    fullname_target = path.join(targetDir, renameTo)
    copy(fullname_source, fullname_target)
    if silent==False:
        print("File "+fullname_source+" copied to "+source_dir)

def copy(source, target):
    """ Copy source file to target. """
    copy(source, target)


def copyFiles(fileNames, sourceDir, targetDir, silent=True):
    """ Copy files whose names are given in the list filenames
    from sourceDir to targetDir.
    """
    fileNames = listR.toList(fileNames)
    for filename in fileNames:
        copyFile(filename, sourceDir, targetDir, None, silent)


def xcopy(namePatterns, sourceDir, targetDir, renameTo=None, flags=None):
    """
        Copy files that match namePatterns in sourceDir to targetDir.
        All files match namePatterns are copied if renameTo is set to None;
        otherwise the first one matches is copied and renamed to renameTo.
        flags are used in listFilesMatch function (indirectly in re module).
    """
    nameL = dirR.listFilesMatch(sourceDir, namePatterns, flags)
    if len(nameL) == 0: return
    if not path.exists(targetDir): makeDir(targetDir)
    if renameTo == None:
        for name in nameL:
            full_source_path = path.join(sourceDir, name)
            full_target_path = path.join(targetDir, name)
            copy(full_source_path, full_target_path)
    else:
        full_source_path = path.join(sourceDir, nameL[0])
        full_target_path = path.join(targetDir, renameTo)
        copy(full_source_path, full_target_path)



def nestedXcopy(namePatterns, sourceDir, targetDir, renameTo=None, flags=None):
    """
        Copy files that match namePatterns under sourceDir to the corresponding
        directories under targetDir (i.e. they have the same relative path).
        All files match namePatterns are copied if renameTo is set to None;
        otherwise the first one matches is copied and renamed to renameTo.
        flags are used in listFilesMatch function (indirectly in re module).
    """
    for aDir in dirR.listNestedDirContainsOneOfFilesM(sourceDir, namePatterns, flags):
        xcopy(namePatterns, aDir, path.join(targetDir, dirR._relativePathString(sourceDir, aDir)), renameTo, flags)



def readData(filename, commentSymbol="#"):
    """
        Read a data file and return a nested list (data block).
        Each line contains data from each row (sub-list).
        All lines containing the commentSymbol (the whole lines)
        are ignored.
    """
    inFile = open(filename, "r")
    data = []
    for aLine in inFile.readlines():
        if aLine.find(commentSymbol)!=-1: continue
        lineData = []
        splited = aLine.split(); # split the line
        if splited == []: continue # skip empty lines
        for piece in splited:
            if listR.isFloat(piece): # is numerical
                lineData.append(float(piece))
            else: # is a string
                lineData.append(piece)
        data.append(lineData)
    inFile.close()
    return data


def readDataI(filename, commentSymbol="#"):
    """
        Read a data file and return a nested list (data block).
        Each line contains data from each row (sub-list).
        All lines containing the commentSymbol (the whole lines)
        are ignored. This is the iterator version.
    """
    inFile = open(filename, "r")
    data = []
    for aLine in inFile.readlines():
        if aLine.find(commentSymbol)!=-1: continue
        lineData = []
        splited = aLine.split(); # split the line
        for piece in splited:
            if listR.isFloat(piece): # is numerical
                lineData.append(float(piece))
            else: # is a string
                lineData.append(piece)
        if lineData: yield lineData # non-empty ones only
    inFile.close()


def readNumericalData(filename, commentSymbol="#"):
    """
        Read a data file and return a nested list (data block).
        Each line contains data from each row (sub-list).
        All lines containing the commentSymbol are ignored.
        Only numerical data are read.
    """
    inFile = open(filename, "r")
    data = []
    for aLine in inFile.readlines():
        if aLine.find(commentSymbol)!=-1: continue
        lineData = []
        splited = aLine.split(); # split the line
        for piece in splited:
            if listR.isFloat(piece): # is numerical
                lineData.append(float(piece))
        if lineData!=[]: data.append(lineData); # add only non-empty lines
    inFile.close()
    return data


def readNumericalDataI(filename, commentSymbol="#"):
    """
        Read a data file and return a nested list (data block).
        Each line contains data from each row (sub-list).
        All lines containing the commentSymbol (the whole lines)
        are ignored. Only numerical data are read.
        This is the iterator version.
    """
    inFile = open(filename, "r")
    data = []
    for aLine in inFile.readlines():
        if aLine.find(commentSymbol)!=-1: continue
        lineData = []
        splited = aLine.split(); # split the line
        for piece in splited:
            if listR.isFloat(piece): # is numerical
                lineData.append(float(piece))
        if lineData: yield lineData; # yield only non-empty lines
    inFile.close()


def writeData(filename, data, seperator="   "):
    """
        Write a nested list (data block) into a file.
        Each line contains data from each row (sub-list).
    """
    try:
        iter(data[0])
    except TypeError:
        data = [[value] for value in data] # impose a sub list
    outFile = open(filename, "w")
    for dataLine in data:
        outFile.write(seperator.join(map(repr, dataLine))+"\n")
    outFile.close()


def readCplxData(filename, commentSymbol="#"):
    """
        Read a data file and return a nested list (data block).
        Each line contains data from each row (sub-list).
        All lines containing the commentSymbol are ignored.
        Only numerical data are read.
        This is the complex version. It is assumed that there
        are even number of columns and the first half of them
        store the real part of the variables and the second
        half store the imaginary part.
    """
    inFile = open(filename, "r")
    data = []
    for aLine in inFile.readlines():
        if aLine.find(commentSymbol)!=-1: continue
        lineData = []
        splited = aLine.split(); # split the line

        # usual read & convert
        for piece in splited:
            if listR.isFloat(piece): # is numerical
                lineData.append(float(piece))

        # skip empty ones
        if not lineData: continue

        # fold columns to form complex data
        number_of_columns = len(lineData)
        if number_of_columns % 2 !=0:
            print("fileRVer2.readCplxData Error: the number of columns in %s must be even!" % filename)
            exit()
        half_columns = int(number_of_columns/2)
        data.append([lineData[col]+1j*lineData[col+half_columns] for col in range(half_columns)]); # add only non-empty lines
    inFile.close()
    return data



def writeCplxData(filename, data, seperator="   "):
    """
        Write a nested list (data block) into a file.
        Each line contains data for each row (sub-list).
        The data can be complex and the store order is
        such that the real part are first stored then
        imaginary part in following columns. For example,
        if the data is [1,2+3j] then the store format
        is "1 2 0 3".
    """
    try:
        iter(data[0])
    except TypeError:
        data = [[value] for value in data] # impose a sub list (also impose column format)
    outFile = open(filename, "w")
    for dataLine in data:
        outFile.write(seperator.join([repr(var.real) for var in dataLine]+[repr(var.imag) for var in dataLine])+"\n")
    outFile.close()



def nestedRenameFiles(dir_path, old_filenames, new_filenames, silent_level=0):
    """
        Rename all files in "old_filenames" under "dir_path" to
        "new_filenames". The first file in "old_filenames" will be
        renamed to the first file in "new_filenames", and similar for
        the rest. If "leaf_only" is specified as "True" then only
        files in leaf subdirectories are modified. "silent_level"=0:
        no output on screen; 1: short output; 2: full output
    """
    old_filenames = listR.toList(old_filenames)
    new_filenames = listR.toList(new_filenames)
    for old_name, new_name in zip(old_filenames, new_filenames):
        dirL = dirR.listNestedDirContainsFiles(dir_path, old_name)
        for aDir in dirL:
            if path.exists(path.join(aDir, new_name)):
                print("File "+path.join(aDir, new_name)+" already exists! skipped.")
            else:
                rename(path.join(aDir, old_name), path.join(aDir, new_name))
                if silent_level==0:
                    pass
                elif silent_level==1:
                    print("File "+_relativePath(dir_path, full_path)
                          +" renamed to "+_relativePath(dir_path,new_full_path))
                else:
                    print("File "+full_path+" renamed to "+new_full_path)




def nestedRenameFilesAdd(dir_path, filenames, str_add, add_to_front=True):
    """
        Rename all files in "filenames" under "dir_path" by adding
        the string "str_add" to the front. If "add_to_front" is
        specified as "False", changes will be made to the end of the
        file name, otherwise (by default) it will be add to the front
        of the file name.
    """
    filenames = listR.toList(filenames)
    for name in filenames:
        if add_to_front == True:
            str_mode = str_add+"%s%s"
        else:
            str_mode = "%s" + str_add + "%s"
        new_filename = str_mode % (path.splitext(name))
        nestedRenameFiles(dir_path, name, new_filename)


def nestedDeleteFiles(dir_path, filenames, silence_level=0, leaf_only=False):
    """
        Delete all files in the list "filenames" under "dir_path".
    """
    filenames = listR.toList(filenames)
    for name in filenames:
        dirL = dirR.listNestedDirContainsFiles(dir_path, name)
        for aDir in dirL:
            remove(path.join(aDir, name))
            if silence_level>0: print("File "+full_path+" deleted.")


def groupingDataOneL(dir_path, data_filename, frontAdd=None, isValid=None, change_name_to=""): # low level function
    """
        Group data files by copying them into one large file with the
        same name placed in the subdirectory one level up. If
        "frontAdd" is given, a string suggested by "frontAdd" function
        using directory name as argument will be added to the file. If
        "isValid" is given, only for those directories that it
        returns True (the input of isValid is the full path of the dir), the
        data file will be combined, and the string "change_name_to"
        will be added to the tail of the combined data
        file name.

    """
    if isValid == None: # build a trivial "isValid" function
        isValid = lambda qvar: True
        cmb_data_filename = data_filename
    else: # modify file name for the combined data file
        if change_name_to == "": change_name_to = data_filename
        cmb_data_filename = change_name_to

    untreated = []
    for aHLDir in dirR.nested_oneL_oneSubDir_hasAll(dir_path, data_filename):
        if __q_debug__: print(aHLDir)
        toWrite = open(path.join(aHLDir, cmb_data_filename), "w") # the combined data file is created here
        is_empty = True # combined data file is empty <==> no treatment
        for aDir in listdir(aHLDir):
            full_path = path.normpath(path.join(aHLDir, aDir))
            if path.isdir(full_path) == False: continue # not a directory
            # if hasNoSubDir(full_path) == False: continue # not a bottom directory
            if path.exists(path.join(full_path, data_filename)) == False: continue # does not contain the data file
            if not isValid(full_path): continue # skip this directory according to "isValid"
            is_empty = False
            if frontAdd != None: # add a string (usually a number) suggested by the name of dir in the combined data file
                toWrite.write(frontAdd(aDir)+" ")
            toRead = open(path.join(full_path, data_filename), "r")
            textBuffer = toRead.read()
            if textBuffer[-1:] != "\n": textBuffer = textBuffer + "\n" # force one and only one return after each data file
            toWrite.write(textBuffer)
            toRead.close()
        toWrite.close()
        if is_empty == True:
            remove(path.join(aHLDir, cmb_data_filename))
            untreated.append(aHLDir)

    return untreated



def makeDir(dir_path, when_conflicts="skip"): # 01-23-2011
    """ Make directory at dir_path. If parent directory does not exist,
    it is created too. The parameter when_conflicts determines the
    action taken when dir_path already exists. If it is "skip" (default)
    or "overwrite", then no action is taken; if it is "new", then a new
    folder with suffix "-#" will be created. The difference between
    opetion "skip" and "overwrite" is that "skip" makes the return value
    an empty string when the directory already exists.
    """
    if path.exists(dir_path):
      if when_conflicts=="skip":
        return "";
      elif when_conflicts=="overwrite":
        return dir_path;
      elif when_conflicts=="new":
        ii = 1;
        while path.exists(dir_path+"-"+str(ii)): ii = ii + 1;
        mkdir(dir_path+"-"+str(ii));
        return dir_path+"-"+str(ii);
      else:
        print("makeDir error: unknown when_conflicts option.");
        return "";
    dir_path = path.realpath(dir_path)
    dir_path = path.normpath(dir_path)
    if path.exists(path.dirname(dir_path)):
        mkdir(dir_path);
        return dir_path;
    else:
        makeDir(path.dirname(dir_path));
        mkdir(dir_path);
        return dir_path;


def removeDir(dir_path):
    """ Remove a directory. """
    if not path.exists(dir_path): return
    for name in listdir(dir_path):
        full_path = path.join(dir_path, name)
        if path.isfile(full_path):
            remove(full_path)
            continue
        if path.isdir(full_path):
            removeDir(full_path)
            continue
    rmdir(dir_path)


def delete(file):
    """ Delete a file. """
    if not path.exists(file): return
    remove(file)


def extractToken(filename, token, numOfLines=2):
    """ Return a list of strings consist of numOfLines lines
    in filename that follow immediately after the first line
    that contains string token.
    """
    file = open(filename, "r")
    lines = file.readlines()
    file.close()
    for aLine in lines:
        if aLine.find(token) != -1:
            return lines[lines.index(aLine)+1:lines.index(aLine)+1+numOfLines]
    return ""


def addColumnsToFile(filename, columns, add_before_original=True):
    """ Add columns of data into a file with filename, before or after the original data
    in each line. The variable "columns" should be a list of strings. Each string will be
    inserted accordingly into the file. This list will be re-used if it is shorter than the
    length of the file.
    """
    tempFile = "TEMP.tmp"
    outFile = open(tempFile, "w")
    inFile = open(filename, "r")
    columns = listR.toList(columns)
    index = 0
    aLine = inFile.readline()
    while aLine:
        if add_before_original:
            outFile.write(columns[index] + aLine)
        else:
            outFile.write(aLine + columns[index])
        index = listR.next(columns, index)
        aLine = inFile.readline()
    inFile.close
    outFile.close()
    copy(tempFile, filename)



if __name__ == "__main__":
    print("Morning!")
