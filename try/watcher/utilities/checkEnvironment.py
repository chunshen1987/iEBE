#! /usr/bin/env python
"""
    Print a list of tests to see whether all required tools for Ebe calculations
    are present.
"""

from os import getcwd, unlink
from subprocess import call

numberOfSpaces = 5

def printWarning(warningString):
    print("-"*(numberOfSpaces-2) + "> " + warningString)

def printMsg(message):
    print(" "*numberOfSpaces + message)

def checkCommand(cmdString, utilityName=None):
    """
        Try to execute "cmdString", then use "utilityName" to echo messages.
    """
    if not utilityName: utilityName=cmdString
    call("%s &> response.txt" % cmdString, shell=True, cwd=getcwd())
    if "command not found" in open("response.txt").readline():
        printWarning("%s *NOT* installed." % utilityName)
        unlink("response.txt")
        return False
    else:
        printMsg("%s installed." % utilityName)
        unlink("response.txt")
        return True

def checkModule(moduleName):
    """
        Try to import "moduleName", then echo messages.
    """
    try:
        __import__(moduleName)
        printMsg("python %s module installed." % moduleName)
        return True
    except:
        printWarning("python %s module *NOT* installed." % moduleName)
        return False

finalMsgs = []

print("Start checking...")
print("-"*80)

# check g++ and icpc
if not checkCommand("g++") and not checkCommand("icpc"):
    finalMsgs.append("You need to install icpc or g++.")

# check gfortran and ifort
if not checkCommand("gfortran") and not checkCommand("ifort"):
    finalMsgs.append("You need to install ifort or gfortran.")

# check make utility
if not checkCommand("make"):
    finalMsgs.append("You need to install the make utility.")

# check gsl
if not checkCommand("gsl-config", "gsl"):
    finalMsgs.append("You need to install gsl library.")

# check zip and unzip
if not checkCommand("zip --help", "zip") or not checkCommand("unzip --help", "unzip"):
    finalMsgs.append("You need both zip and unzip utilities.")

# check numpy
if not checkModule("numpy"):
    finalMsgs.append("You need to install python numpy package.")

# print final messages
print("-"*80)
if not finalMsgs:
    print("All essential packages installed. Test passed.")
else:
    for msg in finalMsgs: print(msg)
