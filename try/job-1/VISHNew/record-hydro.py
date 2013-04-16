#!/usr/bin/env python

# system-wide
from sys import argv
from os import path
from os import getcwd
from time import strftime

# my libraries
from runR import run
from fileR import makeDir

def record_hydro(argv):
  """ Run the commands specified by argv[2:] then copy the generated data files to a uniquely generated directory under argv[1]. """

  properEndingString1 = "Early" # if this string is not found, the program must stoped in the middle
  properEndingString2 = "Decoupling"

  if len(argv) > 1:
    # generate necessary strings
    unique_signature = "-".join(argv[2:])+"--"+strftime("%Y-%m-%d@%H:%M:%S")
    # for debug: print unique_signature
    store_rel_path = unique_signature
    logfile = unique_signature+".log"
    exe_line = " ".join(argv[2:])+" | tee "+logfile

    # clean up first
    run("bash ./del.sh")

    # execute hydro code (presumably)
    print "Executing " + exe_line + "..."
    run(exe_line)

    # copy data files
    tmpfile = open(path.join(getcwd(), logfile))
    tmpBuffer =  tmpfile.read()
    if ((properEndingString1 in tmpBuffer) or (properEndingString2 in tmpBuffer)): # program normally exited
      print("Program quit normally, copy data now.");
      # create target directory
      store_path = path.join(argv[1],store_rel_path);
      store_path = makeDir(store_path, "new"); # get the next available dir
      run("bash ./copy-data-hydro.sh"+" "+str(getcwd())+" "+store_path)
    else:
      print("Error in executation.")
    # clean up again
    run("bash ./del.sh")
  else:
    print "Usage: record base-dir-path command (with parameters)"

if __name__ == "__main__":
  record_hydro(argv)
