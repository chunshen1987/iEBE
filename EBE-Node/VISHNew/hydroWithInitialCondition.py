#! /usr/bin/env python

from sys import argv, exit;
from os import path, getcwd;

from listR import union;
from fileR import copy, makeDir, ls;
from runR import run;


def hydroWithInitialCondition(argv):
  """ Run hydro using intitial conditions specified in a source folder.
  Usage: hydroWithInitialCondition.py direct_init_filename source_folder destination_folder hydro_executable and parameters
  """

  initial_dir = "Initial";

  # check argument and get directories
  if len(argv)<4:
    print("Usage: hydroWithInitialCondition.py direct_init_filename source_folder destination_folder hydro_executable and parameters");
    exit();
  else:
    direct_init_file=argv[1];
    sourceDir=argv[2];
    targetDir=argv[3];

  current_dir = getcwd();

  # get hydro_executable and parameters
  executable = " ".join(argv[4:]);

  # get those initial conditions that already exist
  makeDir(targetDir);
  files_exist = [];
  for aDir in ls(targetDir, "dir"):
    files_exist = union([files_exist, ls(path.join(targetDir, aDir, initial_dir), "file")]);

  # the big loop
  for init_file in ls(sourceDir, "file"):
    # first check if this initial has already be used
    if init_file in files_exist:
      print("\n********************************\n");
      print("Initial file"+" "+init_file+" "+"already exists; skipping it.");
      continue;
    run("rm ./Initial/*");
    print("\n=============================================================\n");
    print("Using"+" "+path.join(sourceDir, init_file)+"...");
    copy(path.join(sourceDir, init_file), path.join(current_dir, initial_dir, init_file));
    copy(path.join(sourceDir, init_file), path.join(current_dir, initial_dir, direct_init_file));
    print("python ./record-hydro.py" + " " + targetDir + " " + executable);
    run("python ./record-hydro.py" + " " + targetDir + " " + executable);


if __name__ == "__main__":
  hydroWithInitialCondition(argv);
