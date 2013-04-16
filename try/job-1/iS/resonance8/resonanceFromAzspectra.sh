#! /usr/bin/env bash

# Perform resonance calculation for all subfolders under the target folder, which is the first argument. Subfolders containing $already_exist will be skipped. The subfolder I refer to is the subfolder that has the structure from the "spectraFromHydro.sh" script from Azspectra; meaning that all the dN/d3p data are stored under "target_folder/subfolder/spectra/".

# Example:
# ./resonanceFromAzspectra.sh ~/Downloads/test/ [executable]

operation_folder="spectra"
already_exist="spectra/spec_211.dat"

current_dir=`pwd`

# check argument and get directories
if [ $# -lt 1 ]
then
  echo "Usage: resonanceFromAzspectra.sh target_folder"
  exit
else
  targetDir=$1
fi

# get executable and parameters
executable="./resonance.e"
while [ $# -gt 1 ]
do
  arg=`echo "$2" | sed 's/ /\?/g'`
  executable="$executable $arg"
  shift
done

for sub_folder in `ls $targetDir`
do
  if [ -e $targetDir/$sub_folder/$already_exist ]
  then
    echo Skipping $sub_folder
  else
    echo Dealing with $sub_folder ...
    mv $targetDir/$sub_folder/$operation_folder/* ./results
    $executable
    mv ./results/* $targetDir/$sub_folder/$operation_folder/
  fi
done
