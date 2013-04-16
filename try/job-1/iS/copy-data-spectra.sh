#! /usr/bin/env bash

if [ $# -eq 0 ]
then
  echo "Usage: copy-data-spectra.sh source destination"
elif [ $# -ge 2 ]
then
  sourceDir=$1
  targetDir=$2
elif [ $# -ge 1 ]
then
  sourceDir=`pwd`
  targetDir=$1
fi

mkdir $targetDir/spectra
cp $sourceDir/*.log $targetDir/spectra
cp -R $sourceDir/results/*vn* $targetDir/spectra/
cp -R $sourceDir/results/dN* $targetDir/spectra/
cp -R $sourceDir/results/spec* $targetDir/spectra/
cp -R $sourceDir/results/v2* $targetDir/spectra/
rm *.log
