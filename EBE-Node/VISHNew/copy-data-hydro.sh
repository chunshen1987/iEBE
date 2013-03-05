#! /usr/bin/env bash

if [ $# -eq 0 ]
then
  echo "Usage: copy-date-hydro.sh source destination"
  exit
elif [ $# -ge 2 ]
then
  curDir=$1
  targetDir=$2
elif [ $# -ge 1 ]
then
  curDir=`pwd`
  targetDir=$1
fi

cp $curDir/*.log $targetDir/
cp $curDir/*.for $targetDir/
cp -R $curDir/movie $targetDir/
cp -R $curDir/Initial $targetDir/
cp -R $curDir/results $targetDir/
cp -R $curDir/*.inp $targetDir/
