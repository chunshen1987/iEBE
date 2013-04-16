#! /usr/bin/env bash

FC=`which ifort`
if [ "$FC" == "" ]; then
   FC=`which gfortran`;
fi
echo $FC
