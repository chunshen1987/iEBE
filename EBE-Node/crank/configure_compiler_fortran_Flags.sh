#! /usr/bin/env bash

FC=`which ifort`;
FFlAGS=" -O3 -fast -heap-arrays -cpp"
if [ "$FC" == "" ]; then
   FC=`which gfortran`;
   FFLAGS=" -O3 -cpp"
fi
echo $FFLAGS
