#! /usr/bin/env bash

# C compiler
CC=`which icc`; CFLAGS=" -O3"
if [ "$CC" == "" ]; then
  CC=`which g++`; CFLAGS=" -O3"
fi
CFLAGS="$CFLAGS `gsl-config --cflags --libs`"

# FORTRAN compiler
FC=`which ifort`; FFLAGS=" -O3 -fast -heap-arrays -cpp"
if [ "$FC" == "" ]; then
  FC=`which gfortran`; FFLAGS=" -O3 -cpp"
fi

$FC VISH2p1V1.10.0.for PhyBdary-1.10.for InputEOS-1.3.for OSCARoutput.for Arsenal-0.7.for Initialization-1.03.for InputFun-1.29RC6.for -o VISHNew.e $FFLAGS
