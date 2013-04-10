#! /usr/bin/env bash

# C compiler
CC=`which icpc`; CFLAGS=" -O3 -fast"
if [ "$CC" == "" ]; then
  CC=`which g++`; CFLAGS=" -O3"
fi
CFLAGS="$CFLAGS `gsl-config --cflags --libs`"

# FORTRAN compiler
FC=`which ifort`; FFLAGS=" -O3 -fast -heap-arrays -cpp"
if [ "$FC" == "" ]; then
  FC=`which gfortran`; FFLAGS=" -O3 -cpp"
fi

$CC decay.c Table.cpp arsenal.cpp functions.c int.c reso.c tools.c -o resonance.e $CFLAGS
