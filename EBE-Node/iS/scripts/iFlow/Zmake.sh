#! /usr/bin/env bash

# C compiler
CC=`which icpc`; CFLAGS=" -fast -O3 -Wall"
if [ "$CC" == "" ]; then
  CC=`which g++`; CFLAGS=" -O3 -Wall"
  fi

# FORTRAN compiler
  FC=`which ifort`; FFLAGS=" -O3 -fast -heap-arrays -cpp"
  if [ "$FC" == "" ]; then
    FC=`which gfortran`; FFLAGS=" -O3 -cpp"
    fi

$CC iFlow.cpp arsenal.cpp Table.cpp -o iFlow.e $CFLAGS

