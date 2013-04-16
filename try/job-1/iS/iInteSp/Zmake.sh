#! /usr/bin/env bash

# C compiler
CC=`which icpc`; CFLAGS=" -O3 -Wall"
if [ "$CC" == "" ]; then
  CC=`which g++`; CFLAGS=" -O3 -Wall"
  fi
# use next line on tranditional system:
#CFLAGS="$CFLAGS `gsl-config --cflags --libs`"
# use next line on newest Ubuntu
  CFLAGS="$CFLAGS `gsl-config --cflags --libs` -B/usr/lib/i386-linux-gnu -I/usr/include/i386-linux-gnub"

# FORTRAN compiler
  FC=`which ifort`; FFLAGS=" -O3 -fast -heap-arrays -cpp"
  if [ "$FC" == "" ]; then
    FC=`which gfortran`; FFLAGS=" -O3 -cpp"
    fi

$CC iInteSp.cpp arsenal.cpp Table.cpp -o iInteSp.e $CFLAGS

