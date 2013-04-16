#! /usr/bin/env bash

# C compiler
CC=`which icpc`; CFLAGS=" -O3 -fast -Wall"
if [ "$CC" == "" ]; then
  CC=`which g++`; CFLAGS=" -O3 -Wall"
  fi
# use next line on tranditional system:
#CFLAGS="$CFLAGS `gsl-config --cflags --libs`"
# use next line on newest Ubuntu
    CFLAGS="$CFLAGS `gsl-config --cflags --libs` -B/usr/lib/i386-linux-gnu -I/usr/include/i386-linux-gnu"

# FORTRAN compiler
FC=`which ifort`; FFLAGS=" -O3 -fast -heap-arrays -cpp"
if [ "$FC" == "" ]; then
  FC=`which gfortran`; FFLAGS=" -O3 -cpp"
fi


$FC osc2u.f pdg2ityp.f vni_procev.f blockres.f dectim.f gnuranf.f -o osc2u.e


