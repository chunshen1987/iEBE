#! /usr/bin/env bash

# C compiler
CC=`which icpc`; CFLAGS=" -O3"
if [ "$CC" == "" ]; then
  CC=`which g++`; CFLAGS=" -O3"
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

$CC main.cpp Bases.cpp MCnucl.cpp GlueDensity.cpp MakeDensity.cpp KLNModel.cpp OverLap.cpp Largex.cpp Regge96.cpp rcBKfunc.cpp MathBasics.cpp ParameterReader.cpp arsenal.cpp EOS.cpp GaussianNucleonsCal.cpp NBD.cpp RandomVariable.cpp TableFunction.cpp Table.cpp -o superMC.e $CFLAGS
