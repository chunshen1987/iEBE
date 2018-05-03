#! /usr/bin/env bash

CC=`which icpc`;
CFLAGS=" -O3"
if [ "$CC" == "" ]; then
    CC=`which g++`;
    CFLAGS=" -O3"
fi
CFLAGS="$CFLAGS -Wall"
echo $CFLAGS
