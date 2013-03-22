#! /usr/bin/env bash

CC=`which icpc`
if [ "$CC" == "" ]; then
   CC=`which g++`;
fi
echo $CC
