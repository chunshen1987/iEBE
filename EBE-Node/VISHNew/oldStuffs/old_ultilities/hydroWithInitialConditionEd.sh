#! /usr/bin/env bash

initial_file=InitialEd.dat
# check argument and get directories
if [ $# -le 2 ]
then
  echo "Usage: hydroWithInitialCondition.sh source_folder destination_folder hydro_executable and parameters"
  exit
fi

para=""
while [ $# -gt 0 ]
do
  arg=`echo "$1" | sed 's/ /\?/g'`
  para="$para $arg"
  shift
done
echo $para
python ./hydroWithInitialCondition.py $initial_file $para

