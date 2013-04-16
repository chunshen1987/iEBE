#! /usr/bin/env bash

if [ $# -lt 2 ]
then
  echo 'Usage: addPrefix prefix "filenames"'
fi

for i in `ls $2`
do
  cp $i $1-$i
done