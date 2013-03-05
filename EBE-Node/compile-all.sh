#! /usr/bin/env bash

for ii in `ls -p | grep /`
    do
    (cd $ii; bash ./Zmake.sh)
done

echo
echo
echo "Compiling finished."
echo "Next generate jobs using generate-jobs-XXX.sh."