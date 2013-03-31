#! /usr/bin/env bash

target_location="PlayGround"

for ii in `ls ./$target_location`
do
	(cd ./$target_location/$ii; qsub *.pbs)
done
