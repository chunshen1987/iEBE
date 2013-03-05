#! /usr/bin/env bash

for ii in `ls -p | grep /`
    do
    (cd $ii; rm *.e)
done
