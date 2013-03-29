#! /usr/bin/env bash

for ii in superMC VISHNew iSS osc2u urqmd
    do
    (cd $ii; make)
done

echo "Compiling finished."
echo "Next generate jobs using generate-jobs-XXX.sh."
