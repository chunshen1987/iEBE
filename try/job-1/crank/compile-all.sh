#! /usr/bin/env bash
(cd ..
for ii in superMC VISHNew iSS iS osc2u urqmd
    do
    (cd $ii; make; make clean)
done

echo "Compiling finished."
echo "Next generate jobs using generate-jobs-XXX.sh."
)
