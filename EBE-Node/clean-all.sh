#! /usr/bin/env bash

for ii in `ls -p | grep /`
    do
    (cd $ii
        rm *.pyc
        if [ -e GNUmakefile ]; then make distclean; fi
    )
done
