#! /usr/bin/env bash
(cd ..
for ii in `ls -p | grep /`
    do
    (cd $ii
        rm *.pyc
        if [ -e GNUmakefile ]; then make distclean; fi
        if [ -e __pycache__ ]; then rm -rf __pycache__; fi
    )
done
)
