#! /usr/bin/env bash

if [ $# -lt 1 ]
then
    echo "Usage: remedy.sh folder"
    exit
fi

cwd=`pwd`
folder=$1

(
cd $folder
for aZip in `ls *.zip`
do
    name=`echo $aZip | cut -d "." -f 1`
    unzip $aZip
    $cwd/EbeCollectorShell_pureHydro.py ./$name
    zip -rm $aZip ./$name/collected.db
    rm -rf $name
done
)
