#! /usr/bin/env bash

if [ $# -eq 1 ]
then
    # pion
    echo "Processing pion integrated flows..."
    python ./regulateVnData.py $1 2 spectra/thermal_211_integrated_vndata.dat
    echo "Processing pion differential flows..."
    python ./regulateDiffVnData.py $1 2 spectra/thermal_211_vndata.dat
    # Kaon
    echo "Processing Kaon integrated flows..."
    python ./regulateVnData.py $1 5 spectra/thermal_321_integrated_vndata.dat
    echo "Processing Kaon differential flows..."
    python ./regulateDiffVnData.py $1 5 spectra/thermal_321_vndata.dat
    # proton
    echo "Processing proton integrated flows..."
    python ./regulateVnData.py $1 18 spectra/thermal_2212_integrated_vndata.dat
    echo "Processing proton differential flows..."
    python ./regulateDiffVnData.py $1 18 spectra/thermal_2212_vndata.dat
else
    echo "Usage: extractThermal3.sh directory_path"
    exit
fi
