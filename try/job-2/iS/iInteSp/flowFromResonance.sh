#! /usr/bin/env bash

# Perform flow calculation for all subfolders under the target folder, which is the first argument. Folders do not contain the skip_if_missing file will be skipped; folders already contain the skip_if_exist file will be skipped. The subfolder I refer to is the subfolder that has the structure from the "spectraFromHydro.sh" script from Azspectra; meaning that all the dN/d3p data are stored under "target_folder/subfolder/spectra/".

# Example:
# ./flowFromResonance.sh ~/Downloads/test/

renamer () # add string $1 in front of those files that have eta in their names
# used to be replace string "eta" by $1 to files whose names contain "eta"
{
    (cd ./results;
    for file in `ls *eta*.dat`
    do
        # file_mod1=`echo $file | sed 's/\(.*\)eta\(.*\)/\1/'`
        # file_mod2=`echo $file | sed 's/\(.*\)eta\(.*\)/\2/'`
        # mv $file "$file_mod1"$1$file_mod2
        mv $file $1$file
    done;
    )
}

mover () # move result files from $1 to $2
{
    mv $1/*Charged* $2
}

operation_folder="spectra"
skip_if_missing="spec_211.dat"
skip_if_exist="R_0.5_Charged_eta_vndata.dat"

table_folder="tables"
executable="./iInteSp.e"

current_dir=`pwd`

# check argument and get directories
if [ $# -lt 1 ]
then
  echo "Usage: flowFromResonance.sh directory"
  exit
else
  targetDir=$1
fi

rm ./results/*

for sub_folder in `ls $targetDir`
do
    echo
    echo Dealing with $sub_folder ...
	if [ -e $targetDir/$sub_folder/$operation_folder/$skip_if_exist ]
    then
        echo "--- Files $skip_if_exist already exists, skipped."
        continue
    fi
    if [ -e $targetDir/$sub_folder/$operation_folder/$skip_if_missing ]
    then
        :
    else
        echo "--- Files $skip_if_missing does not exist, skipped."
        continue
    fi
    cp $targetDir/$sub_folder/$operation_folder/spec* ./results
    mv $current_dir/$table_folder/particle_eta_table.dat $current_dir/$table_folder/particle_eta_table.dat.tmp
    
    # use different eta tables
    
    # |eta|<0.5:
    mv $current_dir/$table_folder/particle_eta_uni_0.5_table.dat $current_dir/$table_folder/particle_eta_table.dat
    $executable
    mv $current_dir/$table_folder/particle_eta_table.dat $current_dir/$table_folder/particle_eta_uni_0.5_table.dat
    renamer "R_0.5_" # R for "restricted to"
    mover $current_dir/results $targetDir/$sub_folder/$operation_folder/
    
    # |eta|<2.5:
    mv $current_dir/$table_folder/particle_eta_uni_2.5_table.dat $current_dir/$table_folder/particle_eta_table.dat
    $executable
    mv $current_dir/$table_folder/particle_eta_table.dat $current_dir/$table_folder/particle_eta_uni_2.5_table.dat
    renamer "R_2.5_"
    mover $current_dir/results $targetDir/$sub_folder/$operation_folder/
    
    # |eta|<2.0:
    mv $current_dir/$table_folder/particle_eta_uni_2.0_table.dat $current_dir/$table_folder/particle_eta_table.dat
    $executable
    mv $current_dir/$table_folder/particle_eta_table.dat $current_dir/$table_folder/particle_eta_uni_2.0_table.dat
    renamer "R_2.0_"
    mover $current_dir/results $targetDir/$sub_folder/$operation_folder/

    # 0.5<|eta|<2.5:
    mv $current_dir/$table_folder/particle_eta_uni_0.5_2.5_table.dat $current_dir/$table_folder/particle_eta_table.dat
    $executable
    mv $current_dir/$table_folder/particle_eta_table.dat $current_dir/$table_folder/particle_eta_uni_0.5_2.5_table.dat
    renamer "R_0.5_2.5_"
    mover $current_dir/results $targetDir/$sub_folder/$operation_folder/

    # 0.5<|eta|<4.8:
    mv $current_dir/$table_folder/particle_eta_uni_0.5_4.8_table.dat $current_dir/$table_folder/particle_eta_table.dat
    $executable
    mv $current_dir/$table_folder/particle_eta_table.dat $current_dir/$table_folder/particle_eta_uni_0.5_4.8_table.dat
    renamer "R_0.5_4.8_"
    mover $current_dir/results $targetDir/$sub_folder/$operation_folder/
    
    # restore default eta table
    mv $current_dir/$table_folder/particle_eta_table.dat.tmp $current_dir/$table_folder/particle_eta_table.dat
    echo Done.
    rm ./results/*
done


