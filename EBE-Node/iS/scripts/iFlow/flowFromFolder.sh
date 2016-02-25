#! /usr/bin/env bash

# Calculate flow using the specified dN/(pT dpT dphi dy) matrix, then copy results to the spcified folder.

if [ $# -lt 3 ]
then
    echo Usage: flowFromFolder.sh dN_dy_matrix_path_and_filename copy_result_to_folder result_name
    exit
fi

dNdy_file=$1
result_folder=$2
result_name=$3

if [ $# -gt 4 ]
then
    pT_cut_min=$4
    pT_cut_max=$5
else
    pT_cut_min=""
    pT_cut_max=""
fi

./iFlow.e $dNdy_file $result_name $pT_cut_min $pT_cut_max
mv ./results/$result_name* $result_folder


