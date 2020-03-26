#!/bin/bash

module load bedtools

#example: sbatch --partition=ccr --time=24:00:00 --export=sample_list=../sample_list.txt ../scripts/calculateTR.sh

data_home='../ChIP_seq/DATA'
script_home='../scripts'

IFS=$'\n' read -d '' -r -a samples < $sample_list

for sample in "${samples[@]}"
do
	echo "$script_home/calculateTR.pl -b $data_home/$sample/$sample.bam -o $sample.trv"
	$script_home/calculateTR.pl -b $data_home/$sample/$sample.bam -o $sample.trv
done	
