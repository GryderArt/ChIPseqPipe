#!/bin/sh
#SBATCH --partition=ccr
#SBATCH --ntasks-per-node=8
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --error=bchn2.err
#SBATCH --output=bchn2.out

module load deeptools
module load bedtools
module load homer

perl /data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts/plotHeatmap.pl -b \
RH4_AllBAFs.bed,RH4_PandBAFonly.bed,RH4_BAFandGBAFonly.bed,RH4_GandPBAF_only.bed,RH4_PBAFonly.bed,RH4_BAFonly.bed,RH4_GBAFonly.bed \
-c BAFheat_rpm.conf -o heatmap_rpm2 -m
