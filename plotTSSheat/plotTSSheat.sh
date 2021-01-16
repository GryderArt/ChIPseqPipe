#!/bin/sh
#SBATCH --partition=norm
#SBATCH --ntasks-per-node=8
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --error=bchn2.err
#SBATCH --output=bchn2.out

module load deeptools
module load bedtools
module load homer

perl /data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts/plotHeatmap.pl -b \
CDK7downGenes.genelist.txt,CDK7nochangeGenes.genelist.txt,CDK7upGenes.genelist.txt \
-c plotTSSheat.conf -o H3K4me3_viaCDK7 -t gene -l THZ1_down,THZ1_unchanged,THZ1_up -m

#Usage:
#
#/data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts/plotHeatmap.pl [options]
#
# required options:
#
#  -b    <string> BED/Gene list (comma seperated)
#  -c    <string> Config file
#
# optional options:
#  -t    type (bed or gene) default (bed)
#  -r    recalculate the matrix file
#  -m    do not use spikeIn data
#  -l    <string> BED label list (comma seperated)
#  -o    <string> Output directory (default: heatmap_out)
#  -e    <int>   extension length (default: 2000)
#  -p    <int>   number of processors (default: 8)
#  -n    <int>   bin size of bigwig file (default: 25)
#  -s    <int>   smooth length of bigwig file (default: 75)
#  -g    <int>   sort by sample (default: 1)
#  -d    <int>   dpi (default: 600)

