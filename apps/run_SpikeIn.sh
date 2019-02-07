#!/bin/sh
#
###############################################################
# Xinyu Wen
# run BWA and analysis pipeline for DM-HG ChiPSeq samples by pipeline automation
# Run Command : 
#       $0 Sample
###############################################################
sample=$1
#echo $sample
#echo "################"
#echo "################"

mkdir /data/khanlab/projects/ChIP_seq/DATA/$sample/SpikeIn
chgrp khanlab /data/khanlab/projects/ChIP_seq/DATA/$sample/SpikeIn
chmod g+w /data/khanlab/projects/ChIP_seq/DATA/$sample/SpikeIn
workdir=/data/khanlab/projects/ChIP_seq/DATA/$sample/SpikeIn
logdir=$workdir/pbs_log

/usr/local/slurm/bin/sbatch -J SPKI --exclusive --mem=60g --cpus-per-task=30 --partition=ccr --gres=lscratch:400 --workdir=$workdir --output=$logdir/pipe.o --error=$logdir/pipe.e --export=MODULEPATH=/usr/local/lmod/modulefiles,sample=$sample --time=24:00:00 /data/khanlab/apps/ChIPseq/dm3_hg19_v4.sh
