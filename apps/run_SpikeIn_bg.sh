#!/bin/sh
#
###############################################################
# Xinyu Wen and Berkley Gryder
# run BWA and analysis pipeline for DM-HG ChiPSeq for ChIP-Rx SpikeIn
# Run Command : 
#       $0 Sample
###############################################################
sample=$1

mkdir /data/khanlab/projects/ChIP_seq/DATA/$sample/SpikeIn
ln -s /data/khanlab/projects/DATA/$sample/*.fastq.gz /data/khanlab/projects/ChIP_seq/DATA/$sample/SpikeIn
chgrp khanlab /data/khanlab/projects/ChIP_seq/DATA/$sample/SpikeIn
chmod g+w /data/khanlab/projects/ChIP_seq/DATA/$sample/SpikeIn
workdir=/data/khanlab/projects/ChIP_seq/DATA/$sample/SpikeIn
logdir=$workdir/pbs_log
mkdir $logdir
chgrp khanlab $logdir
chmod g+w $logdir

/usr/local/slurm/bin/sbatch -J SPKI --exclusive --mem=60g --cpus-per-task=30 --partition=ccr --gres=lscratch:400 --workdir=$workdir --output=$logdir/pipe.o --error=$logdir/pipe.e --export=MODULEPATH=/usr/local/lmod/modulefiles,sample=$sample --time=24:00:00 /data/khanlab/apps/ChIPseq/dm3_hg19_v4.sh
