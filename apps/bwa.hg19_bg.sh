#!/bin/sh
#
###############################################################
#
# Authors Rajesh Patidar: rajbtpatidar@gmail.com and Berkley Gryder: gryderart@gmail.com
# Run BWA on PE Illumina Reads.
# Output is Bam sample with Bam statistics.
# Run Command : 
#		$0 Sample
###############################################################
echo $sample
echo "################"
echo "################"
echo "################"
echo "################"

module load picard/2.17.11
module load bwa/0.7.17
module load samtools/1.8
module load fastqc/0.11.6
workdir=/lscratch/$SLURM_JOB_ID/
#sample=$1
dir=/data/khanlab/projects/ChIP_seq/DATA/$sample
cd $dir
ref=/data/Clinomics/Ref/khanlab/Index/bwa/ucsc.hg19.fasta
bed=/data/khanlab/projects/Genotyping/ALLSamples/bedFile
RG="@RG\tID:${sample}\tLB:${sample}\tSM:${sample}\tPL:ILLUMINA"

cd $dir
echo $dir
if [ -e "${sample}.flagstat.txt" ]
then
	echo "Taking the presence of a flagstat.txt to assume mapping was previously completed."
	echo "   (if you want to override existing BAMs, delete flagstat and relaunch)"
else
	R1=`/bin/ls ${sample}_R1.fastq.gz`
	if [ -e "${sample}_R2.fastq.gz" ]
	then
		echo "Processing in PE mode"
		R2=`/bin/ls ${sample}_R2.fastq.gz*`
		mkdir -p $dir/Fastqc/
		fastqc -t ${SLURM_CPUS_ON_NODE} $R1 -o $dir/Fastqc/
		fastqc -t ${SLURM_CPUS_ON_NODE} $R2 -o $dir/Fastqc/
		unzip $dir/Fastqc/${sample}_R1_fastqc.zip -d $dir/Fastqc/
		unzip $dir/Fastqc/${sample}_R2_fastqc.zip -d $dir/Fastqc/
		echo "Starting BWA"
		bwa mem -M -t ${SLURM_CPUS_ON_NODE} -R "$RG" $ref $R1 $R2 |samtools view -h -q 30 - >$workdir/${sample}.sam
	else
		echo "Processing in SE mode"
		mkdir -p $dir/Fastqc/
		fastqc -t ${SLURM_CPUS_ON_NODE} $R1 -o $dir/Fastqc/
		unzip $dir/Fastqc/${sample}_R1_fastqc.zip -d $dir/Fastqc/
		echo "Starting BWA"
		bwa mem -M -t ${SLURM_CPUS_ON_NODE} -R "$RG" $ref $R1     |samtools view -h -q 30 - >$workdir/${sample}.sam

	fi
	echo "Mapping finished"
	echo "################"
	java -Xmx60g  -Djava.io.tmpdir=$workdir -jar $PICARDJARPATH/picard.jar SortSam INPUT=$workdir/${sample}.sam OUTPUT=$dir/${sample}.bam SORT_ORDER=coordinate
	echo "Sam to bam conversion finished"
	echo "################"

	samtools index $dir/${sample}.bam
	samtools flagstat ${sample}.bam  >${sample}.flagstat.txt
	#java -Xmx60g -Djava.io.tmpdir=$workdir -jar $PICARDJARPATH/picard.jar MarkDuplicates I=${sample}.bam O=${sample}.dd.bam REMOVE_DUPLICATES=false M=${sample}.matrix.txt AS=true VALIDATION_STRINGENCY=LENIENT
	#samtools index ${sample}.dd.bam

	echo -e "Mapping finished on ${sample}.\n\nRegards,\nKhanLab\nOncogenomics Section\nCCR NCI NIH" |mutt -s "Mapping Job Status" ChIP_user 
fi
echo "Launching single sample downstream pipeline."
echo "perl /data/khanlab/apps/ChIPseq/automate_pipeline_bg.pl $sample"
perl /data/khanlab/apps/ChIPseq/automate_pipeline_bg.pl $sample
