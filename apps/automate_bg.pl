#!/usr/local/bin/perl
use strict;
use warnings;
use File::Basename;
#use Excel::Writer::XLSX;
#use Spreadsheet::Read;
#my $DIR = "/data/khanlab/projects/";
my $DIR = "/data/khanlab/projects/ChIP_seq/DATA/";
my $TOOL = "/data/khanlab/apps/ChIPseq/bwa.hg19_bg.sh";
my $TOOL_1="/data/khanlab/apps/ChIPseq/bwa.mm9_xy.sh";
my $FILE = "/data/khanlab/apps/ChIPseq/SampleList";
my $BERKLEY ="berkley.gryder\@nih.gov";
my $batch="/usr/local/slurm/bin/sbatch -J ChIPseq --exclusive --mem=60g --cpus-per-task=30 --partition=ccr --gres=lscratch:100 --time=24:00:00";

open(FH, $FILE);
my %Done;
while(<FH>){
  chomp;
  $Done{"$_"} = "";
}
my @All;
my @FOLDERS = <$DIR*>;
foreach my $folder(@FOLDERS){
  if( -d $folder && -w $folder){
	my $sample=basename($folder);
	opendir(DH, "$folder") || die;
	while(readdir DH){
	  if($_ =~ /(.*)_R1.fastq.gz/ and -R "$folder/$_" ){
		if($1 eq $sample){
	  	  if(! exists $Done{$sample}){
			if($sample =~ /Sample_mm/){
			  print "mm\t\t$sample\n";
			  my $map_job = `$batch --workdir=$folder --output=$folder/pbs_log/pipe.o --export=MODULEPATH=/usr/local/lmod/modulefiles,sample=$sample $TOOL_1`;
			  `echo \"ChIP-Seq Mapping job for sample $sample started\" \|mutt -s \"ChIP-Seq Mapping job started \#$map_job\" ChIP_user`;
			}
			else{
	 		  print "human\t\t$sample\n";
			  my $map_job = `$batch --workdir=$folder --output=$folder/pbs_log/pipe.o --export=MODULEPATH=/usr/local/lmod/modulefiles,sample=$sample $TOOL`;
             `echo \"ChIP-Seq Mapping job for sample $sample started\" \|mutt -s \"ChIP-Seq Mapping job started \#$map_job\" ChIP_user`;
			}
			`echo $sample >>$FILE`;

	   	  }
	    }
	  }
    }
  }
}
