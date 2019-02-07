#!/usr/local/bin/perl
#Xinyu Wen - 201810
use strict;
use warnings;

my $sample = $ARGV[0];
chomp $sample;
print "==>Prepare ChIP-Seq pipeline for sample $sample\n";
my $master_dir = "/data/khanlab/projects/ChIP_seq/manage_samples";
my $master_txt = "ChIP_seq_samples_pipeline.txt";
my $meta_dir = $master_dir."/single_sample_config";
my $meta_txt = $sample."_metadata.txt";
my $conf_template = "Sample_TEMPLATE.conf";
my $conf_txt = $sample.".conf";
my $data_dir = "/data/khanlab/projects/ChIP_seq/DATA";
my $code_dir = "/data/khanlab/apps/ChIPseq";
my $pipe_code = "/data/khanlab/projects/ChIP_seq/scripts/runChipseqPipeline.sh";
my $spikein_code = "run_SpikeIn_bg.sh";
my %ref;
my $line = `head -1 $master_dir\/$master_txt`;
if (defined $line){
  print "==>Parse master file $master_dir\/$master_txt\n";
}
else{
  `echo \"Master file $master_dir\/$master_txt cannot be reached. Please check!\"\|mutt -s \"ChIP-Seq pipeline failed on $sample\" ChIP_user`;
  die "==>Cannot reach master file $master_dir\/$master_txt\! Quit";
}
my $slurm_log_dir = "/data/khanlab/projects/ChIP_seq/manage_samples/slurm_log";

chomp $line;
$line =~s/\s+$//;
$line =~s/^\s+//;
my @ref_header_tmp = split(/\t/,$line);
my @ref_header = ();
foreach my $header (@ref_header_tmp){
  if ($header){
    push @ref_header, $header;
  }
  else{
    push @ref_header, "NA";
  }
}
#set header
foreach my $header (@ref_header){
  $ref{$header} = "";
}
my @a = split(/\_/,$sample);
my $fcid = pop@a;
my $sample_ref = join("_",@a);
$line = "";
$line = `grep $sample_ref $master_dir\/$master_txt 2>/dev/null`;
@a= ();
#assign value to header
if (defined $line){
  print "==>Entry found in master file for $sample\n";
  $line =~s/\s+$//;
  chomp $line;
  $line =~s/TBD/$fcid/g;
  @a = split(/\t/,$line);
  for (my $i = 0;$i<scalar@a;$i++){
    if ($a[$i]){
	  my $entry = $a[$i];
	  $entry =~s/^\s+//;
	  $entry =~s/\s+$//;
      chomp $entry;
      $ref{$ref_header[$i]} = $entry; 
    }
    else{
      $ref{$ref_header[$i]} = "";
    }
  }
}
print "==>Making meta data file $meta_dir\/$meta_txt\n"; 
open my $fo,">$meta_dir/$meta_txt";
for (my $i = 0; $i<scalar@ref_header; $i++){
  if ($i == (scalar@ref_header -1)){
    print $fo $ref_header[$i];
  }
  else{
    print $fo $ref_header[$i]."\t";
  }
}
print $fo "\n";
for (my $i = 0; $i<scalar@ref_header; $i++){
  if ($i == (scalar@ref_header - 1)){
    print $fo $ref{$ref_header[$i]};
  }
  else{
    print $fo $ref{$ref_header[$i]}."\t";
  }
}
print $fo "\n";
print "==>Making conf file $meta_dir\/$conf_txt\n";
open $fo,">$meta_dir/$conf_txt";
open my $fi,"$meta_dir/$conf_template";
while (<$fi>){
  if (m/^sample\_desc\_file\=manage\_samples\/single\_sample\_config\//){
    print $fo "sample_desc_file=manage_samples/single_sample_config/".$meta_txt."\n";
  }
  else{
    print $fo $_;
  }
}
close $fo;
close $fi;
print "==>Change meta and conf file group to khanlab\n";
`chgrp khanlab $meta_dir\/$meta_txt`;
`chgrp khanlab $meta_dir\/$conf_txt`;
print "==>Check SpikeIn status (yes/no): ".$ref{SpikeIn}."\n";
if ($ref{SpikeIn} eq "yes"){
  if (-e "$data_dir\/$sample\/SpikeIn/spike_map_summary"){
    print "==>$data_dir\/$sample\/SpikeIn/spike_map_summary present\n";
  }
  else{					
	  my $spikein_job_id = `$code_dir\/$spikein_code $sample`;
	  chomp $spikein_job_id;
	  print "==>No previous spike_map_summary detected.\n Launched SpikeIn Job ID is: $spikein_job_id\n";
	  `touch $data_dir\/$sample\/SpikeIn/spikein_job_started.txt`;
	  `echo "$spikein_job_id" > $data_dir\/$sample\/SpikeIn/spikein_job_started.txt`;
		my $ct_spkn = 1;
	    while (1){
	      sleep (1800);
		  if (-e "$data_dir\/$sample\/SpikeIn/spike_map_summary"){
              print "==>$data_dir\/$sample\/SpikeIn/spike_map_summary completed!\n";
			  last;
			}
		  else{
            if ($ct_spkn > 14){
              print "==>No $data_dir\/$sample\/SpikeIn/spike_map_summary present, after waiting 7 hours for \#$spikein_job_id\. Give up!\n";
			  `echo \"Cannot find $data_dir\/$sample\/SpikeIn/spike_map_summary, after pipeline job \#$spikein_job_id\. Please check!\"\|mutt -s "ChIP-Seq pipeline alert" ChIP_user`;
			  last;
			}
		  }
	      $ct_spkn++;	print "==> waiting on $data_dir\/$sample\/SpikeIn/spike_map_summary, $ct_spkn\n";
		}
	  }
	}
print "==>Check BAM for ".$sample."\n";
if (-e "$data_dir\/$sample\/$sample\.bam" and -e "$data_dir\/$sample\/$sample\.bam.bai"){
  print "==>BAM present for sample ".$sample."\n";
  print "==>Check PairedInput status: ".$ref{PairedInput}."\n";
  if ($ref{PairedInput} eq "."){
    `chgrp khanlab $data_dir\/$sample\/*.bam*`;
	print "==>Run single sample pipeline\n";
    chdir $meta_dir;
    my $pipeline_job = `sbatch -J singleChIP -e $slurm_log_dir\/$sample\.e -o $slurm_log_dir\/$sample\.o --partition=ccr --time=24:00:00 --mem=121g --cpus-per-task=4  --gres=lscratch:200 --export=config_file=$meta_dir\/$conf_txt $pipe_code`;
	`echo \"ChIP-Seq pipeline launched for sample $sample\" \| mutt -s \"ChIP-Seq pipeline started \#$pipeline_job\" ChIP_user`;
  }
  else{
	print "==>Check BAM for ".$ref{PairedInput}."\n";
	if (-e "$data_dir\/$ref{PairedInput}\/$ref{PairedInput}\.bam" and "$data_dir\/$ref{PairedInput}\/$ref{PairedInput}\.bam.bai"){
	  print "==>BAM present for paired sample ".$ref{PairedInput}."\n";
      `chgrp khanlab $data_dir\/$sample\/*.bam*`;
	  print "==>Run paired sample pipeline\n";
      chdir $meta_dir;
      my $pipeline_job = `sbatch -J singleChIP -e $slurm_log_dir\/$sample\.e -o $slurm_log_dir\/$sample\.o --partition=ccr --time=24:00:00 --mem=121g --cpus-per-task=4  --gres=lscratch:200 --export=config_file=$meta_dir\/$conf_txt $pipe_code`;
	  `echo \"ChIP-Seq pipeline launched for sample $sample\" \| mutt -s \"ChIP-Seq pipeline started \#$pipeline_job\" ChIP_user`;
	}
	else{
	  print "==>No BAM present, wait.\n"; # 1 hour\n";
	  my $ct = 0;
	  while (1){
        sleep(1800);
        if (-e "$data_dir\/$ref{PairedInput}\/$ref{PairedInput}\.bam" and "$data_dir\/$ref{PairedInput}\/$ref{PairedInput}\.bam.bai"){
		  print "==>BAM present for paired sample ".$ref{PairedInput}."\n";
          `chgrp khanlab $data_dir\/$sample\/*.bam*`;
		  print "==>Run paired sample pipeline\n";
          chdir $meta_dir;
          my $pipeline_job = `sbatch -J singleChIP -e $slurm_log_dir\/$sample\.e -o $slurm_log_dir\/$sample\.o --partition=ccr --time=24:00:00 --mem=121g --cpus-per-task=4  --gres=lscratch:200 --export=config_file=$meta_dir\/$conf_txt $pipe_code`;
		`echo \"ChIP-Seq pipeline launched for sample $sample\" \| mutt -s \"ChIP-Seq pipeline started \#$pipeline_job\" ChIP_user`;
	    }
		if ($ct > 10){
          `echo \"No paired input bam for $ref{PairedInput} found, after 5 hours of waiting. The pipeline quit!\n\" \| mutt -s \"ChIP-Seq pipeline process for sample $sample ended with error\" ChIP_user`;
		  die "No paired input bam for $ref{PairedInput} found!";
		}
		$ct++;
	  }
	  #else{
	  #  die "No paired input bam for $ref{PairedInput} found!";
	  #}
	}
  }
}
else{
  `echo \"BAM file for sample $sample is not present, please check!\"\|mutt -s \"ChIP-Seq pipeline failed on $sample\" ChIP_user`;
  die "No Bam found, give up!"; 
}
