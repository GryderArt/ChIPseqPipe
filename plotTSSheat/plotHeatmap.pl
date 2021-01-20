#!/usr/bin/perl -w

#example:
#
# perl /data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts/plotHeatmap.pl -b /data/khanlab/projects/ChIP_seq/DATA/Sample_J0123_T1_ATAC_1_SCH_C_HMGFJBGX7/MACS_Out_p_1e-07/Sample_J0123_T1_ATAC_1_SCH_C_HMGFJBGX7_summits.bed,/data/khanlab/projects/ChIP_seq/DATA/Sample_J0193_T1_ATAC_1_SCH_C_HMGFJBGX7/MACS_Out_p_1e-07/Sample_J0193_T1_ATAC_1_SCH_C_HMGFJBGX7_summits.bed -c heatmap.conf -l J0123,J0193
#
use strict;
use Getopt::Long qw(GetOptions);
use Time::Piece;
use File::Basename;
use Cwd 'abs_path';


my $bed_list;
my $bed_labels;
my $config_file;
my $output_dir="heatmap_out";
my $ext_bp=2000;
my $recalculate=0;
my $not_spikein=0;
my $num_processors=8;
my $bin_size = 25;
my $smooth_length = 75;
my $dpi = 600;
my $type = "bed";
my $sort_sample = 1;
my $tss_file = "/data/khanlab/projects/ChIP_seq/data_by_file_type/bed/genes_refseq.hg19.bed";
my $usage = <<__EOUSAGE__;

Usage:

$0 [options]

required options:

  -b	<string> BED/Gene list (comma seperated)
  -c	<string> Config file  
 
 optional options:
  -t	type (bed or gene) default ($type)
  -r	recalculate the matrix file
  -m	do not use spikeIn data
  -l	<string> BED label list (comma seperated)  
  -o	<string> Output directory (default: $output_dir)
  -e	<int>	extension length (default: $ext_bp)
  -p	<int>	number of processors (default: $num_processors)
  -n	<int>	bin size of bigwig file (default: $bin_size)
  -s	<int>	smooth length of bigwig file (default: $smooth_length)
  -g	<int>	sort by sample (default: $sort_sample)
  -d	<int>	dpi (default: $dpi)
  
__EOUSAGE__



GetOptions (
  'b=s' => \$bed_list,
  'l=s' => \$bed_labels,
  't=s' => \$type,
  'c=s' => \$config_file,
  'o=s' => \$output_dir,
  'e=i' => \$ext_bp,
  'n=i' => \$bin_size,
  's=i' => \$smooth_length,
  'g=i' => \$sort_sample,
  'r' => \$recalculate ,
  'm' => \$not_spikein
);


unless ($bed_list && $config_file) {
    die "Please input BED/gene list and config files\n$usage";
}

if ($type ne "bed" && $type ne "gene") {
	die "Unknown type: $type. Please use 'bed' or 'gene'\n";
}
if (!$bed_labels) {
	my @labels = ();
	my @arr = split(/,/, $bed_list);
	foreach my $f (@arr) {
		my $fn = basename($f);
		push @labels, $fn;
	}
	$bed_labels = join(",", @labels);
}


#Global variables

$output_dir = &formatDir($output_dir);
&runCommand("mkdir -p $output_dir");
if ($type eq "gene") {
	my @gene_lists = split(/,/, $bed_list);
	my @out_list = ();
	foreach my $gene_list(@gene_lists) {
		my $out_name = $output_dir."/".basename($gene_list)."_tss.bed";
		push @out_list, $out_name;
		&runCommand("dos2unix $gene_list");
		my %genes = ();
		open(GENE_FILE, $gene_list) or die "Cannot open file $gene_list";
		while(<GENE_FILE>) {
			chomp;
			$genes{$_} = "";
		}
		close(GENE_FILE);
		open(TSS_OUT_FILE, ">$out_name") or die "Cannot open file $out_name";
		open(TSS_FILE, $tss_file) or die "Cannot open file $tss_file";
		while(<TSS_FILE>) {
			chomp;
			my @fields = split(/\t/);
			if (exists $genes{$fields[3]}) {
				print TSS_OUT_FILE $_."\n";
			}			
		}
		close(TSS_FILE);
		close(TSS_OUT_FILE);		
	}
	$bed_list = join(",", @out_list);
}
$bed_list =~ s/,/ /g;
$bed_labels =~ s/,/ /g;
#read config file
my @bws = ();
my @bwlabels = ();
my @colors = ();
my @min_values = ();
my @max_values = ();

open(CONFIG_FILE, "$config_file") or die "Cannot open file $config_file";
my $has_max = 0;
while (<CONFIG_FILE>) {
	chomp;
	my @tokens = split(/\t/);
	next if ($#tokens < 3);
	my $bw_file = $tokens[0];
	my $bw_label = $tokens[1];
	my $ext_len = $tokens[2];
	my $color = $tokens[3];
	my $max_value = 0;
	my $min_value = 0;
	if ($#tokens >= 4) {
		$has_max = 1;
		my @values = split(/,/, $tokens[4]);
		if ($#values == 1) {
			$min_value = $values[0];
			$max_value = $values[1];
		} else {
			$max_value = $tokens[4];
		}
	}
	my ($main, $ext) = $bw_file =~ /(.*)\.([^.]+)$/;
	if ($ext eq "bam") {
		my $bam_file = $bw_file;
		$bw_file = "${main}.${bin_size}.RPM.bw";
		if (!-s $bw_file) {
			my $cmd = "bamCoverage -e $ext_len -b $bam_file -o $bw_file --outFileFormat bigwig --smoothLength $smooth_length --binSize ${bin_size} --normalizeUsing CPM -p $num_processors";
			&runCommand($cmd);
			&runCommand("chgrp khanlab $bw_file");
		}
		my $spike_in_file = dirname($bw_file)."/SpikeIn/spike_map_summary";
		if (-s $spike_in_file && !$not_spikein) {			
			open(SPIKE_IN_SUMMARY, $spike_in_file) or die "Cannot open file $spike_in_file";
			<SPIKE_IN_SUMMARY>;
			my $summary_line = <SPIKE_IN_SUMMARY>;
			my @summary_fields = split(/\t/, $summary_line);
			my $scale_factor_raw = 1000000/$summary_fields[2];
			my $scale_factor  = sprintf "%.4f", $scale_factor_raw;
			$bw_file = "${main}.${bin_size}.scaled.bw";
			if (!-s $bw_file) {
				my $cmd = "bamCoverage -e $ext_len -b $bam_file -o $bw_file --outFileFormat bigwig --smoothLength $smooth_length --binSize ${bin_size} --scaleFactor $scale_factor -p $num_processors";
				&runCommand($cmd);
				&runCommand("chgrp khanlab $bw_file");
			}
		}
	}
	push @bws, $bw_file;
	push @bwlabels, $bw_label;
	push @colors, $color;
	push @min_values, $min_value;
	push @max_values, $max_value;
}

my $bw_list = join(' ', @bws);
my $bwlabel_list = join(' ', @bwlabels);
my $color_list = join(' ', @colors);
my $ymax_option = "";
my $zmax_option = "";
if ($has_max) {
	$ymax_option = "--yMin 0 --yMax ".join(' ', @max_values);
	$zmax_option = "--zMin ".join(' ', @min_values)." --zMax ".join(' ', @max_values);
}
my $out_matrix_file = "${output_dir}matrix.txt";
my $matrix_zfile = "${output_dir}matrix.gz";
my $matrix_file = "${output_dir}matrix.tab";
my $heatmap_file = "${output_dir}heatmap.pdf";
my $heatmap_bed_file = "${output_dir}heatmap.bed";
my $profile_file = "${output_dir}profile.pdf";
if (!-e $matrix_zfile || $recalculate) {
	my $ref_point = ($type eq "gene")? "TSS" : "center";
	my $cmd = "computeMatrix reference-point -R $bed_list -S $bw_list --referencePoint $ref_point -a $ext_bp -b $ext_bp --samplesLabel $bwlabel_list -p $num_processors -o $matrix_zfile --outFileNameMatrix $matrix_file";	
	#print ("$cmd\n");
	&runCommand($cmd);
}
my $ref_label = ($type eq "gene")? "TSS" : "Peak";
my $cmd = "plotHeatmap -m $matrix_zfile --missingDataColor white --colorList $color_list --regionsLabel $bed_labels -o $heatmap_file --refPointLabel $ref_label --dpi $dpi --sortUsingSamples $sort_sample --outFileSortedRegions $heatmap_bed_file --interpolationMethod nearest $ymax_option $zmax_option";
#print("$cmd\n");
&runCommand($cmd);
$cmd = "plotProfile -m $matrix_zfile --numPlotsPerRow 3 --regionsLabel $bed_labels -o $profile_file --refPointLabel $ref_label $ymax_option";
&runCommand($cmd);
&runCommand("chgrp khanlab ${output_dir}*");
sub formatDir {
    my ($dir) = @_;
    if ($dir !~ /\/$/) {
        $dir = $dir."/";
    }
    return $dir;
}

sub runCommand {
    my ($cmd) = @_;
	print $cmd."\n";
    system($cmd);
}

sub runCommandToSTDOUT {
    my ($cmd, $stdout_file, $stderr_file) = @_;
    open OUT_FILE, ">$stdout_file" || print STDERR "cannot create file: $stdout_file";
    open ERR_FILE, ">$stderr_file" || print STDERR "cannot create file: $stderr_file";
    open (CMD, "($cmd | sed 's/^/STDOUT:/') 2>&1 |");
    while (<CMD>) {
         if (s/^STDOUT://)  {
             print OUT_FILE $_;
         } else {
             print ERR_FILE $_;
         }
    }
    close(OUT_FILE);
    close(ERR_FILE);

}
