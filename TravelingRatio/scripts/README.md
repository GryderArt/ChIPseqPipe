# Protocol for calculating RNA Pol2 Pause Index
by Berkley E. Gryder, March 2020

### (1)	Generate high quality ChIP-seq data for RNA-Pol2, with and without a perturbation of interest (a drug or knockout of a gene of interest).
    [use on a DATA directory with the /DATA/SampleName/SampleName.bam file structure as here https://github.com/CBIIT/ChIP_seq.]
### (2)	Generate .TRV files containing read counts along each gene body region, all genes.
    a.	Using codebuilder, go to the “TR” (for Traveling Ratio) tab.  
    b.	Create a list of samples for which you want to calculate the TR (as a 1 column text file, one sample per row).
    c.	On the command line (biowulf unix login node), navigate to output folder with the TR text file
    d.	run the wrapper calculateTR.sh shell script as a sbatch, which executes calculateTR.pl perl script once per sample.  Example:
      i.	sbatch --partition=ccr --export=sample_list=TRsample.list.txt ../scripts/calculateTR.sh
### (3)	Calculate and plot cumulative TR across % of genes
    a.	Load samples into R script (plotTRV.R) using R studio.  This allows for filtering to smaller TR files, based on either a gene set of interest, or based on the presence of Pol2 in the promoter (to exclude calculating a TR for genes the have no Pol2).  
    b.	Write new filtered *.trv files from R.
    c.	On the command line, run the cumulateValue.pl script on all *.trv files as follows:
      i.	sinteractive
      ii.	find -name '*.trv' -exec /data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts/cumulateValue.pl -i {} -o {}.out \;
    d.	Once complete, return to the R studio session and load in the .trv.out files
    e.	Plot with ggplot, save as PDF file, clean up in Adobe Illustrator
