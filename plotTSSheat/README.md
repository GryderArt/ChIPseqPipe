# Plot Heatmaps by Genesets
added by Berkley E. Gryder, Janurary 2021

Scripts, code protocol and video tutorial for plotting ChIP-seq data around gene sets of interest
https://github.com/GryderArt/ChIPseqPipe/tree/master/plotTSSheat

<img src="THZ1_genes_heatmap.png" width="700"/></a>

# Basic steps to follow 

1. Download scripts and example files, edit file paths in scripts to suit your own

2. Plot gene set heatmap with plotBEDheat.sh 

	build a configuration using codebuilder.xlsx (tab: heatmap)
	
	copy into config file (filename.conf), defining BAM locations, name and color
	
	edit plotBEDheat.sh to use the gene sets of interest, parameters, and config file name
	
3. Extract and separate BED coordinates with filterGeneset_byHeatmapRank.R

# Tutorial video for example usage

https://drive.google.com/file/d/19CGKYSLPgO8lhpsEJNvl0ItHgUV4gv3v/view
