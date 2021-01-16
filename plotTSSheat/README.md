# Plot Heatmaps by Genesets

1. Plot gene set heatmap with plotBEDheat.sh 

	build a configuration using codebuilder.xlsx (tab: heatmap)
	
	copy into config file (filename.conf), defining BAM locations, name and color
	
	edit plotBEDheat.sh to use the gene sets of interest, parameters, and config file name
	
2. Extract and separate BED coordinates with filterGeneset_byHeatmapRank.R
