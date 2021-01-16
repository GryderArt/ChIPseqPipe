##############################################
## subset gene set by heatmap.bed output    ##
## Berkley Gryder, gryderart@gmail.com 2020 ##
##############################################

#### 1. Prepare environmental settings
  setwd("K:/projects/ChIP_seq/projects/KDM5A/CDK7_v_KDM5A/")

#### 2. Load in heatbed file, define rank based cut off
  
  heat.bed <- read.table("heatmap.bed",sep="\t", header=F)
  colnames(heat.bed) = c("chrom","start","end","name","score","strand","thickStart","thickEnd","itemRGB","blockCount","blockSizes","blockStart","deepTools_group")
  
  #cut by each group
  deepToolsgroups = unique(heat.bed$deepTools_group)
  cutbottomprcntile = 0.43 #express percentile as a fraction of 1
  lapply(deepToolsgroups, function(x) {
    heat.temp <- subset(heat.bed, heat.bed$deepTools_group %in% x)
    heat.temp$rank = as.numeric(row.names(heat.temp))
    heat.temp$rank = heat.temp$rank - heat.temp$rank[1] +1
    heat.temp.rankcutoff <- max(heat.temp$rank)-max(heat.temp$rank)*cutbottomprcntile
    heat.temp <- subset(heat.temp, heat.temp$rank < heat.temp.rankcutoff)
    groupname <<- x
    write.table(heat.temp, file = paste("heatmap.",groupname,".cut",cutbottomprcntile,".bed",sep=""),sep = "\t",col.names = F,row.names = F)
    write.table(heat.temp, file = paste("heatmap.allgroupscut",cutbottomprcntile,".bed",sep=""),sep = "\t", append = T,col.names = F,row.names = F)
  })
  
  
#### 3. Read in genesets, remove by heatmap defined cut off, write out results
  genesets <- list.files(pattern = "_tss.bed", full.names=F, recursive=F)
  
  heatcuts = read.table(paste("heatmap.allgroupscut",cutbottomprcntile,".bed",sep=""))
    colnames(heatcuts) = c("chrom","start","end","name","score","strand","thickStart","thickEnd","itemRGB","blockCount","blockSizes","blockStart","deepTools_group","rank")
  
  lapply(genesets, function(x) {
    #read in geneset bed
    geneset.bed <- read.table(x,sep="\t", header=F)
    geneset.bed$name = paste(geneset.bed$V1,":",geneset.bed$V2,"-",geneset.bed$V3,sep="")
    geneset.bed.cut = subset(geneset.bed, geneset.bed$name %in% heatcuts$name)
    geneset.name = x; geneset.cut.name = gsub(pattern = ".txt_tss.bed",replacement = ".heatcut.txt",x = geneset.name)
    write.table(geneset.bed.cut$V4,file = geneset.cut.name, row.names = F, col.names = F, quote = F)
  })
  
    
    
    
    