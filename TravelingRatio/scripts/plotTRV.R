###########################
# Traveling ratio analysis
# Berkley Gryder
# gryderart@gmail.com
###########################

##### 1. load TR files
  setwd("K:/projects/ChIP_seq/projects/")
  project.folder = "KDM5A/TravelingRatio/"
  PromLength = 0.770 #in kilobases
  TSSRLength = 0.330 #in kilobases

  #### Gene set to split TSS's by 
  Geneset = read.table(file = paste(project.folder,"KDM5AandMYC_JQKD82notdown.genelist.txt",sep=""), sep="\t", header = F)
  geneset.name = "GenesNotDown_JQKD82"

  ### TR files output from scripts/calculateTR.sh
  TRV1.file = "Sample_MM1S_DMSO_Pol2_DFCI_C_JQ.trv"
  TRV1.name = "MM1S_Pol2_DMSO"
  TRV1 = read.table(paste(project.folder,TRV1.file,sep=""), header=F, sep="\t", fill=T)
    colnames(TRV1) = c("Chromosome","TSSR_Start", "TSSR_End", "GeneBody_Start", "GeneBody_End", "Strand", "Symbol", "IsOverlapped", "Traveling_Ratio", "reads.Promoter", "reads.TSSR", "reads.GeneBody","reads.TESR")
    TRV1$iRatio = ((TRV1$reads.Promoter+1)/PromLength)/((TRV1$reads.TSSR+1)/TSSRLength)
    TRV1.pBound = subset(TRV1, TRV1$reads.TSSR > quantile(TRV1$reads.TSSR, 0.6))
    TRV1.Geneset = subset(TRV1, TRV1$Symbol %in% Geneset$V1)
        library(ggplot2)
        ggplot(TRV1, aes(x=log10(reads.TSSR+1)))+geom_density()+geom_density(data=TRV1.pBound, color = "red")
    
    colnames(TRV1) = paste(colnames(TRV1),TRV1.name,sep="_")

  #load condition 2
  TRV2.file = "Sample_MM1S_JQKD82_Pol2_DFCI_C_JQ.trv"
  TRV2.name = "MM1S_Pol2_JQKD82"
  TRV2 = read.table(paste(project.folder,TRV2.file,sep=""), header=F, sep="\t", fill=T)
    colnames(TRV2) = c("Chromosome","TSSR_Start", "TSSR_End", "GeneBody_Start", "GeneBody_End", "Strand", "Symbol", "IsOverlapped", "Traveling_Ratio", "reads.Promoter", "reads.TSSR", "reads.GeneBody", "reads.TESR")
    TRV2$iRatio = ((TRV2$reads.Promoter+1)/PromLength)/((TRV2$reads.TSSR+1)/TSSRLength)
    TRV2.pBound = subset(TRV2, TRV1$reads.TSSR > quantile(TRV1$reads.TSSR, 0.6))
    TRV2.Geneset = subset(TRV2, TRV2$Symbol %in% Geneset$V1)
    colnames(TRV2) = paste(colnames(TRV2),TRV2.name,sep="_")

  TRV.all = cbind(TRV1,TRV2[,c(9,10,11,12,13,14)])
  TRV.geneset.all = cbind(TRV1.Geneset,TRV2.Geneset[,c(9,10,11,12,13,14)])
  

#comparison plots
  #TRV.all$Pol2_D48_GeneDensity = TRV.all$reads.TSSR_Pol2_D48/TSSRLength+TRV.all$reads.Promoter_Pol2_D48/PromLength+TRV.all$reads.GeneBody_Pol2_D48/((TRV.all$GeneBody_End_Pol2_D48-TRV.all$GeneBody_Start_Pol2_D48)/1000)
  #TRV.all = subset(TRV.all, TRV.all$reads.Promoter_scr_Pol2 > quantile(TRV.all$reads.Promoter_scr_Pol2, 0.6))
  library(LSD)
  comparisonplot(log2(TRV.all$reads.Promoter_MM1S_Pol2_DMSO/PromLength), log2(TRV.all$reads.Promoter_MM1S_Pol2_JQKD82/PromLength), colpal="ylgnbu",... = abline(a=0, b=1),xlim = c(0,12),ylim = c(0,12))
  comparisonplot(log2(TRV.all$reads.Promoter_MM1S_Pol2_DMSO/PromLength), log2(TRV.all$reads.Promoter_MM1S_Pol2_JQKD82/PromLength), colpal="ylgnbu",... = abline(a=0, b=1),xlim = c(0,12),ylim = c(0,12))
  comparisonplot(log2(TRV.all$Traveling_Ratio_MM1S_Pol2_DMSO), log2(TRV.all$Traveling_Ratio_MM1S_Pol2_JQKD82), colpal="ylgnbu")
  
  ###write TRV files, filtered for Pol2 bound in the promoter####
  write.table(TRV.all,paste(project.folder,TRV1.name,"_all.txt",sep=""),sep="\t",row.names=F, quote=F, col.names=T)
  write.table(TRV1.pBound,paste(project.folder,TRV1.name,"_pBound.trv",sep=""),sep="\t",row.names=F, quote=F, col.names=F)
  write.table(TRV2.pBound,paste(project.folder,TRV2.name,"_pBound.trv",sep=""),sep="\t",row.names=F, quote=F, col.names=F)
  write.table(TRV.geneset.all,paste(project.folder,TRV1.name,geneset.name,"_all.txt",sep=""),sep="\t",row.names=F, quote=F, col.names=T)
    
  ######################calculate TR cumulative percentages (using scripts/cumulateValue.pl):
  #find /xxx -name '*.trv' -exec /data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts/cumulateValue.pl -i {} -o {}.out \;
  #where /xxx is the folder name of your trv files
  
  #load .trv.out files
  TRV1.file.out = read.table(paste(project.folder,TRV1.name,"_pBound.trv.out",sep=""), header=F, sep="\t", fill=T)
  TRV1.file.out$condition = TRV1.name
  
  TRV2.file.out = read.table(paste(project.folder,TRV2.name,"_pBound.trv.out",sep=""), header=F, sep="\t", fill=T)
  TRV2.file.out$condition = TRV2.name
  
  TRV.out.all = rbind(TRV1.file.out,TRV2.file.out)
  colnames(TRV.out.all) = c("Traveling_Ratio","GeneCount","Percentage_of_Genes","condition")
  
  ##### Treated v untreated plots

  library(ggplot2)
  ggplot(TRV.out.all, aes(x=log10(Traveling_Ratio),y=Percentage_of_Genes,color=condition))+geom_line()+theme_bw()
  
### 3. Build TRV out and plot again on a gene set of interest
  
  ###write subset of TRV files, filtered for Pol2 bound in the promoter and for gene set of interest
  
  write.table(subset(TRV1.pBound, TRV1.pBound$Symbol %in% Geneset$V1),paste(project.folder,TRV1.name,"_pBound.",geneset.name,".trv",sep=""),sep="\t",row.names=F, quote=F, col.names=F)
  write.table(subset(TRV2.pBound, TRV2.pBound$Symbol %in% Geneset$V1),paste(project.folder,TRV2.name,"_pBound.",geneset.name,".trv",sep=""),sep="\t",row.names=F, quote=F, col.names=F)

  ######################calculate TR cumulative percentages (using scripts/cumulateValue.pl):
  #find /xxx -name '*.trv' -exec /data/khanlab/projects/ChIP_seq/data_by_file_type/test/scripts/cumulateValue.pl -i {} -o {}.out \;
  #where /xxx is the folder name of your trv files
  
  #load .trv.out files
  TRV1.geneset.out = read.table(paste(project.folder,TRV1.name,"_pBound.",geneset.name,".trv.out",sep=""), header=F, sep="\t", fill=T)
  TRV1.geneset.out$condition = paste(TRV1.name,geneset.name)
  
  TRV2.geneset.out = read.table(paste(project.folder,TRV2.name,"_pBound.",geneset.name,".trv.out",sep=""), header=F, sep="\t", fill=T)
  TRV2.geneset.out$condition = paste(TRV2.name,geneset.name)
  
  TRV.geneset.out = rbind(TRV1.geneset.out,TRV2.geneset.out)
  colnames(TRV.geneset.out) = c("Traveling_Ratio","GeneCount","Percentage_of_Genes","condition")
  
  TRV.all.subsets.out = rbind(TRV.out.all,TRV.geneset.out)
  ##### Treated v untreated plots
  
  library(ggplot2)
  ggplot(TRV.all.subsets.out, aes(x=log10(Traveling_Ratio),y=Percentage_of_Genes,color=condition))+geom_line()+theme_bw()
  
  t.test(TRV1.geneset.out[,1], TRV2.geneset.out[,1], paired = F)
  t.test(TRV1.pBound$Traveling_Ratio, TRV1.Geneset$Traveling_Ratio, paired = F)
  t.test(TRV1.pBound$Traveling_Ratio, TRV2.pBound$Traveling_Ratio, paired = F)
  t.test(TRV1.Geneset$Traveling_Ratio, TRV2.Geneset$Traveling_Ratio, paired = F)
