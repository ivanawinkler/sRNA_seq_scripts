#####################################################################
#libraries
#####################################################################

library(jsonlite)
library("DESeq2")
library(biomaRt)
library(BiocParallel)
library(ggplot2)
library(RColorBrewer)
library("gplots")

register(MulticoreParam(6))

#---------------------------------------------------------------------------------------------------------------------------------
#Function that reads in files from STAR in one dataframe and metadata-takes as arguments: mypath (path to folder containing 
#read count files); mypath_metadata (path to file containing metadata-dataframe-columns: sample names, conditions and type(SR or PE))
#----------------------------------------------------------------------------------------------------------------------------------
setwd("~/miRNA/sRNA-seq/TCGA")

file_importeR=function(mypath){
  filenames=list.files(mypath, full.names=TRUE) #contains names of files in the folder
  mydataframefile=read.table(filenames[1],header=T,sep="\t",row.names = NULL,stringsAsFactors = F) #reads in the first file
  mydataframefile=mydataframefile[,1:2]
  colnames(mydataframefile)=c("V1",filenames[1])
  for (x in seq(from=2, to=length(filenames))){
    mydataframe=read.table(filenames[x],header=T,sep="\t",row.names = NULL,stringsAsFactors = F) #reads in other files
    mydataframe=mydataframe[,1:2]
    colnames(mydataframe)=c("V1",filenames[x])
    mydataframefile=merge(mydataframefile,mydataframe,by.x="V1",by.y="V1") #merges the other files
  }
  rownames(mydataframefile)=mydataframefile[,1]
  mydataframefile=mydataframefile[,2:length(mydataframefile)]
  assign("countData",mydataframefile,envir = .GlobalEnv) #returns the dataframe to global environment
}

file_importeR("./raw_counts")

#remove name of the folder- DON'T forget (you won't be able to match) 
colnames(countData)=gsub("./raw_counts/","",colnames(countData))

#generating colData
setwd("~/miRNA/sRNA-seq/TCGA/metadata")
metadata=read.table("gdc_sample_sheet.2018-11-09.tsv",header=T,sep = "\t",check.names = F)
metadata_stripped=metadata[,c(2,7)]
colData=as.data.frame(metadata_stripped[1,])
if (substr(as.character(colData[1,2]),14,15)=="11"){colData[1,3]="control"} else {colData[1,3]="tumor"}
for (i in seq(from=2,to=dim(metadata)[1])){
  colData[i,1]=metadata_stripped[i,1]
  colData[i,2]=metadata_stripped[i,2]
  if (substr(as.character(colData[i,2]),14,15)=="11"){colData[i,3]="control"} else {colData[i,3]="tumor"}}

#changing names in countdata table to barcodes
names_list=setNames(as.list(colData[,1]),colData[,2])
for (i in seq(from=1, to=length(countData))){
  for (j in seq(from=1, to=length(names_list))){
    if (colnames(countData)[i]==names_list[[j]]){colnames(countData)[i]=names(names_list[j])}
  }
}
countData=countData[,order(names(countData))]
colData=colData[order(colData$`Sample ID`),]
colnames(colData)=c("file_name","barcode","condition")
colData=colData[,2:3]
rownames(colData)=make.names(colData$barcode,unique = T)
colnames(countData)=rownames(colData)
countData[]=lapply(countData, function (x) as.numeric(x))

#---------------------------------------------------------------------------------------------------------------------------------
#DEG function - takes as arguments: counData (dataframe containing gene counts with sample_names as columns and genes as rownames); 
#               colData (metadata-dataframe-columns: sample names, conditions and type(SR or PE)); condition_list (list of paired 
#               conditions which are compared); cook (logical: T (default)-cook distance filtering turned on)
#---------------------------------------------------------------------------------------------------------------------------------
multiple_deseq2_analyseR=function(countData,colData,cook=T){
  
  dds=DESeqDataSetFromMatrix(countData=countData,colData=colData,design=~condition)
  dds$conditions=factor(dds$condition,levels = c("control","tumor"))
  dds=dds[rowSums(counts(dds))>1,]
  #differential expression analysis
  dds=DESeq(dds, parallel = TRUE)
  res=results(dds,cooksCutoff = cook, parallel = TRUE)
  return(res) #returns the list to global environment #returns the list to global environment
}

results_tcga=multiple_deseq2_analyseR(countData = countData, colData = colData)
results_tcga_df=as.data.frame(results_tcga)
dds=DESeqDataSetFromMatrix(countData=countData,colData=colData,design=~condition)
dds$conditions=factor(dds$condition,levels = c("control","tumor"))
dds=dds[rowSums(counts(dds))>1,]
#differential expression analysis
dds=DESeq(dds, parallel = TRUE)
rld=varianceStabilizingTransformation(dds,blind=F)
rld_df=assay(rld)
rld_df=as.data.frame(rld_df)
colnames(rld_df)=gsub("\\.","-",colnames(rld_df))
rld_df$miRNA=row.names(rld_df)
#colnames(results_tcga_df)[1]="miRNA"
write.table(rld_df,"TCGA_miRNA_rld",quote = F,row.names = F,sep=",")

#--------------------------------------------------------------------------------------------------------------------------------
#generating normalized gene counts for the tumor samples (normalization to control)
#-------------------------------------------------------------------------------------------------------------------------------
survival_tumor=rld_df[,grep("01A",colnames(rld_df))]
survival_tumor_p=rld_df[,grep("02A",colnames(rld_df))]
survival_tumor=cbind(survival_tumor,survival_tumor_p)
#--------------------------------------------------------------------------------------------------------------------------------
#generating normalized gene counts for the control samples & averageing them
#-------------------------------------------------------------------------------------------------------------------------------
survival_control=rld_df[,grep("11A",colnames(rld_df))]
survival_control_average=data.frame(Means=rowMeans(survival_control[,-1]))

#--------------------------------------------------------------------------------------------------------------------------------
#generating normalized gene counts for the tumor samples normalise to the control
#-------------------------------------------------------------------------------------------------------------------------------

survival_tumor_av=as.data.frame(survival_tumor[1]-survival_control_average$Means)
for (i in seq(from=1,to=dim(survival_tumor)[2])){
  survival_tumor_av[i]=survival_tumor[i]-survival_control_average$Means}
colnames(survival_tumor_av)=colnames(survival_tumor)
rownames(survival_tumor_av)=row.names(rld_df)


#-------------------------------------------------------------------------------------------------------
#genearting table with gene expression which is to be used in heatmap construction
#------------------------------------------------------------------------------------------------------

#genes=c("hsa-let-7a-5p","hsa-let-7c-5p","hsa-let-7g-5p","hsa-miR-29c-3p","hsa-miR-335-3p","hsa-miR-338-3p","hsa-miR-30d-5p","hsa-miR-30e-5p")
genes=c("hsa-let-7a-1","hsa-let-7a-2","hsa-let-7a-3","hsa-let-7c","hsa-let-7g","hsa-mir-29c","hsa-mir-335","hsa-mir-30e","hsa-mir-338","hsa-mir-30d")
genes_heatmap=data.frame()
genes_heatmap=survival_tumor_av[grep(genes[1],row.names(survival_tumor_av)),]
for(i in seq(from=2,to=length(genes))){
  genes_heatmap[i,]=survival_tumor_av[grep(genes[i],row.names(survival_tumor_av)),]}

#---------------------------------------------------------------------------------
#genearting heatmaps
#---------------------------------------------------------------------------------
genes_heatmap=as.matrix(genes_heatmap)
genes_heatmap=apply(genes_heatmap,c(1,2),as.numeric)
colors <- colorRampPalette(c("#08306B","#377EB8","white","#E41A1C","#67000D"))(150)
heatmap.2(genes_heatmap,trace="none",col=colors,labCol = F,offsetRow = 0.00001,Rowv = F,margins = c(5,8))

#---------------------------------------------------------------------------------
#generating waterfall graphs
#---------------------------------------------------------------------------------
genes_heatmap_individual=as.data.frame(genes_heatmap[4,])
#genes_heatmap_individual=as.data.frame(t(genes_heatmap_individual))
genes_heatmap_individual[,2]=rownames(genes_heatmap_individual)
for (i in seq(from=1,to=dim(genes_heatmap_individual)[1])){
  if (genes_heatmap_individual[i,1]<(-0.6)){genes_heatmap_individual[i,3]="low"}
  if (genes_heatmap_individual[i,1]>(0.6)){genes_heatmap_individual[i,3]="high"}
  if ((genes_heatmap_individual[i,1]>(-0.6))&(genes_heatmap_individual[i,1]<(0.6))){genes_heatmap_individual[i,3]="no_change"}
}
colnames(genes_heatmap_individual)=c("Tumor_expression","Tumors","Expression")
ggplot(data=genes_heatmap_individual,mapping=aes(x=reorder(Tumors,Tumor_expression),y=Tumor_expression,fill=Expression))+
  geom_bar(stat="identity",position="identity",colour="white",size=0.4)+labs(x="")+ theme(axis.title.x=element_blank(),
                                                                                           axis.text.x=element_blank(),axis.ticks.x=element_blank())+scale_fill_manual(values=c("#E41A1C","#377EB8","grey"),guide=F)+ 
  ggtitle(rownames(genes_heatmap)[4])+ theme(axis.text=element_text(size=16), axis.title=element_text(size=16))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +theme(axis.text=element_text(size=20), axis.title=element_text(size=24))



#----------------------------------------------------------------------------------------
#making dot plots
#----------------------------------------------------------------------------------------
rld_df=ensembleR(rld_df)
rld_df[,1]=rld_df$external_gene_name
rld_df=rld_df[,1:425]
survival_control_gene=rld_df[grep(genes[7],rld_df[,1]),]
survival_control_gene=as.data.frame(survival_control_gene[1,])
survival_control_gene_tum=survival_control_gene[,grep("01A",colnames(survival_control_gene))]
survival_control_gene_tum_2=survival_control_gene[,grep("02A",colnames(survival_control_gene))]
survival_control_gene_tum=cbind(survival_control_gene_tum,survival_control_gene_tum_2)
survival_control_gene_con=survival_control_gene[,grep("11A",colnames(survival_control_gene))]    
survival_control_gene_tum=t(survival_control_gene_tum)
survival_control_gene_tum=as.data.frame(survival_control_gene_tum)
survival_control_gene_tum[,2]="tumor"
survival_control_gene_con=t(survival_control_gene_con)
survival_control_gene_con=as.data.frame(survival_control_gene_con)
survival_control_gene_con[,2]="control"
survival_control_gene=rbind(survival_control_gene_con,survival_control_gene_tum)
colnames(survival_control_gene)=c("Normalized_read_counts","tissue")
ggplot(survival_control_gene,aes(x=tissue,y=Normalized_read_counts,fill=tissue))+geom_boxplot(outlier.color = NA,width=.3)+labs(x="")+
  geom_jitter(position=position_jitter(width=.05, height=0))+scale_fill_manual(values=c("#377EB8","#E41A1C"))+
  scale_color_manual(values=c("#377EB8","#E41A1C"))+ theme(panel.background = element_rect(fill = 'white', colour = "white"))+
  theme( axis.line = element_line(colour = "black", size = 1, linetype = "solid"))+ theme(axis.text=element_text(size=24), axis.title=element_text(size=26))+
  theme(legend.text=element_text(size=22),legend.title = element_text(size=24))
