#----------------------------------------------
#libraries
#---------------------------------------------
library(jsonlite)
library("DESeq2")
library(biomaRt)
library(BiocParallel)
library(ggplot2)
library(RColorBrewer)
library("gplots")
library(dplyr)


register(MulticoreParam(6))

#importing count files (countData)
#------------------------------------------------------------
setwd("~/RNA_seq/TCGA/lihc/miRNA_reads")

file_importeR=function(mypath){
  filenames=list.files(mypath, full.names=TRUE) #contains names of files in the folder
  mydataframefile=read.table(filenames[1],header=T,sep="\t",row.names = NULL,stringsAsFactors = F) #reads in the first file
  mydataframefile=mydataframefile[,c(6,3)]
  mydataframefile=mydataframefile[!(mydataframefile$miRNA_region=="precursor")&!(mydataframefile$miRNA_region=="unannotated")&!(mydataframefile$miRNA_region=="stemloop"),]
  mydataframefile=as.data.frame(tapply(mydataframefile$read_count, mydataframefile$miRNA_region, FUN=sum))
  mydataframefile$miRNA=rownames(mydataframefile)
  row.names(mydataframefile)=NULL
  colnames(mydataframefile)=c(filenames[1],"V1")
  mydataframefile=mydataframefile[,c(2,1)]
  mydataframefile[,2]=as.numeric(mydataframefile[,2])
  for (x in seq(from=2, to=length(filenames))){
    mydataframe=read.table(filenames[x],header=T,sep="\t",row.names = NULL,stringsAsFactors = F) #reads in the first file
    mydataframe=mydataframe[,c(6,3)]
    mydataframe=mydataframe[!(mydataframe$miRNA_region=="precursor")&!(mydataframe$miRNA_region=="unannotated")&!(mydataframe$miRNA_region=="stemloop"),]
    mydataframe=as.data.frame(tapply(mydataframe$read_count, mydataframe$miRNA_region, FUN=sum))
    mydataframe$miRNA=rownames(mydataframe)
    row.names(mydataframe)=NULL
    colnames(mydataframe)=c(filenames[x],"V1")
    mydataframe=mydataframe[,c(2,1)]
    mydataframe[,2]=as.numeric(mydataframe[,2])
    mydataframefile=merge(mydataframefile,mydataframe,by.x="V1",by.y="V1") #merges the other files
  }
  rownames(mydataframefile)=mydataframefile[,1]
  mydataframefile=mydataframefile[,2:length(mydataframefile)]
  assign("countData",mydataframefile,envir = .GlobalEnv) #returns the dataframe to global environment
}

file_importeR("./raw_reads")
#remove name of the folder- DON'T forget (you won't be able to match) 
colnames(countData)=gsub("./raw_reads/","",colnames(countData))
#---------------------------------------------------------------------------
  
#generating colData
setwd("~/RNA_seq/TCGA/lihc/miRNA_reads/metadata")
metadata=read.table("gdc_sample_sheet.2018-11-09.tsv",header=T,sep = "\t",check.names = F)
metadata_stripped=metadata[,c(2,7)]
colData=as.data.frame(metadata_stripped[1,])
if (substr(as.character(colData[1,2]),14,15)=="11"){colData[1,3]="control"} else {colData[1,3]="tumor"}
for (i in seq(from=2,to=dim(metadata)[1])){
  colData[i,1]=metadata_stripped[i,1]
  colData[i,2]=metadata_stripped[i,2]
  if (substr(as.character(colData[i,2]),14,15)=="11"){colData[i,3]="control"} else {colData[i,3]="tumor"}}
#---------------------------------------------------------------------------------------
  
#changing names in countdata table to barcodes (preparing countData and colData)
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
#----------------------------------------------------------------------

#doing differential expression analysis  
multiple_deseq2_analyseR=function(countData,colData,cook=T){
  
  dds=DESeqDataSetFromMatrix(countData=countData,colData=colData,design=~condition)
  dds$conditions=factor(dds$condition,levels = c("control","tumor"))
  dds=dds[rowSums(counts(dds))>1,]
  #differential expression analysis
  dds=DESeq(dds, parallel = TRUE)
  res=results(dds,cooksCutoff = cook, parallel = TRUE)
  return(res) #returns the list to global environment #returns the list to global environment
}

results_tcga=multiple_deseq2_analyseR(countData = countData, colData = colData)#
results_tcga_df=as.data.frame(results_tcga)
results_tcga_df$miRNA=substr(row.names(results_tcga_df),8,19)
#----------------------------------------------------------------------------------
#mapping miRNA IDs
ensembleR=function(res){
    ensembl=useMart(biomart="ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset="hsapiens_gene_ensembl")
    filters=listFilters(ensembl)
    ensembl_list=getBM(attributes=c("mirbase_id","mirbase_accession"),filters="mirbase_accession",values = res$miRNA, mart=ensembl)
    full_list=merge(res,ensembl_list,by.x="miRNA",by.y="mirbase_accession")
    assign("full_list",full_list,envir = .GlobalEnv) #returns the list to global environment #returns the list to global environment
  }

#ensembleR(results_tcga_df)
miRNA_id=read.delim("mirna_mature.txt",header=F,sep="\t",blank.lines.skip = F)
miRNA_id=miRNA_id[,c(2,4)]
full_list=merge(results_tcga_df,miRNA_id,by.x="miRNA",by.y="V4")
#-----------------------------------------------------------------------------------  
#producing normalized reads
dds=DESeqDataSetFromMatrix(countData=countData,colData=colData,design=~condition)
dds$conditions=factor(dds$condition,levels = c("control","tumor"))
dds=dds[rowSums(counts(dds))>1,]
#differential expression analysis
dds=DESeq(dds, parallel = TRUE)
rld=varianceStabilizingTransformation(dds,blind=F)
rld_df=assay(rld)
rld_df=as.data.frame(rld_df)
colnames(rld_df)=gsub("\\.","-",colnames(rld_df))
rld_df$miRNA=substr(row.names(results_tcga_df),8,19)
rld_df=merge(rld_df,miRNA_id,by.x="miRNA",by.y="V4")


#--------------------------------------------------------------------------------------------------------------------------------
#generating normalized gene counts for the tumor samples
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
rownames(survival_tumor_av)=rld_df$V2


#-------------------------------------------------------------------------------------------------------
#genearting table with gene expression which is to be used in heatmap construction
#------------------------------------------------------------------------------------------------------

#genes=c("hsa-let-7a-5p","hsa-let-7c-5p","hsa-let-7g-5p","hsa-miR-29c-3p","hsa-miR-335-3p","hsa-miR-338-3p","hsa-miR-30d-5p","hsa-miR-30e-5p")
genes=read.table("miRNAs_DLK1_DIO3.txt",header=F)
genes$V1=gsub("mmu-miR","hsa-mir",genes$V1)
genes$V1=gsub("-3p*$","",genes$V1)
genes$V1=gsub("-5p*$","",genes$V1)
common=merge(genes,results_tcga_df,by.x="V1",by.y="mirna")
genes=t(genes)
genes=as.vector(genes)
genes_heatmap=data.frame()
genes_heatmap=survival_tumor_av[grep(genes[6],row.names(survival_tumor_av)),]
for(i in seq(from=2,to=length(genes))){
  genes_heatmap[i,]=survival_tumor_av[grep(genes[i],row.names(survival_tumor_av)),]}



#---------------------------------------------------------------------------------
#generating waterfall graphs
#---------------------------------------------------------------------------------
genes_heatmap_individual=as.data.frame(genes_heatmap[8,])
genes_heatmap_individual=as.data.frame(t(genes_heatmap_individual))
genes_heatmap_individual[,2]=rownames(genes_heatmap_individual)
for (i in seq(from=1,to=dim(genes_heatmap_individual)[1])){
  if (genes_heatmap_individual[i,1]<(-0.6)){genes_heatmap_individual[i,3]="low"}
  if (genes_heatmap_individual[i,1]>(0.6)){genes_heatmap_individual[i,3]="high"}
  if ((genes_heatmap_individual[i,1]>(-0.6))&(genes_heatmap_individual[i,1]<(0.6))){genes_heatmap_individual[i,3]="no_change"}
}
#high=as.numeric(count(genes_heatmap_individual, Expression)[1,2])/372*100
#low=as.numeric(count(genes_heatmap_individual, Expression)[2,2])/372*100
#no_change=as.numeric(count(genes_heatmap_individual, Expression)[3,2])/372*100
colnames(genes_heatmap_individual)=c("log2_fold_change","Tumors","Expression")

#gene=c(paste(round(high,1),"%"), paste(round(no_change,1),"%"),paste(round(low,1),"%"))
#text_height=max(genes_heatmap_individual$log2_fold_change)+(max(genes_heatmap_individual$log2_fold_change)*0.1)
#line_height=max(genes_heatmap_individual$log2_fold_change)-(max(genes_heatmap_individual$log2_fold_change)*0.00005)
col_b=brewer.pal(9,"Blues")
col_r=brewer.pal(9,"Reds")
ggplot(data=genes_heatmap_individual,mapping=aes(x=reorder(Tumors,log2_fold_change),y=log2_fold_change,fill=Expression))+
  geom_bar(stat="identity",position="identity")+labs(x="")+ theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_fill_manual(values=c(col_r[7],col_b[7],"grey"),guide=F)+ 
  ggtitle(rownames(genes_heatmap)[8])+ theme(axis.text=element_text(size=36), axis.title=element_text(size=38))+
  theme(axis.line.y = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(plot.title = element_text(size = 42))+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(y=expression("log"[2]*" fold change"),x=expression(""))
  #geom_segment(aes(x = 0, y = line_height, xend = 75, yend = line_height),size=0.8)+
  #geom_segment(aes(x = 80, y = line_height, xend = 310, yend = line_height),size=0.8)+
  #geom_segment(aes(x = 315, y = line_height, xend = 370, yend = line_height),size=0.8)+
#annotate("text",x=40,y=text_height+(0.005*text_height),label=gene[1],size=8)+
#annotate("text",x=190,y=text_height+(0.005*text_height),label=gene[2],size=8)+
#annotate("text",x=350,y=text_height+(0.005*text_height),label=gene[3],size=8)


