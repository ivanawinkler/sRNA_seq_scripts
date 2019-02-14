#####################################################################
#libraries
#####################################################################

library(jsonlite)
library("DESeq2")
library(biomaRt)
library(BiocParallel)

register(MulticoreParam(6))

#---------------------------------------------------------------------------------------------------------------------------------
#Function that reads in files from STAR in one dataframe and metadata-takes as arguments: mypath (path to folder containing 
#read count files); mypath_metadata (path to file containing metadata-dataframe-columns: sample names, conditions and type(SR or PE))
#----------------------------------------------------------------------------------------------------------------------------------
setwd("~/sRNA-seq/analysis_30_10_2018")


file_importeR=function(mypath,mypath_metadata){
  colData=read.csv(mypath_metadata,row.names = 1,)#loading countData
  filenames=list.files(mypath, full.names=TRUE) #contains names of files in the folder
  mydataframefile=read.table(filenames[1],header=T,sep="\t",row.names = NULL) #reads in the first file
  mydataframefile=mydataframefile[,c(1,7)]
  colnames(mydataframefile)=c("V1",filenames[1])
  for (x in seq(from=2, to=length(filenames))){
    mydataframe=read.table(filenames[x],header=T,sep="\t",row.names = NULL) #reads in other files
    mydataframe=mydataframe[, c(1,7)]
    colnames(mydataframe)=c("V1",filenames[x])
    mydataframefile=merge(mydataframefile,mydataframe,by.x="V1",by.y="V1") #merges the other files
  }
  colnames(mydataframefile)=sub(paste(mypath,"/",sep=""),"",colnames(mydataframefile))
  rownames(mydataframefile)=mydataframefile[,1]
  countData=as.matrix(mydataframefile[2:NCOL(mydataframefile)])
  colnames(countData)=rownames(colData)
  countData=countData[2:NROW(countData),]
  assign("countData",countData,envir = .GlobalEnv) #returns the dataframe to global environment
  assign("colData",colData,envir = .GlobalEnv) 
}

file_importeR("./DESeq_input/","./colData_m.csv")


#---------------------------------------------------------------------------------------------------------------------------------
#DEG function - takes as arguments: counData (dataframe containing gene counts with sample_names as columns and genes as rownames); 
#               colData (metadata-dataframe-columns: sample names, conditions and type(SR or PE)); condition_list (list of paired 
#               conditions which are compared); cook (logical: T (default)-cook distance filtering turned on)
#---------------------------------------------------------------------------------------------------------------------------------
multiple_deseq2_analyseR=function(countData,colData,condition_list,cook=T){
  library("DESeq2")
  dds=DESeqDataSetFromMatrix(countData=countData,colData=colData,design=~condition)
  
  #filtering out genes with 0 counts
  dds=dds[rowSums(counts(dds))>5,]
  
  #differential expression analysis
  dds=DESeq(dds,parallel = TRUE)
  rld=varianceStabilizingTransformation(dds,blind=F)
  rld_df=assay(rld)
  rld_df=as.data.frame(rld_df)
  res_list=list()
  for (i in seq(from=1, to=length(condition_list))){
    res=results(dds,cooksCutoff=T,contrast=c("condition",condition_list[[i]]))
    colnames(res)=sapply(colnames(res),function(x) paste(x,sapply(condition_list[i], "[", c(1)),sep="_"))
    res_list[[i]]=res
  }
  
  #merging all result tables in one
  res_temp=merge(as.data.frame(res_list[[1]]),as.data.frame(res_list[[2]]),by.x="row.names",by.y="row.names")
  colnames(res_temp)[1]="Gene_id"
  #for (i in seq(from=3, to=length(condition_list))){
   # res_temp=merge(res_temp,as.data.frame(res_list[[i]]),by.x="Gene_id",by.y="row.names")
  #}
  assign("res_temp",res_temp,envir = .GlobalEnv)  #returns the list to global environment #returns the list to global environment
  assign("rld_df",rld_df,envir = .GlobalEnv)
}

condition_list=list(c("nodul","control"),c("tumor","control"))
multiple_deseq2_analyseR(countData = countData, colData = colData,condition_list = condition_list)
results_df=as.data.frame(res_temp)
#-----------------------------------------------------------------------------------------------------------------------------------
#Function that returns dataframe with ensembl_ids, gene description,entrez_id and official gene symbol- it takes as an argument filter_id
#(identifier used for mapping;default ensembl_id) and res- DEG dataframe
#-------------------------------------------------------------------------------------------------------------------------------------
ensembleR=function(res){
  library(biomaRt)
  ensembl=useMart(biomart="ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset="mmusculus_gene_ensembl")
  filters=listFilters(ensembl)
  ensemnl_list=getBM(attributes=c("ensembl_gene_id","description","entrezgene","external_gene_name"),filters="ensembl_gene_id",values = gsub("\\..*","",as.vector(res[,1])),mart=ensembl)
  ensemnl_list=ensemnl_list[!duplicated(ensemnl_list$ensembl_gene_id),]
  res[,1]=gsub("\\..*","",res[,1])
  full_list=merge(res,ensemnl_list,by.x="Gene_id",by.y="ensembl_gene_id")
  assign("full_list",full_list,envir = .GlobalEnv) #returns the list to global environment #returns the list to global environment
}

#ensembleR(results_df)

