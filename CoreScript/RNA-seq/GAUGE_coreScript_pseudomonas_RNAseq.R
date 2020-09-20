## install required packages first 
#install.packages("stringdist")
#install.packages("dynamicTreeCut")
#install.packages("ape")
#install.packages("parallel")
#install.packages("data.table")
#install.packages("rjson")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("edgeR")

#load required packages 
library(stringdist)
library(dynamicTreeCut)
library(ape)
library(parallel)
library(data.table)
library(edgeR)
library(rjson)

#in this demonstration of GAUGE core script, 
#analyzing RNAseq datasets processed by refine.bio,
#I randomly chose three experiments with accession#: SRP228531, SRP055410, SRP090296
#The selected studies were processed by refine.bio with parameters: 
#a). aggregate by experiment. b).no transformation 
#c). Skip quantile normalization

#the output file from refine.bio was downloaded and unzipped
#All files, even those in the sub-directories, were placed at this r project directory 

#aggregated_metadata.json was a file created and provided by refine.bio
#for detail metadata structure of this file, please visit:
#http://docs.refine.bio/en/latest/main_text.html#the-download-folder-structure-for-data-aggregated-by-species

#clean and create the metadata I want from the aggregated_metadata.json
 
metadata_list <- fromJSON(file = "aggregated_metadata.json")

metadata_samples<-sapply(metadata_list$experiments,function(x){
  samples<-x$sample_accession_codes
})
self_metadata_df<-stack(metadata_samples)
colnames(self_metadata_df)<-c("sample_accession","experiment_accession")

metadata_sample_tech<-lapply(metadata_list$experiments,function(x){
  sample_accession<-x$sample_accession_codes
  technology<-x$technology
  technology<-rep(technology,length(sample_accession))
  data.frame(sample_accession,technology)
})


metadata_sample_tech<-rbindlist(metadata_sample_tech)
dim(self_metadata_df)    ##30*2
dim(metadata_sample_tech) ##30*4
length(unique(self_metadata_df$sample_accession)) ##No samples are used in multiple experiments!!!
length(unique(metadata_sample_tech$sample_accession))
self_metadata_df<-merge(self_metadata_df,metadata_sample_tech,by="sample_accession")

## now, extract refinebio_title from each object of metadata_list$samples

metadata_samples_title<-sapply(metadata_list$samples,function(x){
  title<-x$refinebio_title
})
self_metadata_sample_title<-stack(metadata_samples_title)
colnames(self_metadata_sample_title)<-c("refinebio_title","sample_accession")
self_metadata_df<-merge(self_metadata_df,self_metadata_sample_title,by="sample_accession")
dim(self_metadata_df) ##30*4
str(self_metadata_df) ##this metadata dataframe has all information that GAUGE needs, 
                      ##besides the expression data
length(unique(self_metadata_df$sample_accession))


#expGSE is a vecotor containing all unique experiment accession
expGSE<-unique(self_metadata_df$experiment_accession)
length(expGSE)

pdf("allPAdend_RNAseq.pdf") #this will create a pdf file of dendrograms for validation
library(parallel)
## the following lapply function will read in each unique experiment accession#.tsv,
## which stores the expression data, perform the GAUGE analysis to caculate the Mantel test p-value
## and assign sample groups, if available
report<-mclapply(expGSE,mc.cores = 1,mc.preschedule = FALSE,function(y){
  
  expMat<-fread(paste(y,".tsv", sep = ""))#read in the expression matrix
  expMat<-expMat[,-1] #remove the first EMSMBL ID column
  
  
  if(ncol(expMat)==1){return("NULL_only1Sample")}
  if(ncol(expMat)!=length(as.character( 
    self_metadata_df$refinebio_title[match(colnames(expMat),self_metadata_df$sample_accession)] 
  ))) {return("NULL_#ofGSM!=#ofSampleName")}
  
 
  colnames(expMat)<-self_metadata_df$refinebio_title[match(colnames(expMat),self_metadata_df$sample_accession)]
  x=colnames(expMat)  
  # Now make a distance matrix
  myMat <- matrix(nrow = length(x),
                  ncol = length(x),
                  data = 1)
  
  for (i in 1:dim(myMat)[1]) {
    for (j in 1:dim(myMat)[2]) {
      myMat[i, j] <- stringdist(x[i], x[j])
    }
  }
  
  AllSampleNamesDist <-as.dist(myMat)
  
  sampleGroups <- 
    cutreeDynamic(hclust(AllSampleNamesDist),minClusterSize = 1, method = "hybrid", deepSplit=4,distM = as.matrix(AllSampleNamesDist))
  x<-table(sampleGroups) 
  d = as.data.frame(x)
  
  if(length(d$sampleGroups) == 1){return("NULL_onlyonegroup")}
  if(min(d$Freq)==1){return("NULL_oneSamplepergroup")}
  
  
  #### pull out experiments have same dimension of stringdist matrix and data distance matrix
  
  
  if (dim(as.matrix(AllSampleNamesDist))[1] != dim(as.matrix(dist(t(expMat), upper = TRUE, diag = TRUE)))[1]) {return("NULL_dimdifferent")}
  
  ###use edgeR to filter low abundant count and normalization
  ###since some experiments have same refinebio_sample title for multiple samples, which cause error for edgeR analysis, 
  ###set colume names to NULL before edgeR analysis.
  namesForDendrogram<-colnames(expMat)
  
  colnames(expMat)<-NULL
  DGE1 <- DGEList(counts = expMat)
  keep <- filterByExpr(DGE1)
  DGE1<- DGE1[keep, , keep.lib.sizes=FALSE]
  DGE1 <- calcNormFactors(DGE1) # adding information to DGEList object!
 
  #cluster analysis
  CPM1 <- cpm(DGE1, normalized.lib.sizes = TRUE, log = TRUE)
  expMat <- as.data.frame(CPM1)
  
  
  MantelPval <- 
    mantel.test(as.matrix(AllSampleNamesDist),
                as.matrix(dist(t(expMat), upper = TRUE, diag = TRUE)))$p
  
  ##add colnames and rownames into two distant matrix and then draw dengrogram
  namesForDendrogram<-paste(namesForDendrogram,sampleGroups,sep = "/")
  
  rownames(myMat)<-namesForDendrogram
  colnames(myMat)<-namesForDendrogram
  
  colnames(expMat)<-namesForDendrogram
  
  
  myMatPlotMain<-paste(y,";string distance",";pValue=",MantelPval,sep = "")
  expMatPlotMain<-paste(y,";data distance",sep = "")
  plot(hclust(as.dist(myMat)),main = myMatPlotMain,cex=0.3,cex.main=0.7)
  plot(hclust(dist(t(expMat))),main = expMatPlotMain,cex=0.3,cex.main=0.7)
  return(list(MantelPval,sampleGroups))
  
})

names(report)<-expGSE
report
str(report)
# the report object is a list, which has three elements for three unique experiment accession.
# Two sub-elements can be found for each experiment. The first one is the Mantel test p value
# The second is a vector storing the sample group assignment.
dev.off() # close the graphic machine to get the pdf file of dendrograms
