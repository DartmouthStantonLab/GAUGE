
## install packages first 
#install.packages("stringdist")
#install.packages("dynamicTreeCut")
#install.packages("ape")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("GEOquery")
#BiocManager::install("Biobase")

library(stringdist)
library(dynamicTreeCut)
library(GEOquery)
library(Biobase)
library(ape)

#three random accession numbers are selected.
#they are all Pseudomonas studies using microarray technique (GPL84)
#you can feed what ever GSE# corresponding to microarray studies into object "expname"
expname<-c("GSE10030","GSE10065","GSE28429")

pdf("allPAdend_microarray.pdf") #this will create a pdf file of dendrograms for validation
report<-lapply(expname, function(y){

  print(y)
  print(which(expname ==y))
  
  #use getGEO to download data matrix, this would need internet connection
  gset <- getGEO(y, GSEMatrix =TRUE, AnnotGPL=TRUE)
  #if there are more than one microarray platform used for this study, select "GPL84"
  if (length(gset) > 1) idx <- grep("GPL84", attr(gset, "names")) else idx <- 1
  gset <- gset[[idx]]
  
  
  # make proper column names to match toptable 
  fvarLabels(gset) <- make.names(fvarLabels(gset))
  
  # log2 transform the expression data
  ex <- exprs(gset)
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  if (LogC) { ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) }
  
  
  
  #put transformed expression data into expMat
  expMat<-exprs(gset)
  colnames(expMat)<-as.character(gset$title)
  expMat<-na.omit(expMat)
  
  if(ncol(expMat)==1){return("NULL_only1Sample")}
  
  namesForDendrogram<-colnames(expMat)
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
  
  #assign sample groups 
  sampleGroups <- 
    cutreeDynamic(hclust(AllSampleNamesDist),minClusterSize = 1, method = "hybrid", deepSplit=4,distM = as.matrix(AllSampleNamesDist))
  x<-table(sampleGroups) 
  d = as.data.frame(x)
  
  if(length(d$sampleGroups) == 1){return("NULL_onlyonegroup")}
  if(min(d$Freq)==1){return("NULL_oneSamplepergroup")}
  
  
  #### pull out experiments have same dimension of stringdist matrix and data distance matrix
  
  
  if (dim(as.matrix(AllSampleNamesDist))[1] != dim(as.matrix(dist(t(expMat), upper = TRUE, diag = TRUE)))[1]) {return("NULL_dimdifferent")}
  
  
  #calculate Matel test p-value
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
dev.off() # close the graphic device to get the pdf file of dendrograms

