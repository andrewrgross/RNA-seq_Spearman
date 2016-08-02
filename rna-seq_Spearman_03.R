### RNA Seq preliminary analysis, Andrew R Gross, 2016-05-16
### This script is intended to upload normalized expression data and plot a variety of features in order to make general assessments of the data

########################################################################
### Header
########################################################################

library(gplots)
library(RColorBrewer)
library(biomaRt)
library(DESeq2)
library(Hmisc)

#ensembl = useMart(host="www.ensembl.org")
ensembl = useMart(host="www.ensembl.org",dataset="hsapiens_gene_ensembl")
listMarts(host="www.ensembl.org")
ensembl = useMart(host='www.ensembl.org',biomart='ENSEMBL_MART_ENSEMBL')
listDatasets(ensembl)
ensembl = useMart(host='www.ensembl.org',biomart='ENSEMBL_MART_ENSEMBL',dataset="hsapiens_gene_ensembl")
#filters = listFilters(ensembl)
#attributes = listAttributes(ensembl)

########################################################################
### Functions
########################################################################

addMedSD <- function(dataframe) {
  median <- apply(dataframe,1,median)
  sd <- apply(dataframe,1,sd)
  return(data.frame(dataframe,median,sd))  
}

sortByMed <- function(dataframe) {
  order <- order(dataframe$median,decreasing=TRUE)
  return(dataframe[order,])
}

convertIDs <- function(dataframe) {
  ensemblIDs <- c()
  for (rowName in row.names(dataframe)) {
    ensemblID <- strsplit(rowName,"\\.")[[1]][1]
    ensemblIDs <- c(ensemblIDs,ensemblID)
  }
  row.names(dataframe)<-ensemblIDs
  return(dataframe)
}
########################################################################
### Import Data
########################################################################

# Normalized
TPMdata <- read.csv("z://Data/RNAseq HT neurons and tissue/Andrews_files/20160317_Sareen_rerun-29299281.tpm.csv", row.names=1)
metaData <- read.csv("z://Data/RNAseq HT neurons and tissue/2nd rerun/samples.csv", row.names=1)

# Import references

references <- read.table("c://Users/grossar/Bioinform/DATA/rna_seq/Reference_transcriptomes/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct",sep="\t",header=TRUE,row.names=1)

########################################################################
### Formating
########################################################################

sampleNames <- c("iHT_03iCTR","iHT_90iOBS","iHT_77iOBS","iHT_02iOBS","iMN_87iCTR","iMN_201iCTR","aHT_1662_S6","aHT_1838_S7","aHT_1843_S8","aHT_2266_S9","iHT_02iCTR_S13","aHT_2884_S11","21-Reference","iHT_87iCTR","iHT_201iCTR","iHT_25iCTR_S16","iHT_688iCTR_","iHT_80iCTR","iHT_74iOBS","iHT_03iOBS")
names(TPMdata) <- sampleNames
TPMdata <- convertIDs(TPMdata)

sampleTranscriptomeList <- list()
for (sampleNumber in 1:length(TPMdata)) {
  currentDataframe <- TPMdata[sampleNumber]
  order <- order(currentDataframe[1], decreasing=TRUE)
  sampleTranscriptomeList[[sampleNumber]] <- data.frame(currentDataframe[order,],row.names=row.names(currentDataframe)[order])
}
names(sampleTranscriptomeList) <- sampleNames

#references <- references[2:48]
references <- convertIDs(references)
referenceTranscriptomeList <- list()
for (tissueNumber in 2:length(references)) {
  currentDataframe <- references[tissueNumber]
  order <- order(currentDataframe[1],decreasing = TRUE)
  referenceTranscriptomeList[[(tissueNumber-1)]] <- data.frame(currentDataframe[order,],row.names=row.names(currentDataframe)[order])
}
names(referenceTranscriptomeList) <- names(references[2:length(references)])

########################################################################
### Add samples to full data list
########################################################################

transcriptomeList <- append(referenceTranscriptomeList,sampleTranscriptomeList)
transcriptomeList <- referenceTranscriptomeList
transcriptomeList <- sampleTranscriptomeList

########################################################################
### Subsample
########################################################################

specifiedEnd <- 15000
for (elementNumber in 1:length(transcriptomeList)) {
  currentDatafame <- transcriptomeList[[elementNumber]]
  #selectedRows <- currentDatafame[,1] > 1 
  selectedRows <- 1:specifiedEnd
  transcriptomeList[[elementNumber]] <- data.frame(currentDatafame[1][selectedRows,],row.names=row.names(currentDatafame)[selectedRows])
}

info <- print(str(transcriptomeList))

########################################################################
### Generate list of unique ids
########################################################################

allIDs <- c()
for (df in transcriptomeList) {
  allIDs <- c(allIDs,row.names(df))
}
uniqueIDs <- unique(allIDs)
length(uniqueIDs)

### Generate a dataframe with each of the unique IDs

transcriptsDF <- data.frame(uniqueIDs)
for (df in transcriptomeList) {
  unusedIDs <- setdiff(uniqueIDs,row.names(df))
  newDF <- data.frame(c(df[,1],rep(0,length(unusedIDs))))
  row.names(newDF) <- c(row.names(df),unusedIDs)
  newDF <- newDF[order(row.names(newDF)),]
  transcriptsDF[length(transcriptsDF)+1] <- newDF
}

row.names(transcriptsDF) <- transcriptsDF[,1]
transcriptsDF <- transcriptsDF[2:length(transcriptsDF)]
names(transcriptsDF) <- names(transcriptomeList)

transcriptsMatrix <- as.matrix(transcriptsDF)

########################################################################
### Compare, v3
########################################################################

comparisonMatrix <- rcorr(transcriptsMatrix, type="spearman")[[1]]
comparisonMatrix <- round(comparisonMatrix*100,0)

######################################################################################
### Plot heatmap
######################################################################################

nColors=99
start = 50
end = 100
my_palette <- colorRampPalette(c("yellow","white","blue"))(n = nColors)
my_palette <- colorRampPalette(c("red","black","green"))(n = nColors)
my_palette <- colorRampPalette(c("black","green"))(n = nColors)

col_breaks <- seq(start,end,(end-start)/(nColors))

heatmap.2(comparisonMatrix,
          main = paste("Spearman comparison",specifiedEnd), # heat map title
          cellnote = comparisonMatrix,
          notecol = "black",
          density.info="none",  # turns off density plot inside color legend
          notecex=0.6,
          trace="none",         # turns off trace lines inside the heat map
          col=my_palette,       # use on color palette defined earlier 
          distfun=dist,
          #dendrogram="column",     # only draw a row dendrogram
          #Rowv="NA",
          margins =c(12,9),     # widens margins around plot
          breaks=col_breaks    # enable color transition at specified limits
          
          #cellnote = comparisonMatrix,  # same data set for cell labels
          #notecol="black",      # change font color of cell labels to black
          #Colv="NA"             # turn off column clustering

          )


