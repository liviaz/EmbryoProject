# Process output from HTSeq-count and turn it into a count table

source("http://bioconductor.org/biocLite.R")
library(Rsamtools)
library(DESeq)
library(DESeq2)
library(edgeR)
library(rtracklayer)
library(GenomicRanges)
library(ape)
library(gplots)
library(ggplot2)
library(RColorBrewer)

if (Sys.info()["sysname"] == "Linux") {
  source("/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/Scripts/R/UsefulFunctions.R")
  baseDataDirectory <- "/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/MouseHiSeq/HTSeq_Count_Out"
} else {
  source("C:/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/Scripts/R/UsefulFunctions.R")
  baseDataDirectory <- "C:/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/MouseHiSeq/HTSeq_Count_Out"
}

##################################################################################################
# Prepare count table (output from HTSeq-count)                                                  #
##################################################################################################

fileNameList <- c("bad1", "bad2", "good1", "good2")

# read in count tables for all samples
bad_1_count <- read.table(paste(baseDataDirectory, "/", fileNameList[1], "_paired_count.txt", sep = ""), header = FALSE, row.names=1, col.names=c("name", "bad_1"))
bad_2_count <- read.table(paste(baseDataDirectory, "/", fileNameList[2], "_paired_count.txt", sep = ""), header = FALSE, row.names=1, col.names=c("name", "bad_2"))
good_1_count <- read.table(paste(baseDataDirectory, "/", fileNameList[3], "_paired_count.txt", sep = ""), header = FALSE, row.names=1, col.names=c("name", "good_1"))
good_2_count <- read.table(paste(baseDataDirectory, "/", fileNameList[4], "_paired_count.txt", sep = ""), header = FALSE, row.names=1, col.names=c("name", "good_2"))


# combine tables into one
nRowToUse = (nrow(bad_1_count) - 5)
allCounts <- data.frame(matrix(ncol = 4, nrow = nRowToUse))
names(allCounts) <- c("bad_1", "bad_2", "good_1", "good_2")
row.names(allCounts) <- row.names(bad_1_count)[1:nRowToUse]

allCounts[,"bad_1"] <- bad_1_count[1:nRowToUse,"bad_1"]
allCounts[,"bad_2"] <- bad_2_count[1:nRowToUse,"bad_2"]
allCounts[,"good_1"] <- good_1_count[1:nRowToUse,"good_1"]
allCounts[,"good_2"] <- good_2_count[1:nRowToUse,"good_2"]

exptDesign <- data.frame(row.names = colnames(allCounts), condition = c("bad", "bad", "good", "good"), 
                         libType = c("paired-end", "paired-end", "paired-end", "paired-end"))

condition = exptDesign$condition
