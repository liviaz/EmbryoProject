
## Starts off with DESeq, variance stabilization, etc

rm(list = ls())


if (Sys.info()["sysname"] == "Linux") {
  source("/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/Scripts/R/LoadHTSeqCountData_Human_HiSeq.R")
  # read in RPKM data
  conditionCombined = c(condition, 5,5,5)
  allRPKM = read.table("/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/HumanHiSeq/allHumanEmbryoRPKM.txt")
  yanData = read.table("/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/HumanHiSeq/Yan2013.csv", header = TRUE, row.names = 1)
  combinedRPKMtable = read.table("/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/HumanHiSeq/combinedEmbryoRPKM.txt")
} else {
  source("C:/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/Scripts/R/LoadHTSeqCountData_Human_HiSeq.R")
  # read in RPKM data
  conditionCombined = c(condition, 5,5,5)
  allRPKM = read.table("C:/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/HumanHiSeq/allHumanEmbryoRPKM.txt")
  yanData = read.table("C:/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/HumanHiSeq/Yan2013.csv", header = TRUE, row.names = 1)
  combinedRPKMtable = read.table("C:/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/HumanHiSeq/combinedEmbryoRPKM.txt")
}

# run edgeR
yAll = DGEList(counts = allCounts, group = condition)
yAll = yAll[rowSums(allCounts > 15) > 1,]
yAll = calcNormFactors(yAll)
yAll = estimateDisp(y = yAll, design = model.matrix(~0 + factor(condition)))

# assemble count data set
cdsAll = newCountDataSet(allCounts[rowSums(allCounts > 15) > 1,], condition)
# 2. normalize
cdsAll = estimateSizeFactors(cdsAll)
# 3. estimate dispersions of each gene across all samples
cdsAll = estimateDispersions(cdsAll)
vsdAll = varianceStabilizingTransformation(cdsAll)
# add in decision boundary distances
phenoData(vsdAll)$decDist = c(4.36, 4, 1, 1.1)
phenoData(vsdAll)$exptNum = c(1,2,1,2)
phenoData(vsdAll)$condition = condition
embryosSamePatient = c(1,2,3,4)
exptNum = c(1,2,1,2)

# create model matrix and null matrix
batch = exptNum - 1
patient = c(0, 1, 0, 1)
mod = model.matrix(~as.factor(condition), data = pData(vsdAll))
mod0 = model.matrix(~1, data = pData(vsdAll))

# use ComBat
combat.edata = as.data.frame(ComBat(dat = log2(cpm(yAll)+1), batch = c(0, 1, 0, 1), mod = mod, 
                      par.prior = TRUE, prior.plots = TRUE))
# filter out genes with < 1 cpm average expression
combat.edata = combat.edata[2^(rowMeans(combat.edata) - 1) > 1,]
# calculate pvalues and qvalues
pValuesComBat = f.pvalue(as.matrix(combat.edata), mod, mod0)
qValuesComBat = p.adjust(pValuesComBat, method = "BH")
write.table(as.character(names(qValuesComBat[qValuesComBat < .2])), paste(baseDataDirectory, "/edgeR/ComBat_SVA/qVals_lt_pt2_ComBat.txt", sep = ""), row.names=FALSE, quote=FALSE)
# write.table(as.character(names(qValuesSv[qValuesSv < .1])), paste(baseDataDirectory, "/edgeR/ComBat_SVA/qVals_lt_pt05_SVA.txt", sep = ""), row.names=FALSE, quote=FALSE)
# write.table(combat.edata, paste(baseDataDirectory, "/edgeR/ComBat_SVA/combat_edata.txt", sep = ""))

# make PCA plot
## find weightings of principal components
x = as.data.frame(combat.edata) #[qValuesComBat < .01,]
x = x[!is.na(rowSums(x)),]
pca = prcomp(t(x))
fac = factor(condition)
colours = brewer.pal(nlevels(fac), "Paired") 
#colours = hmcol
xyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$x), pch = 16, cex = 2, aspect = "fill", 
       col = colours, main = draw.key(key = list(rect = list(col = colours), text = list(levels(fac)), rep = FALSE)))#, 
   #    xlim = c(-20,20), ylim = c(-30, 15))

colours = colours[c(2,4,1,3)]


scatterplot3d(pca$x[,"PC1"], pca$x[,"PC2"], pca$x[,"PC3"], color = colours[as.numeric(as.factor(conditionAll[embryosSamePatient]))],
              xlab = "PC1", ylab = "PC2", zlab = "PC3", angle = 135)

colours = c("#00009b", "#009b00")
plot3d(pca$x[,"PC1"], pca$x[,"PC2"], pca$x[,"PC3"], col = colours[as.numeric(as.factor(conditionAll[embryosSamePatient]))], 
       size = 20, box = FALSE, xlab = "", ylab = "", zlab = "")
decorate3d(xlim = c(-30, 30), xlab = "", ylab = "", zlab = "", box = FALSE)

ggplot(as.data.frame(pca$x), aes(PC1,PC2, color=conditionAll[embryosSamePatient])) + geom_point(size=5) + geom_text(hjust=-.5, vjust=-.5,aes(label=as.character(1:22)[embryosSamePatient]))


# re-plot phylo
hmcol = c("#009b00", "#009b00", "#00009b", "#00009b", "#00009b", "#00009b", "#00009b")
hc = hclust(dist(t(combat.edata)))
# hc = hclust(dist(t(log2(cpm(yDiff)+1))))
#plot(as.phylo(hc), tip.color=hmcol[lane-1])
plot(as.phylo(hc), tip.color=hmcol[ceiling(decDist[embryosSamePatient])])



# calculate average logFC
avgGood = rowMeans(combat.edata[,conditionAll[embryosSamePatient] == "good"])
avgBad = rowMeans(combat.edata[,conditionAll[embryosSamePatient] == "bad"]) 
logFC = rowMeans(combat.edata[,conditionAll[embryosSamePatient] == "good"]) - 
                        rowMeans(combat.edata[,conditionAll[embryosSamePatient] == "bad"])
names(logFC) = row.names(combat.edata)
DEnames = names(qValuesComBat[qValuesComBat < .01])
combat.edata.withq = as.data.frame(cbind(combat.edata, qValuesComBat, avgGood, avgBad, logFC))
combat.edata.DE = combat.edata[DEnames,]

# plot logFC

exprMeans = as.data.frame((rowMeans(combat.edata[namesSymbol,conditionAll[embryosSamePatient] == "good"]) + 
                            rowMeans(combat.edata[namesSymbol,conditionAll[embryosSamePatient] == "bad"]))/2)

exprsToPlot = cbind(exprMeans, logFC[namesSymbol])
colorMat = rep(1, length(DEnames.qvals[namesSymbol]))
colorMat[DEnames.qvals[namesSymbol] < .01] = colorMat[DEnames.qvals[namesSymbol] < .01] + 1
colorsForPlot = c("#636363", "#3182bd")
plot(exprsToPlot, main="Smear Plot", xlab = "mean expr", ylab = "log fold change (viable vs nonviable) ", col = colorsForPlot[colorMat])
abline(h = 0, lwd = 2, col = "black", lty = "dashed")
abline(h = 1, lwd = 1, col = "black")
abline(h = -1, lwd = 1, col = "black")


# classify DE genes into categories
# For the reverse map:
x <- org.Hs.egSYMBOL2EG
# Get the entrez gene identifiers that are mapped to a gene symbol
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])


# DEnames.logical = logical vector signifying whether a given gene is DE or not. row.names contains ALL genes 
DEnames.qvals = qValuesComBat
DEnames.logical = as.integer(DEnames.qvals < .01)
names(DEnames.logical) = names(DEnames.qvals)
namesSymbol = names(DEnames.qvals)

# get entrezID for each gene symbol name
# remove gene symbol names for which no entrez ID is found
namesEntrez = unlist(sapply(namesSymbol, getEntrezID, mapping = xx))

DEnames.logical = DEnames.logical[names(namesEntrez)]
names(DEnames.logical) = names(namesEntrez)
DEnames.logical = DEnames.logical[!is.na(DEnames.logical)]
namesEntrez = namesEntrez[names(DEnames.logical)]
names(DEnames.logical) = namesEntrez

DEnames.logFC = logFC[names(namesEntrez)]
names(DEnames.logFC) = names(namesEntrez)
namesSymbol = names(namesEntrez)
DEnames.qvals = DEnames.qvals[names(namesEntrez)]


# package goseq
pwf = nullp(DEnames.logical, "hg18", "knownGene")
GO.wall = goseq(pwf, "hg18", "knownGene")
enriched.GO = GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method = "BH") < .9]
head(enriched.GO)
go = getgo(names(DEnames.logical), "hg18", "knownGene")


# get lists of upregulated and downregulated genes
namesSymbol.DE = names(qValuesComBat[qValuesComBat < .01])
logFC.DE = DEnames.logFC[namesSymbol.DE]
names(logFC.DE) = namesSymbol.DE
upGenesDE = names(logFC.DE[logFC.DE > 0])
downGenesDE = names(logFC.DE[logFC.DE < 0])
#write.table(as.character(upGenesDE), paste(baseDataDirectory, "/edgeR/ComBat_SVA/upGenes_qVals_pt01.txt", sep = ""), row.names=FALSE, quote=FALSE)
#write.table(as.character(downGenesDE), paste(baseDataDirectory, "/edgeR/ComBat_SVA/downGenes_qVals_pt01.txt", sep = ""), row.names=FALSE, quote=FALSE)
# write.table(names(DEnames.qvals[DEnames.qvals < .02]), paste(baseDataDirectory, "/edgeR/ComBat_SVA/allGenes_qvals_pt02.txt", sep = ""), 
#              row.names=FALSE, quote=FALSE, col.names = FALSE)

# upGenes cell cycle genes
c6 = c("GO:0051209")#, "GO:0006325")#, "GO:0031062")#, "GO:0007127", "GO:0007135")#, "GO:0072687", "GO:0000212", "GO:0044779")

# finds indices in list of ALL genes which contain desired GO terms
# so we take the intersection of genes containing desired GO terms, and genes that are differentially expressed
genes.ind = sapply(go, findGenesWithGoTerms, desiredGoList = c6)
genes = namesSymbol[genes.ind]#intersect(namesSymbol[genes.ind], namesSymbol.DE)
genes

write.table(genes, paste(baseDataDirectory, "/edgeR/Cytoscape/releaseOfCalciumIonGenes_Symbols.txt", sep = ""), row.names = FALSE, col.names = FALSE, quote = FALSE)

length(DEnames.qvals[genes][DEnames.qvals[genes] < .01])
length(DEnames.qvals[genes])
length(DEnames.qvals[genes][DEnames.qvals[genes] < .01])/length(DEnames.qvals[genes])
DEnames.qvals[genes][DEnames.qvals[genes] < .01]


# now take average and standard deviation logFC for all those genes
logFC = rowMeans(combat.edata[genes,conditionAll[embryosSamePatient] == "good"]) - 
  rowMeans(combat.edata[genes,conditionAll[embryosSamePatient] == "bad"])
logFC[DEnames.qvals[genes] < .05]


logFC.mean = mean(rowMeans(combat.edata[genes,conditionAll[embryosSamePatient] == "good"]) - 
                       rowMeans(combat.edata[genes,conditionAll[embryosSamePatient] == "bad"]))
logFC.std = sqrt(var(rowMeans(combat.edata[genes,conditionAll[embryosSamePatient] == "good"]) - 
                          rowMeans(combat.edata[genes,conditionAll[embryosSamePatient] == "bad"])))
plot(logFC)
abline(h = 0, lwd = 2, col = "blue", lty = "dashed")



# go through all gene ontology terms (1 through 99,999) and find genes associated with that term
# see if the DE genes among those genes are mostly up or down regulated 

percent.misregulated = rep(0,75500)
numGenesInCategory = rep(0,75500)

for (i in 1:32000){
  
  iterChar = as.character(i)
  addString = paste(as.character(rep(0,5-nchar(iterChar))), sep = "", collapse = "")
  currGOstring = paste("GO:00", addString, iterChar, sep="")
  
  genes.ind = sapply(go, findGenesWithGoTerms, desiredGoList = currGOstring)
  genes = namesSymbol[genes.ind]
  numGenesInCategory[i] = length(genes)
  names(percent.misregulated)[i] = currGOstring
  percent.misregulated[i] = length(DEnames.qvals[genes][DEnames.qvals[genes] < .01])/length(DEnames.qvals[genes])
  
#   # now take average and standard deviation logFC for all those genes
#   logFC = rowMeans(combat.edata[genes,conditionAll[embryosSamePatient] == "good"]) - 
#     rowMeans(combat.edata[genes,conditionAll[embryosSamePatient] == "bad"])

  if (i %% 1000 == 0){
    print(i)
    print(currGOstring)
    print(percent.misregulated[i])
  }

#   if (length(logFC[logFC > 0]) > (.7*length(logFC)) && length(logFC) > 5){
#     print(currGOstring)
#     print("up")
#     print(length(logFC))
#     
#   }
#   
#   if (length(logFC[logFC > 0]) < (.3*length(logFC)) && length(logFC) > 5){
#     print(currGOstring)
#     print("down")
#     print(length(logFC))
#     
#   }
  
}

# sort percent.misregulated
names(numGenesInCategory) = names(percent.misregulated)
percent.mis.sort = sort(percent.misregulated[numGenesInCategory > 50], decreasing = TRUE)
numGenes.sort = numGenesInCategory[names(percent.mis.sort)]

# scatterplot with 2 genes
allGeneNames = row.names(combat.edata)
embryoViability = conditionAll[embryosSamePatient]
namesSymbol[DEnames.qvals < .05 & namesSymbol < "MB" & namesSymbol > "MA"]
allGeneNames[allGeneNames < "J" & allGeneNames > "I"]
ggplot(as.data.frame(t(combat.edata)), aes(AURKA,CCNA2, color=embryoViability)) + 
  geom_point(size=5) + geom_text(hjust=-.5, vjust=-.5,aes(label=as.character(1:22)[embryosSamePatient])) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size=20), 
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))



# bar chart with names of all genes in a particular list / GO category
genesToPlot = c("IGF1R", "MAPK1", "MAP2K1", "PRKACA", "YWHAG", "CDK1","MAD2L1", "SKP2", "CCNA2", "CCNE1",
                "PPP2R5A", "ESPL1", "RAD21",  "EME1", "PTTG1", "BTRC", "PPP1CC", "SLK", "PIM1",
                "CDC25C", "FBXO5", "ANAPC4", "SMC2", "SMC3",   "BUB1")
genesToPlot = c("B4GALT1", "CALML4", "CCT2", "CCT5", "CCT6A", "CD9", "CDK1", "CRCP", "ITPR1",
                 "NOX5", "OOEP", "PARK7", "PLCB1", "PRKCA", "PRKCG", "PRKCI", "PLCZ1", "STX2", "WEE2", "ZP3", "ZP4")

goodMeans = rowMeans(combat.edata[genesToPlot,embryoViability == "good"])
badMeans = rowMeans(combat.edata[genesToPlot,embryoViability == "bad"])
goodStds = sqrt(rowVars(combat.edata[genesToPlot, embryoViability == "good"]))
badStds = sqrt(rowVars(combat.edata[genesToPlot, embryoViability == "bad"]))

edata.means = as.data.frame(cbind(goodMeans, badMeans))
edata.stds = as.data.frame(cbind(goodStds, badStds))

mp = barplot2(as.matrix(t(edata.means)), beside=TRUE, las = 2, col = c("#009b00", "#00009b"), ylim = c(0,15),
              plot.ci = TRUE, ci.l =t(edata.means) - t(edata.stds), ci.u = t(edata.means) + t(edata.stds))

legend(1,14, c("viable", "nonviable"), fill = c("#009b00", "#00009b"))
box()


########
# Filter genes by highest within-group mean/variance
conditionNum = c(0,2,0,0,0,2,2,2,2,0,1,1,1,1,1,1,3,3,3,3,3,1)
BCV_b1 = as.data.frame(sqrt(rowVars(cpm(yDiff)[,conditionNum == 0])) / rowMeans(cpm(yDiff)[,conditionNum == 0]))
BCV_g1 = as.data.frame(sqrt(rowVars(cpm(yDiff)[,conditionNum == 2])) / rowMeans(cpm(yDiff)[,conditionNum == 2]))
BCV_b2 = as.data.frame(sqrt(rowVars(cpm(yDiff)[,conditionNum == 1])) / rowMeans(cpm(yDiff)[,conditionNum == 1]))
BCV_g2 = as.data.frame(sqrt(rowVars(cpm(yDiff)[,conditionNum == 3])) / rowMeans(cpm(yDiff)[,conditionNum == 3]))

# average difference between good and bad divided by standard deviation of that batch
d1 = as.data.frame(abs(rowMeans(cpm(yDiff)[,conditionNum == 0]) - rowMeans(cpm(yDiff)[,conditionNum == 2])) / 
                     (rowMeans(cpm(yDiff)[,conditionNum == 0 | conditionNum == 2])))
d2 = as.data.frame(abs(rowMeans(cpm(yDiff)[,conditionNum == 1]) - rowMeans(cpm(yDiff)[,conditionNum == 3])) / 
                     (rowMeans(cpm(yDiff)[,conditionNum == 1 | conditionNum == 3])))

meanBCV = as.data.frame(rowMeans(cbind(BCV_b1, BCV_g1, BCV_b2, BCV_g2)))
x = as.data.frame(log2(cpm(yDiff)[order(meanBCV, decreasing = FALSE),]+1))
x = x[1:1000,]
pca = prcomp(t(x))
fac = factor(conditionDiff)
colours = brewer.pal(nlevels(fac), "Paired") 
ggplot(as.data.frame(pca$x), aes(PC1,PC2, color=conditionDiff)) + geom_point(size=5) + geom_text(hjust=-.5, vjust=-.5,aes(label=as.character(1:22)))




scatterplot3d(pca$x[,"PC1"], pca$x[,"PC2"], pca$x[,"PC3"], color = colours)


colorVec = 1:22
colorVec[conditionDiff == "bad1"] = colours[1]
colorVec[conditionDiff == "bad2"] = colours[2]
colorVec[conditionDiff == "good1"] = colours[3]
colorVec[conditionDiff == "good2"] = colours[4]
scatterplot3d(pca$x[,"PC1"], pca$x[,"PC2"], pca$x[,"PC3"], color = colorVec)




#######################################################
# Make Venn Diagram from Grondahl and Fragouli data
########################################################

genesG = c("ACVR1B", "ESPL1", "MAD2L1", "EME1", "ANAPC4", "HAS3", "UBQLN1", "UBE2D1", 
           "USP2", "USP42", "USP9X", "BUB1")

logFC.G1 = rowMeans(combat.edata[genesG,conditionAll[embryosSamePatient] == "bad"]) - 
                  rowMeans(combat.edata[genesG,conditionAll[embryosSamePatient] == "good"])
logFC.G2 = log2(c(1.5, 1.6, 1.6, 1/1.6, 1/1.6, 1.6, 1/1.6, 1/1.9, 2.9, 1.6, 1/1.8, 1/1.5))
logFC.G = as.data.frame(cbind(logFC.G1, logFC.G2))

mp = barplot2(as.matrix(t(logFC.G)), beside=TRUE, las = 2, col = c("#d7191c", "#2b83ba"))

legend(1,1.5, c("our data", "Grondahl 2010"), fill = c("#d7191c", "#2b83ba"))
box()



genesF = c("ATG7", "ASPM", "LNX2", "ELF1", "PCSK5", "RCCD1", "TMEM87A", "PDHB")

logFC.F1 = rowMeans(combat.edata[genesF,conditionAll[embryosSamePatient] == "bad"]) - 
  rowMeans(combat.edata[genesF,conditionAll[embryosSamePatient] == "good"])
logFC.F2 = log2(c(1/1.5, 1/1.5, 1/1.5, 1/1.5, 1.5, 1.5, 1.5, 1/1.5))
logFC.F = as.data.frame(cbind(logFC.F1, logFC.F2))

mp = barplot2(as.matrix(t(logFC.F)), beside=TRUE, las = 2, col = c("#d7191c", "#2b83ba"), ylim = c(-.7, 1))

legend(1,.9, c("our data", "Fragouli 2010"), fill = c("#d7191c", "#2b83ba"))
box()



