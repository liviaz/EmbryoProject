
## Starts off with DESeq, variance stabilization, etc

rm(list = ls())


if (Sys.info()["sysname"] == "Linux") {
  source("/host/Users/Livia/Desktop/IVF/Code/EmbryoProject/RNA_seq_analysis/R/LoadHTSeqCountData_Human_HiSeq.R")
} else {
  source("C:/Users/Livia/Desktop/IVF/Code/EmbryoProject/RNA_seq_analysis/R/LoadHTSeqCountData_Human_HiSeq.R")
}

# assemble count data set
conditionAll = c("bad", "good", "bad", "bad", "bad", "good", "good", "good", "good", "bad", "bad", "bad", "bad", "bad", 
                 "good", "bad", "good", "good", "good", "good", "good", "bad")

conditionDiff = c("bad1", "good1", "bad1", "bad1", "bad1", "good1", "good1", "good1", "good1", "bad1", "bad2", "bad2", "bad2", "bad2", 
                  "good2", "bad2", "good2", "good2", "good2", "good2", "good2", "bad2")
decDist = c(4.36, 0.7, 2.09, 3.19, 3.38, 1.08, 1.63, 0.50, 1.12, 3.90, 6.31, 4.04, 4.63, 2.87, 1.72, 2.55, 1.07, 0.94, 1.08, 1.46, 1.82, 2.57)

# 1. run edgeR
yAll = DGEList(counts = allCounts, group = conditionAll)
yAll = yAll[rowSums(allCounts > 20) > 4,]
yAll = calcNormFactors(yAll)
yAll = estimateDisp(y = yAll, design = model.matrix(~0 + factor(conditionAll)))


cdsAll = newCountDataSet(allCounts[rowSums(allCounts > 20) > 4,], conditionAll)
# 2. normalize
cdsAll = estimateSizeFactors(cdsAll)
# 3. estimate dispersions of each gene across all samples
cdsAll = estimateDispersions(cdsAll)
vsdAll = varianceStabilizingTransformation(cdsAll)

# add in decision boundary distances measured using Matlab
phenoData(vsdAll)$decDist = c(4.36, 0.7, 2.09, 3.19, 3.38, 1.08, 1.63, 0.50, 1.12, 3.90, 6.31, 4.04, 4.63, 2.87, 1.72, 2.55, 1.07, 0.94, 1.08, 1.46, 1.82, 2.57)
phenoData(vsdAll)$exptNum = c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2)
phenoData(vsdAll)$condition = conditionAll

# note which embryos had siblings
embryosSamePatient = c(3,4,6,7,8,10,11,12,13,14,15,16,17,18,20,21,22)

# create model matrix and null matrix
batch = exptNum - 1
patient = c(1,2,3,3,4,5,5,6,7,6,8,9,8,8,9,8,10,10,11,12,10,12) - 1
#mod = model.matrix(~as.factor(condition), data = pData(vsdAll))
mod = model.matrix(~as.factor(conditionAll)[embryosSamePatient], data = pData(vsdAll)[embryosSamePatient,])
mod0 = model.matrix(~1, data = pData(vsdAll)[embryosSamePatient,])

mod1 = model.matrix(~as.factor(patient[embryosSamePatient]))


# use ComBat
combat.edata = as.data.frame(ComBat(dat = log2(cpm(yAll)[,embryosSamePatient]+1), batch = patient[embryosSamePatient], mod = mod, 
                      par.prior = TRUE, prior.plots = TRUE))
# filter out genes with < 1 cpm average expression
combat.edata = combat.edata[rowMeans(combat.edata) > 0,]
# calculate pvalues and qvalues
pValuesComBat = f.pvalue(as.matrix(combat.edata), mod, mod0)
qValuesComBat = p.adjust(pValuesComBat, method = "BH")
write.table(as.character(names(qValuesComBat[qValuesComBat < .01])), paste(baseDataDirectory, "/edgeR/ComBat_SVA/qVals_lt_pt01_ComBat.txt", sep = ""), row.names=FALSE, quote=FALSE)
# write.table(as.character(names(qValuesSv[qValuesSv < .1])), paste(baseDataDirectory, "/edgeR/ComBat_SVA/qVals_lt_pt05_SVA.txt", sep = ""), row.names=FALSE, quote=FALSE)
# write.table(combat.edata, paste(baseDataDirectory, "/edgeR/ComBat_SVA/combat_edata.txt", sep = ""))
# write.table(DEnames.qvals, paste(baseDataDirectory, "/edgeR/ComBat_SVA/DEnames_qvals.txt", sep = ""))

# make PCA plot
## find weightings of principal components
x = as.data.frame(combat.edata) #[qValuesComBat < .01,]
x = x[!is.na(rowSums(x)),]
pca = prcomp(t(x))
fac = factor(conditionDiff[embryosSamePatient])
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



# calculate average logFC (fold-change)
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
c6 = c("GO:0030163", "GO:0042176", "GO:0070271",  "GO:0030163", "GO:0008104")

# finds indices in list of ALL genes which contain desired GO terms
# so we take the intersection of genes containing desired GO terms, and genes that are differentially expressed
genes.ind = sapply(go, findGenesWithGoTerms, desiredGoList = c6)
genes = namesSymbol[genes.ind]#intersect(namesSymbol[genes.ind], namesSymbol.DE)
genes


length(DEnames.qvals[genes][DEnames.qvals[genes] < .05])
length(DEnames.qvals[genes])
length(DEnames.qvals[genes][DEnames.qvals[genes] < .05])/length(DEnames.qvals[genes])
DEnames.qvals[genes][DEnames.qvals[genes] < .01]

# write.table(DEnames.qvals[genes][DEnames.qvals[genes] < .01], paste(baseDataDirectory, "/edgeR/Cytoscape/proteinFunctions_DE_q_lt01.txt", sep = ""), row.names = TRUE, col.names = FALSE, quote = FALSE)

# now take average and standard deviation logFC for all those genes
logFC = rowMeans(combat.edata[genes,conditionAll[embryosSamePatient] == "good"]) - 
  rowMeans(combat.edata[genes,conditionAll[embryosSamePatient] == "bad"])
logFC[DEnames.qvals[genes] < .01]


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


