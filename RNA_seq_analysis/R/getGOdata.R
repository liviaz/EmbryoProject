
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
go = getgo(names(DEnames.logical), "hg18", "knownGene")
