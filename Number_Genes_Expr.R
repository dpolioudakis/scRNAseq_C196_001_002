# Damon Polioudakis
# 2016-01-03
# Graph of number of genes expressed in 2015-12 scRNAseq run C196-001
# SxaQSEQsXbp060L2 and 2015-12 scRNAseq additional sequencing C196-002
rm(list=ls())
sessionInfo()

# Input gene expression table (counts from HTseq)
exDat1DF <- read.csv("../data/Exprs_HTSCexon_Run1.csv", header = TRUE, row.names = 1)
exDat2DF <- read.csv("../data/Exprs_HTSCexon_Run2.csv", header = TRUE, row.names = 1)
fpm1DF <- read.csv("../data/FPM.csv", header = TRUE)
dim(exDat1DF)
dim(exDat2DF)

# Picard Sequencing Statistics - scRNAseq
picStatsDF <- read.csv("../metadata/PicardToolsQC.csv")

# Combine counts from scRNAseq Run1 and Run2
exDatDF <- Reduce('+', list(exDat2DF, exDat1DF))

# Calculate CPM
exDatDF <- exDatDF[ ,order(names(exDatDF))]
cpmDF <- exDatDF / (picStatsDF$PF_READS_ALIGNED / 10^6)

# Subset to cells remaining after Jason's ERCC filters
fpmCells <- gsub("X", "Cell", colnames(fpm1DF))
cpmDF <- cpmDF[ , colnames(cpmDF) %in% fpmCells]
ex1fTdF <- exDat1DF[ , colnames(exDat1DF) %in% fpmCells]
ex2fTdF <- exDat2DF[ , colnames(exDat2DF) %in% fpmCells]

# Number of genes expressed
gExprd <- apply(cpmDF, 2, function(genes) length(subset(genes, genes >= 1)))
mgExprd <- mean(gExprd)
mean(gExprd)
median(gExprd)
pdf("../analysis/graphs/Number_Genes_Detected_Histogram.pdf")
hist(gExprd, breaks = 20, col = "blue", xlab = "Number of genes with >=1 CPM"
	, main = paste("Number of Genes Detected Per Cell", "\nMean: ", mgExprd, sep=""))
dev.off()






gExprd <- apply(exDat2DF, 2, function(genes) length(subset(genes, genes > 1)))
mean(gExprd)
median(gExprd)
hist(gExprd, breaks = 20)

gExprd <- apply(exDat1DF, 2, function(genes) length(subset(genes, genes >= 1)))
mean(gExprd)
median(gExprd)
hist(gExprd, breaks = 20)




corR1R2 <- cor(exDat1DF, exDat2DF, method = "spearman")
corR1R2 <- as.matrix(corR1R2)

pdf("../analysis/graphs/Histogram_Genes_Expressed.pdf")
png("../analysis/graphs/Histogram_Genes_Expressed.png")
heatmap(as.matrix(corR1R2), Rowv = NA, Colv = NA)
dev.off()


corR1R2 <- cor(ex1fTdF, ex2fTdF, method = "spearman")
corR1R2 <- as.matrix(corR1R2)

pdf("../analysis/graphs/Histogram_Genes_Expressed_Ft.pdf")
png("../analysis/graphs/Histogram_Genes_Expressed_Ft.png")
heatmap(as.matrix(corR1R2), Rowv = NA, Colv = NA)
dev.off()




