# Damon Polioudakis
# 2016-01-03
# Graph of number of genes expressed in 2015-12 scRNAseq run SxaQSEQsXbp060L2
rm(list=ls())
sessionInfo()

# Input gene expression table (counts from HTseq)
exDat1DF <- read.csv("../data/Exprs_HTSCexon_Run1.csv", header = TRUE, row.names = 1)
exDat2DF <- read.csv("../data/Exprs_HTSCexon_Run2.csv", header = TRUE, row.names = 1)
fpm1DF <- read.csv("../data/FPM.csv", header = TRUE)
dim(exDat1DF)
dim(exDat2DF)


exDatDF <- Reduce('+', list(exDat2DF, exDat1DF))

fpmCells <- gsub("X", "Cell", colnames(fpm1DF))
exDatDF <- exDatDF[ , colnames(exDatDF) %in% fpmCells]
colnames() <- gsub("_.*", "", colnames(exDat1DF))
ex1fTdF <- exDat1DF[ , colnames(exDat1DF) %in% fpmCells]
ex2fTdF <- exDat2DF[ , colnames(exDat2DF) %in% fpmCells]


# Number of genes expressed
gExprd <- apply(exDatDF, 2, function(genes) length(subset(genes, genes > 5)))
gExprd <- sort(gExprd)
mgExprd <- mean(gExprd)

pdf("../analysis/graphs/Alignment_Stats_Graph_Genes_Expressed.pdf")
barplot(gExprd, xlab = "Cells", ylab = "Number of genes TPM > 0", xaxt = 'n'
        , main = paste("Kriegstein 2015: Number of Genes Expressed"
                       ,"\nMean: ", mgExprd, sep=""))
Axis(side = 1, labels = FALSE)
dev.off()


gExprd <- apply(exDatDF, 2, function(genes) length(subset(genes, genes >= 1)))
mgExprd <- mean(gExprd)
mean(gExprd)
median(gExprd)
pdf("../analysis/graphs/Histogram_Genes_Expressed.pdf")
png("../analysis/graphs/Histogram_Genes_Expressed.png")
hist(gExprd, breaks = 20, xlab = "Genes Detected with >=1 reads"
	, main = paste("Number of Genes Expressed", "\nMean: ", mgExprd, sep=""))
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




