# Damon Polioudakis
# 2016-08-01
# Plot Y chromosome and Xist expression for Fluidigm LT C196
################################################################################

rm(list=ls())
sessionInfo()

library(ggplot2)
library(reshape2)
require(biomaRt)

### Load data and assign variables

## Load data

scExDatDF <- read.csv("../data/htseq/merged/Exprs_HTSCexon.csv", row.names = 1)

# Picard Sequencing Statistics - scRNAseq
picStatsScDF <- read.csv("../metadata/PicardToolsQC.csv")


## Variables

outGraph <- "../analysis/graphs/Ychr_Xist_Expression_"
graphTitle <- "Ychr_Xist_Expression.R"
################################################################################

# Remove ERCCs
scExDatDF <- head(scExDatDF, -97)

# Identify Cells IDs with number of mapped reads > 1.5*10^6
pfCells <- picStatsScDF$X[picStatsScDF$PF_READS_ALIGNED > 1.5*10^6]
# Filter expression dataframe for cell IDs with number of mapped reads > 1.5*10^6
ftScExDatDF <- scExDatDF[ ,colnames(scExDatDF) %in% pfCells]

## Chomp . and everything after from ensembl IDs
row.names(ftScExDatDF) <- gsub("\\..*", "", row.names(ftScExDatDF))

## Convert Ensembl IDs to Gene Symbols

AddChromosomeHuman <- function (ensemblList) {
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  moduleGenes <- data.frame(ensemblList)
  # bioMart manual:
  #http://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf
  # Attribute: values to retrieve
  # Filters: input query type
  # Values: input query
  #ensembl <- useMart("ensembl")
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
  ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
  # Data frame of module Ensembl IDs and gene symbols
  ensemblGeneSymDF <- getBM(  attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name")
                              , filters = "ensembl_gene_id"
                              , values = moduleGenes
                              , mart = ensembl
  )
  ensemblGeneSymDF
}

# Dataframe of ensembl IDs and chromosome
hsEnsemblChrDF <- AddChromosomeHuman(row.names(ftScExDatDF))

# Add chromosome to expression dataframe
hsFtExDF <- merge(hsEnsemblChrDF, ftScExDatDF
                  , by.x = "ensembl_gene_id", by.y = "row.names")

## Y expression
pdf(paste0(outGraph, "Mean_Ychr_Expression.pdf"))
yHsFtExDF <- hsFtExDF[hsFtExDF$chromosome_name == "Y", ]
barplot(sort(apply(yHsFtExDF[ ,-c(1:3)], 2, mean))
        , xaxt = "n", xlab = "Cells"
        , ylab = "Expression (counts)"
        , main = paste0(graphTitle
                        , "\n", "Fluidigm LT mean Y chromosome expression"
                        , "\n", "> 1.5e6 human reads aligned"))
dev.off()

## Xist expression
pdf(paste0(outGraph, "Xist_Expression.pdf"))
# Human
xistHsM <- as.matrix(hsFtExDF[hsFtExDF$hgnc_symbol == "XIST", -c(1:3)])
barplot(xistHsM
        , xaxt = "n", xlab = "Cells"
        , ylab = "Expression (counts)"
        , main = paste0(graphTitle
                        , "\n", "Fluidigm LT Xist expression"
                        , "\n", "> 1.5e6 human reads aligned"))
dev.off()

## Log2 (expression + 1) by chromosome
ggDF <- melt(hsFtExDF)
# Remove biomart extra chromosomes
rmvChr <- c(grep("CHR*", unique(ggDF$chromosome_name), value = TRUE)
            , grep("KI27*", unique(ggDF$chromosome_name), value = TRUE))
ggDF <- ggDF[! ggDF$chromosome_name %in% rmvChr, ]
ggDF$value <- log(ggDF$value + 1, 2)
aggregate(ggDF$value, list(chr = ggDF$chromosome_name), mean)
ggplot(ggDF, aes(y = value, x = chromosome_name)) +
  geom_jitter(alpha = 0.25) +
  ylab("Log2(counts + 1)") +
  xlab("Chromosome") + 
  ggtitle(paste0(graphTitle
                 , "\n", "Fluidigm LT expression by chromosome"
                 , "\n", "> 1.5e6 human reads aligned"
                 , "\n"))
ggsave(paste0(outGraph, "Expression_By_Chromosome.png"), width = 9, height = 6)  