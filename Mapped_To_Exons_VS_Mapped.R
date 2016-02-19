# Damon Polioudakis
# 2016-02-10
# Plot bulk RNAseq of VZ and CP from Luis and Jason's ATAC versus pooled
# scRNAseq VZ and CP

# Inputs
#   HTseq counts for bulk RNAseq VZ and CP from Luis and Jason's ATAC
#   HTseq counts for scRNAseq C196-001_002

# Outputs

################################################################################

rm(list=ls())
sessionInfo()

library(ggplot2)

# Load data and assign variables

exDatDF <- read.csv("../data/htseq/merged/Exprs_HTSCexon.csv", row.names = 1)

# Picard Sequencing Statistics
picStatsDF <- read.csv("../metadata/PicardToolsQC.csv")

colnames(scExDatDF)
picStatsDF$X

# Remove ERCCs
exDatDF <- head(exDatDF, -97)


# Order
exDatDF <- exDatDF[ ,order(colnames(exDatDF))]
colnames(exDatDF)

picStatsDF <- picStatsDF[order(picStatsDF$X), ]
picStatsDF$X

nMapToExons <- apply(exDatDF, 2, sum)
cor(picStatsDF$PF_READS_ALIGNED, nMapToExons)
plot(picStatsDF$PF_READS_ALIGNED, nMapToExons)
