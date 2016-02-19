
# Damon Polioudakis
# 2016-02-16
# Plot % Unmapped Reads vs Mapped Reads and ERCC % vs Mapped Reads

# Inputs
  # Metadata
  # Sequencing Picard Statistics
# Outputs

################################################################################

rm(list=ls())
sessionInfo()

library(ggplot2)

# Picard Sequencing Statistics
picStatsDF <- read.csv("../metadata/PicardToolsQC.csv")

# Metadata
metDatDF <- read.table("../metadata/combined_DP_2016-01-28.txt", header = TRUE)

exDatDF <- read.csv("../data/htseq/merged/Exprs_HTSCexon.csv", row.names = 1)

colnames(picStatsDF)

exDatDF <- exDatDF[ ,order(colnames(exDatDF))]
colnames(exDatDF)

picStatsDF <- picStatsDF[order(picStatsDF$X), ]
picStatsDF$X

# pctUM <- picStatsDF$PF_READS_ALIGNED/picStatsDF$TOTAL_READS
nMapGenes <- apply(head(exDatDF, -97), 2, sum)
nMapERCC <- apply(head(tail(exDatDF, 97), -5), 2, sum)

pctERCC <- nMapERCC/picStatsDF$PF_READS_ALIGNED
round(pctERCC, 2)

plot(picStatsDF$PF_READS_ALIGNED/(10^6), pctERCC)
