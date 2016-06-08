# Damon Polioudakis
# 2016-06-07
# Number of reads per cell
################################################################################

rm(list=ls())

library(xlsx)

# Load data and assign variables

# RNA STAR Stats Lane 1 
stStatsL1dF <- read.table("../metadata/RNAstar_Stats_SxaQSEQsXap089L2.txt"
                          , sep = "\t", header = TRUE)
# RNA STAR Stats Lane 2
stStatsL2dF <- read.table("../metadata/RNAstar_Stats_SxaQSEQsXbp060L2.txt"
                          , sep = "\t", header = TRUE)
################################################################################

# Combine total reads stats from Lane 1 and Lane 2
totReads <- (stStatsL1dF$Number.of.input.reads
             + stStatsL2dF$Number.of.input.reads)

mnReads <- mean(totReads)

# Plot as histogram
pdf("../analysis/graphs/Number_Reads_Per_Cell_Histogram.pdf")
hist(totReads, breaks = 20, col = "blue", xlab = "Number of reads per cell"
     , main = paste("Number of Reads Per Cell", "\nMean: ", round(mnReads,1), sep=""))
dev.off()