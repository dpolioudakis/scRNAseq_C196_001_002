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

# Metadata - has total reads
metDatDF <- read.table("../metadata/Compiled_Metadata_20160218.txt"
                       , header = TRUE, row.names = 1, fill = TRUE, sep = "\t")

## Variables

graphCodeTitle <- "Mapped_To_Exons_VS_Mapped.R"
outGraph <- "../analysis/graphs/Mapped_To_Exons_VS_Mapped_"
################################################################################

# Remove ERCCs
exDatDF <- head(exDatDF, -97)

# Order
exDatDF <- exDatDF[ ,order(colnames(exDatDF))]
colnames(exDatDF)

picStatsDF <- picStatsDF[order(picStatsDF$X), ]

# Calc number mapped to exons
nMapToExons <- apply(exDatDF, 2, sum)
################################################################################

pdf(paste0(outGraph, "Aligned_Vs_Exons.pdf"))
plot(picStatsDF$PF_READS_ALIGNED, nMapToExons
     , xlab = "Aligned Reads"
     , ylab = "Mapped to Exons"
     , main = paste0(graphCodeTitle
                     , "\nComparing Number of Reads Aligned to Number Mapped to Exons"
                     , "\nPearson Correlation: "
                     , signif(cor(picStatsDF$PF_READS_ALIGNED, nMapToExons), 2))
)
dev.off()

pdf(paste0(outGraph, "Total_Vs_Exons.pdf"))
plot(metDatDF$Total_Reads, nMapToExons
     , xlab = "Total Reads"
     , ylab = "Mapped to Exons"
     , main = paste0(graphCodeTitle
                     , "\nComparing Total Number of Reads to Number Mapped To Exons"
                     , "\nPearson Correlation: "
                     , signif(cor(metDatDF$Total_Reads, nMapToExons), 2))
)
dev.off()

pdf(paste0(outGraph, "Aligned_Vs_Total.pdf"))
plot(picStatsDF$PF_READS_ALIGNED, metDatDF$Total_Reads
     , xlab = "Aligned Reads"
     , ylab = "Total Reads"
     , main = paste0(graphCodeTitle
                     , "\nComparing Number of Reads Aligned to Total Number of Reads"
                     , "\nPearson Correlation: "
                     , signif(cor(picStatsDF$PF_READS_ALIGNED, metDatDF$Total_Reads), 2))
)
dev.off()