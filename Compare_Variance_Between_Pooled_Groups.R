# Damon Polioudakis
# 2016-02-10
# Compare variance between groups of pooled scRNAseq

################################################################################

rm(list=ls())
sessionInfo()

library(ggplot2)

### Load data and assign variables

## Load data

scExDatDF <- read.csv("../data/htseq/merged/Exprs_HTSCexon.csv", row.names = 1)
buExDatDF <- read.csv("../data/htseq/bulk_VZ_CP_from_ATAC/Exprs_HTSCexon.csv"
                      , row.names = 1)

# Picard Sequencing Statistics - scRNAseq
picStatsScDF <- read.csv("../metadata/PicardToolsQC.csv")

## Variables

graphCodeTitle <- "Compare_Bulk_to_scRNAseq.R"
outGraph <- "../analysis/graphs/Compare_Bulk_to_scRNAseq_"

################################################################################

## Format data

# Remove ERCCs
buExDatDF <- head(buExDatDF, -97)
scExDatDF <- head(scExDatDF, -97)
################################################################################

### Filter, read depth normalize, and pool in groups bulk and scRNAseq

## Filter for number of mapped reads > 1.5*10^6

# Identify Cells IDs with number of mapped reads > 1.5*10^6
pfCells <- picStatsScDF$X[picStatsScDF$PF_READS_ALIGNED > 1.5*10^6]
# Filter expression dataframe for cell IDs with number of mapped reads > 1.5*10^6
ftScExDatDF <- scExDatDF[ ,colnames(scExDatDF) %in% pfCells]


## Randomly split into 2 pooled groups

set.seed(11)
ncol(ftScExDatDF)
rNums <- sample(1:130, 130, replace = FALSE)
rndmGroups <- split(rNums, ceiling(seq_along(rNums) / (length(rNums) / 2)))
pScExDF <- data.frame(lapply(rndmGroups
                             , function (group) {apply(ftScExDatDF[ ,group], 1, sum)}))
head(pScExDF, 20)
dim(pScExDF)


## Read depth normalize

# Read depth normalize pooled scRNAseq
rDep <- (apply(pScExDF, 2, sum) / 10^6)
pScRdnExDatDF <- pScExDF / rDep
head(pScRdnExDatDF)
################################################################################


# Spearman correlation
sprCor <- round(cor(pScRdnExDatDF$X1, pScRdnExDatDF$X2, method = "spearman"), 2)

# Format for ggplot2
ggDF <- data.frame(Pool1 = log(pScRdnExDatDF$X1 + 1, 2)
                   , Pool2 = log(pScRdnExDatDF$X2 + 1, 2))

ggplot(ggDF, aes(x = Pool1, y = Pool2)) +
  geom_point(alpha = 0.5, shape = 1) +
  stat_smooth() +
  theme_bw(base_size = 18) +
  ylab("Bulk: log2(Mean TPM + 1)") +
  xlab("Pooled: log2(Mean TPM + 1)") +
  ggtitle(paste0(graphCodeTitle
                 , "\nFive Pools scRNAseq vs Bulk RNAseq - Human Fetal Brain VZ and CP"
                 , "\nMean of TPM across samples"
                 , "\nSpearman correlation: ", sprCor))
ggsave(paste0(outGraph, "Bulk_vs_Pools_TPM_log2p1.pdf"))
