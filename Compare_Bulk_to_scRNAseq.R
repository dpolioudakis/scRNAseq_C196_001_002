
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

### Load data and assign variables

## Load data

buExDatDF <- read.csv("../data/htseq/bulk_VZ_CP_from_ATAC/Exprs_HTSCexon.csv"
                     , row.names = 1)

scExDatDF <- read.csv("../data/htseq/merged/Exprs_HTSCexon.csv", row.names = 1)

# Picard Sequencing Statistics - scRNAseq
picStatsScDF <- read.csv("../metadata/PicardToolsQC.csv")

# Picard Sequencing Statistics - scRNAseq
picStatsBuDF

# GC and Gene Length
load("../../source/ENSEMBLhg19_UnionAnno.rda")
lengthGCdF <- ENSEMBLhg19.70UnionAnno
rm(ENSEMBLhg19.70UnionAnno)


## Variables

graphCodeTitle <- "Compare_Bulk_to_scRNAseq.R"
outGraph <- "../analysis/graphs/Compare_Bulk_to_scRNAseq_"

################################################################################

### Functions


################################################################################

# Remove ERCCs
buExDatDF <- head(buExDatDF, -97)
scExDatDF <- head(scExDatDF, -97)
################################################################################

# Average counts for bulk for each transcript and counts of scRNAseq for each
# transript

# Pool scRNAseq
pScEx <- apply(scExDatDF, 1, sum)
head(pScEx)
tail(pScEx)

# Mean expression for bulk for each gene
mnBuEx <- apply(buExDatDF, 1, mean)

# Median for bulk for each gene
mdBuEx <- apply(buExDatDF, 1, median)
################################################################################

## Log2 (counts + 1) Pooled scRNAseq vs bulk

# Spearman correlation
sprCor <- round(cor(pScEx, mnBuEx, method = "spearman"), 2)

# Format for ggplot2
ggDF <- data.frame(Pooled = log(pScEx + 1, 2), Bulk = log(mnBuEx + 1, 2))

ggplot(ggDF, aes(x = Pooled, y = Bulk)) +
  geom_point(alpha = 0.5, shape = 1) +
  theme_bw(base_size = 18) +
  ylab("Bulk: log2(Mean Counts + 1)") +
  xlab("Pooled: log2(Mean Counts + 1)") +
  ggtitle(paste0(graphCodeTitle
          , "\nPooled scRNAseq vs Bulk RNAseq - Human Fetal Brain VZ and CP"
          , "\n Mean of counts across samples"
          , "\nSpearman correlation: ", sprCor))
ggsave(paste0(outGraph, "Bulk_vs_Pooled_log2p1.pdf"))
################################################################################

### Filter, read depth normalize, and pool in groups bulk and scRNAseq

## Filter for number of mapped reads > 1.5*10^6

# Identify Cells IDs with number of mapped reads > 1.5*10^6
pfCells <- picStatsScDF$X[picStatsScDF$PF_READS_ALIGNED > 1.5*10^6]
# Filter expression dataframe for cell IDs with number of mapped reads > 1.5*10^6
ftScExDatDF <- scExDatDF[ ,colnames(scExDatDF) %in% pfCells]


## Randomly split into 5 pooled groups

set.seed(11)
ncol(ftScExDatDF)
rNums <- sample(1:130, 130, replace = FALSE)
rndmGroups <- split(rNums, ceiling(seq_along(rNums) / (length(rNums) / 5)))
pScExDF <- data.frame(lapply(rndmGroups
                        , function (group) {apply(ftScExDatDF[ ,group], 1, sum)}))
head(pScExDF, 20)


## Read depth normalize

# Read depth normalize pooled scRNAseq
rDep <- (apply(pScExDF, 2, sum) / 10^6)
pScRdnExDatDF <- pScExDF / rDep

# Read depth normalize bulk
rDep <- (apply(buExDatDF, 2, sum) / 10^6)
buRdNExDatDF <- buExDatDF / rDep
################################################################################

### Plot read depth normalized

# Boxplot read depth normalized bulk log2 (counts)
boxplot(log(data.frame(buRdNExDatDF, pScRdnExDatDF) + 1, 2), range = 0)


## Log2 mean bulk vs pooled TPM

# Mean counts for bulk RNAseq
mnBuEx <- apply(buRdNExDatDF, 1, mean)

# Mean counts for pooled scRNAseq groups
mnPdScEx <- apply(pScRdnExDatDF, 1, mean)

# Spearman correlation
sprCor <- round(cor(mnPdScEx, mnBuEx, method = "spearman"), 2)

# Format for ggplot2
ggDF <- data.frame(Pooled = log(mnPdScEx + 1, 2), Bulk = log(mnBuEx + 1, 2))

ggplot(ggDF, aes(x = Pooled, y = Bulk)) +
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


## Log2 median bulk vs pooled TPM

# Median counts for bulk RNAseq
mdBuEx <- apply(buRdNExDatDF, 1, median)

# Median counts for pooled scRNAseq groups
mdPdScEx <- apply(pScRdnExDatDF, 1, median)

# Spearman correlation
sprCor <- round(cor(mdPdScEx, mdBuEx, method = "spearman"), 2)

# Format for ggplot2
ggDF <- data.frame(Pooled = log(mdPdScEx + 1, 2), Bulk = log(mdBuEx + 1, 2))

ggplot(ggDF, aes(x = Pooled, y = Bulk)) +
  geom_point(alpha = 0.5, shape = 1) +
  stat_smooth() +
  theme_bw(base_size = 18) +
  ylab("Bulk: log2(Median TPM + 1)") +
  xlab("Pooled: log2(Median TPM + 1)") +
  ggtitle(paste0(graphCodeTitle
                 , "\nFive Pools scRNAseq vs Bulk RNAseq - Human Fetal Brain VZ and CP"
                 , "\nMedian of TPM across samples"
                 , "\nSpearman correlation: ", sprCor))
ggsave(paste0(outGraph, "Bulk_vs_Pools_TPM_median_log2p1.pdf"))


## MA Plot

# Added +1 to all counts to prevent Inf values
ggDF <- data.frame(Avg = (0.5*log((mnPdScEx + 1) * (mnBuEx + 1), 2))
                   , Log2Ratio = log(((mnPdScEx + 1) / (mnBuEx + 1)), 2))
head(ggDF)
ggplot(ggDF, aes(x = Avg, y = Log2Ratio)) +
  geom_point(shape = 1, alpha = 0.5) +
  stat_smooth() +
  geom_vline(xintercept = 3.5) +
  geom_hline(yintercept = -0.5) +
  geom_hline(yintercept = 0.5) +
  theme_bw(base_size = 14) +
  xlab("0.5*log2((Mean Pools TPM + 1) * (Mean Bulk TPM + 1))") +
  ylab("log2((Mean Pools TPM + 1) / (Mean Bulk TPM + 1))") +
  ggtitle(paste0(graphCodeTitle
                 , "\nMA Plot: Five Pools scRNAseq vs Bulk RNAseq - Human Fetal Brain VZ and CP"
                 , "\nMean of TPM across samples"))
ggsave(paste0(outGraph, "Bulk_vs_Pools_TPM_MAPlot_TPM_mean_p1.pdf"))


## MA Plot: Mean Average > 3.5 (~10 TPM)

# Added +1 to all counts to prevent Inf values
ggDF <- data.frame(Avg = (0.5*log((mnPdScEx + 1) * (mnBuEx + 1), 2))
                   , Log2Ratio = log(((mnPdScEx + 1) / (mnBuEx + 1)), 2))
ggDF <- ggDF[ggDF$Avg > 3.5, ]

head(ggDF)
ggplot(ggDF, aes(x = Avg, y = Log2Ratio)) +
  geom_point(shape = 1, alpha = 0.5) +
  stat_smooth() +
  theme_bw(base_size = 14) +
  xlab("0.5*log2((Mean Pools TPM + 1) * (Mean Bulk TPM + 1))") +
  ylab("log2((Mean Pools TPM + 1) / (Mean Bulk TPM + 1))") +
  ggtitle(paste0(graphCodeTitle
                 , "\nMA Plot: Five Pools scRNAseq vs Bulk RNAseq - Human Fetal Brain VZ and CP"
                 , "\nMean of TPM across samples"
                 , "\nFiltered for Mean Average > 3.5 (~10 TPM)"))
ggsave(paste0(outGraph, "Bulk_vs_Pools_TPM_MAPlot_AG3.5_TPM_mean_p1.pdf"))

  # stat_smooth(method = lm)
# + coord_cartesian(xlim = c(0, 500000), ylim = c(0, 2500))



# Standard deviation for bulk RNAseq
sdBuEx <- apply(buRdNExDatDF, 1, sd)
################################################################################


###### Artifactual looking points at top of ScCov are most likely due to 4 of the
# pooled groups having 0 counts and 1 group having 1 count, not sure what the
# 1 normalizes to after TPM, but that math should carry through


rDep <- (apply(ftScExDatDF, 2, sum) / 10^6)
scRdnExDatDF <- ftScExDatDF / rDep



Coef_Of_Var <- function(x) {
  sd(x) / mean(x)
}

buCov <- apply(buRdNExDatDF, 1, function(tpm) Coef_Of_Var(tpm))
pdScCov <- apply(pScRdnExDatDF, 1, function(tpm) Coef_Of_Var(tpm))
scCov <- apply(scRdnExDatDF, 1, function(tpm) Coef_Of_Var(tpm))

# Mean counts for bulk RNAseq
mnBuEx <- apply(buRdNExDatDF, 1, mean)
lg2mnBuEx <- log(mnBuEx, 2)

# Mean counts for pooled scRNAseq groups
mnPdScEx <- apply(pScRdnExDatDF, 1, mean)
lg2mnPdScEx <- log(mnPdScEx, 2)

# Mean counts for scRNAseq
mnScEx <- apply(scRdnExDatDF, 1, mean)
lg2mnScEx <- log(mnScEx, 2)

plot(lg2mnBuEx, buCov)

sel <- ! is.na(pdScCov)
pdScCov <- pdScCov[sel]
lg2mnPdScEx <- lg2mnPdScEx[sel]
plot(lg2mnPdScEx, pdScCov)


plot(lg2mnScEx, scCov)


################################################################################

### Subset transcripts

## Use MA plot values to subset

# Mean counts for bulk RNAseq
mnBuEx <- apply(buRdNExDatDF, 1, mean)
lg2mnBuEx <- log(mnBuEx, 2)

# Mean counts for pooled scRNAseq groups
mnPdScEx <- apply(pScRdnExDatDF, 1, mean)
lg2mnPdScEx <- log(mnPdScEx, 2)

# Added +1 to all counts to prevent Inf values
maDF <- data.frame(Avg = (0.5*log((mnPdScEx + 1) * (mnBuEx + 1), 2))
                   , Log2Ratio = log(((mnPdScEx + 1) / (mnBuEx + 1)), 2))

# Low expressed bias towards bulk
# Transcripts > 0.5 M (log ratio) towards bulk and > 3.5 A (mean average)
# (expression > ~11 TPM)
biasLExBuDF <- maDF[maDF$Avg < 3.5 & maDF$Log2Ratio < -0.5, ]
dim(biasLExBuDF) # 4100
biasLExBuDF <- merge(biasLExBuDF, lengthGCdF, by.x = "row.names", by.y = "row.names")
sum(biasLExBuDF$UnionExonLength) / nrow(biasLExBuDF)
sum(biasLExBuDF$UnionGCcontent, na.rm = TRUE) / nrow(biasLExBuDF)

# Low expressed bias towards scRNAseq
# Transcripts > 0.5 M (log ratio) towards pooled scRNAseq and > 3.5 A (mean
# average) (expression > ~11 TPM)
biasLExScDF <- maDF[maDF$Avg < 3.5 & maDF$Log2Ratio > 0.5, ]
dim(biasLExScDF) # 2461
biasLExScDF <- merge(biasLExScDF, lengthGCdF, by.x = "row.names", by.y = "row.names")
sum(biasLExScDF$UnionExonLength) / nrow(biasLExScDF)
sum(biasLExScDF$UnionGCcontent) / nrow(biasLExScDF)

# High expressed bias towards bulk
# Transcripts > 0.5 M (log ratio) towards bulk and > 3.5 A (mean average)
# (expression > ~11 TPM)
biasHExBuDF <- maDF[maDF$Avg > 3.5 & maDF$Log2Ratio < -0.5, ]
dim(biasHExBuDF) # 4100
biasHExBuDF <- merge(biasHExBuDF, lengthGCdF, by.x = "row.names", by.y = "row.names")
sum(biasHExBuDF$UnionExonLength) / nrow(biasHExBuDF)
sum(biasHExBuDF$UnionGCcontent) / nrow(biasHExBuDF)

# How expressed bias towards scRNAseq
# Transcripts > 0.5 M (log ratio) towards pooled scRNAseq and > 3.5 A (mean
# average) (expression > ~11 TPM)
biasHExScDF <- maDF[maDF$Avg > 3.5 & maDF$Log2Ratio > 0.5, ]
dim(biasHExScDF) # 2461
biasHExScDF <- merge(biasHExScDF, lengthGCdF, by.x = "row.names", by.y = "row.names")
sum(biasHExScDF$UnionExonLength) / nrow(biasHExScDF)
sum(biasHExScDF$UnionGCcontent) / nrow(biasHExScDF)

boxplot(biasLExBuDF$UnionExonLength, biasLExScDF$UnionExonLength
        , biasHExBuDF$UnionExonLength, biasHExScDF$UnionExonLength
        , col = c(2:5), outline = FALSE)

boxplot(biasLExBuDF$UnionGCcontent, biasLExScDF$UnionGCcontent
        , biasHExBuDF$UnionGCcontent, biasHExScDF$UnionGCcontent
        , col = c(2:5), outline = FALSE)

ggDF <- rbind(data.frame(Length = biasLExBuDF$UnionExonLength
                         , Gene_Set = "Low Expression, Higher Bulk")
             , data.frame(Length = biasLExScDF$UnionExonLength
                          , Gene_Set = "Low Expression, Lower Bulk")
             , data.frame(Length = biasHExBuDF$UnionExonLength
                          , Gene_Set = "Higher Expression, Higher Bulk")
             , data.frame(Length = biasHExScDF$UnionExonLength
                          , Gene_Set = "Higher Expression, Lower Bulk"))
ggplot(ggDF, aes(y = Length, x = Gene_Set, col = Gene_Set)) +
  geom_boxplot(outlier.shape = NA) +
  # Adjust limits after outlier removal
  coord_cartesian(ylim = range(boxplot(ggDF$Length, plot = FALSE)$stats) * c(.9, 1.1)) +
  ggtitle(paste0("Compare_Bulk_to_scRNAseq.R"
                 , "\nGene Length for Subsets of Genes with Biased Expression"
                 , "\nP-value Low Expression, Higher Bulk vs Higher Pooled scRNA-seq: ",
                 signif(t.test(x = biasLExBuDF$UnionExonLength, y = biasLExScDF$UnionExonLength)$p.value, 3)
                 ,"\nP-value High Expression, Higher Bulk vs Higher Pooled scRNA-seq: ",
                 signif(t.test(x = biasHExBuDF$UnionExonLength, y = biasHExScDF$UnionExonLength)$p.value, 3)
))
###### Add number of genes to graph (n)
#### Repeat for GC