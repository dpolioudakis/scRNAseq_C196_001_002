# Damon Polioudakis
# 2016-03-17
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

graphCodeTitle <- "Compare_Variance_Between_Pooled_Groups.R"
outGraph <- "../analysis/graphs/Compare_Variance_Between_Pooled_Groups_"

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
pScIDs <- lapply(rndmGroups
                 , function (group) colnames(ftScExDatDF)[group])
dim(pScExDF)


## Read depth normalize

# Read depth normalize pooled scRNAseq
sel <- lapply(pScIDs, function(group) picStatsScDF$X %in% group)
splitStats <- lapply(sel, function(group) picStatsScDF[group, ])
mappedReads <- sapply(splitStats, function(group) sum(as.numeric(group$PF_READS_ALIGNED)))
pScTpmExDatDF <- pScExDF / (mappedReads / 10^6)
################################################################################

## Plot Pool 1 vs Pool 2 TPM

# Pearson correlation
peaCor <- round(cor(pScTpmExDatDF$X1, pScTpmExDatDF$X2, method = "pearson"), 2)

# Format for ggplot2
ggDF <- data.frame(Pool1 = log(pScTpmExDatDF$X1 + 1, 2)
                   , Pool2 = log(pScTpmExDatDF$X2 + 1, 2))

ggplot(ggDF, aes(x = Pool1, y = Pool2)) +
  geom_point(alpha = 0.5, shape = 1) +
  stat_smooth() +
  theme_bw(base_size = 18) +
  ylab("Pool 2: log2(Mean TPM + 1)") +
  xlab("Pool 1: log2(Mean TPM + 1)") +
  ggtitle(paste0(graphCodeTitle
                 , "\nHuman Fetal Brain VZ and CP"
                 , "\nscRNAseq Split Into Two Groups and Pooled"
                 , "\nPearson correlation: ", peaCor))
ggsave(paste0(outGraph, "Pool1_Vs_2_TPM_log2p1.pdf"))


## Split Pool 1 TPM into 10 even ranges from max to min TPM and calculate
## correlation between Pool 1 and Pool 2 for that range of TPM values

df <- data.frame(Lg2Exp = log(pScTpmExDatDF + 1, 2))

splits <- seq(0, range(df$Lg2Exp.X1)[2], length.out = 11)

df$split <- with(df, cut(Lg2Exp.X1, 
                            breaks = splits, 
                            include.lowest = TRUE))
splitDF <- split(df, df$split)

# Pearson correlation for each range
pcorDF <- data.frame(sapply(splitDF, function(subgroup) cor(subgroup$Lg2Exp.X1
                                                        , subgroup$Lg2Exp.X2)))
colnames(pcorDF) <- "Pearson_Correlation"

# Number of genes with TPM for Pool 1 that falls in range
nGenes <- data.frame(Number_of_Genes = sapply(splitDF, nrow))

# Plot pearson correlation for each range
# Set factor to set order on x-axis
ggplot(pcorDF, aes(x = factor(row.names(pcorDF), levels = row.names(pcorDF))
                   , y =  Pearson_Correlation)) +
  geom_bar(stat = "identity") +
  geom_text(data = nGenes, aes(label = paste("n =", Number_of_Genes)
                                 , y = pcorDF$Pearson_Correlation + 0.05)) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("TPM ranges") +
  ylab("Pearson Correlation Between Pool 1 and Pool 2") +
  ggtitle(paste0(graphCodeTitle
                 , "\nHuman Fetal Brain VZ and CP"
                 , "\nscRNAseq Split Into Two Groups and Pooled"
                 , "\nPearson Correlation Between Pool 1 and Pool 2"))
ggsave(paste0(outGraph, "Ranges_PearCor_TPM_log2p1.pdf"))
