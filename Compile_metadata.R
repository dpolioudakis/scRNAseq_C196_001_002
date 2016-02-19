
# Damon Polioudakis
# 2016-02-18
# Add some RNA STAR statistics and Picard Statistics to sample metadata
# Gene Length and GC bias were calculated with Calc_Gene_Length_and_GC.R
################################################################################

rm(list=ls())

library(xlsx)

# Load data and assign variables

# Picard Sequencing Statistics
picStatsDF <- read.csv("../metadata/PicardToolsQC.csv")

# RNA STAR Stats for Run 1 
stStats1dF <- read.table("../metadata/RNAstar_Stats_SxaQSEQsXap089L2.txt"
                         , sep = "\t", header = TRUE)
# RNA STAR Stats for Run 2
stStats2dF <- read.table("../metadata/RNAstar_Stats_SxaQSEQsXbp060L2.txt"
                         , sep = "\t", header = TRUE)

# Metadata
metDatDF <- read.xlsx("../metadata/PercentageReadsMapping.xlsx", 1)
libBatch1DF <- read.xlsx("../metadata/LibraryBatch_Well_ID.xlsx", 1, header = FALSE)
libBatch2DF <- read.xlsx("../metadata/LibraryBatch_Well_ID.xlsx", 2, header = FALSE)

# GC content and average length for each sample
load("../analysis/tables/Avg_Gene_Length_and_GC.rda")
avgLengthDF <- data.frame(AvgGeneLength = avgLength)
avgGCdF <- data.frame(GCcontent = avgGC)
################################################################################

### Format data and add to metadata

##
# Remove statistics from Run 1 from metadata
metDatDF <- subset(metDatDF, select = -c(NA., uniquelymappedreads
                              , mappedtomultipleloci, NumReads))

##
# Calculate and add percent duplication data to metadata table

# Calculate percent duplication
pctDupDF <- data.frame(Pct_Duplicates =
                  picStatsDF$READ_PAIR_DUPLICATES / picStatsDF$PF_READS_ALIGNED
                  , row.names = picStatsDF$X)
# And add percent duplication data to metadata table
metDatDF <- merge(x = metDatDF, y = pctDupDF, by.x = "CellID", by.y = "row.names")

##
# Calculate and add total reads, percent uniquely mapped reads, percent mapped
# to multiple loci, and percent unmapped reads to metadata table

# Total Reads
totReads <- (stStats1dF$Number.of.input.reads
                + stStats2dF$Number.of.input.reads)
metDatDF$Total_Reads <- totReads

# Uniquely Mapped - Convert % to numeric and multiply by number of reads in lane
umReads <- (((as.numeric(sub("%", "", stStats1dF$Uniquely.mapped.reads..))/100) * stStats1dF$Number.of.input.reads)
         + (((as.numeric(sub("%", "", stStats2dF$Uniquely.mapped.reads..))/100) * stStats2dF$Number.of.input.reads)))
# Percent Uniquely Mapped
pctUMRds <- umReads / totReads
metDatDF$Pct_Uniquely_Mapped <- pctUMRds

# Total mapped to multiple loci
totMpMtLoci <- (stStats1dF$Number.of.reads.mapped.to.multiple.loci
              + stStats2dF$Number.of.reads.mapped.to.multiple.loci)
# Percent mapped to multiple loci
pctMpMtLoci <- totMpMtLoci / totReads
metDatDF$Pct_Mapped_Multiple_Loci <- pctMpMtLoci

# Percent unmapped reads
Calc_Unmapped <- function (data) {
  # Convert % to numeric and multiply by number of reads in lane
  df <- apply(data[ ,c("X..of.reads.unmapped..too.many.mismatches"
                                      , "X..of.reads.unmapped..too.short"
                                      , "X..of.reads.unmapped..other")]
              , 2, function(column) {as.numeric(sub("%", "", column)) / 100})
  # Multiple by number of reads in lane to convert to read number
  df <- apply(df, 2, function(column) {column * data$Number.of.input.reads})
  apply(df, 1, sum)
  }
# Lane 1
unMap1 <- Calc_Unmapped(stStats1dF)
# Lane 2
unMap2 <- Calc_Unmapped(stStats2dF)
# Percent unmapped reads
pctUnMap <- (unMap1 + unMap2) / totReads
metDatDF$Pct_Unmapped <- pctUnMap

##
# Add batch ID to metadata

batc2ID <- gsub("^ID:2-", "", as.vector(t(libBatch2DF)))
batc2ID <- gsub("^ID:1-", "1_", batc2ID)
batc2ID <- gsub("_[[:digit:]].*$", "", batc2ID)
batc1ID <- gsub("^ID:2-", "", as.vector(t(libBatch1DF)))
batc1ID <- gsub("^ID:1-", "1_", batc1ID)
batc1ID <- gsub("_[[:digit:]].*$", "", batc1ID)
batchIDdF <- rbind(data.frame(LibraryBatch = "1", CaptureWellID = batc1ID)
                   , data.frame(LibraryBatch = "2", CaptureWellID = batc2ID))
metDatDF <- merge(x = metDatDF, y = batchIDdF,
                  by.x = "CaptureWellID", by.y = "CaptureWellID")

##
# Add average gene length to metadata

row.names(avgLengthDF) <- gsub("_e.*", "", row.names(avgLengthDF))
avgLengthDF$CellID <- gsub("_e.*", "", row.names(avgLengthDF))
avgLengthDF <- avgLengthDF[2:(nrow(avgLengthDF)-1), ]
metDatDF <- merge(x = metDatDF, y = avgLengthDF, by.x = "CellID", by.y = "CellID")

##
# Add average GC content to metadata

row.names(avgGCdF) <- gsub("_e.*", "", row.names(avgGCdF))
avgGCdF$CellID <- gsub("_e.*", "", row.names(avgGCdF))
avgGCdF <- avgGCdF[2:(nrow(avgGCdF)-1), ]
metDatDF <- merge(x = metDatDF, y = avgGCdF, by.x = "CellID", by.y = "CellID")

##
# Write out to tab separated table
head(metDatDF)
write.table(metDatDF, file = paste0("../metadata/Compiled_Metadata_"
                                , format(Sys.time(), "%Y%m%d"), ".txt")
          , quote = FALSE, row.names = FALSE, sep = "\t")
