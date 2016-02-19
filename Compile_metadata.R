

rm(list=ls())




# RNA STAR Stats for Run 1 
stStats1dF <- read.table("../metadata/RNAstar_Stats_SxaQSEQsXap089L2.txt"
                         , sep = "\t", header = TRUE)
# RNA STAR Stats for Run 2
stStats2dF <- read.table("../metadata/RNAstar_Stats_SxaQSEQsXbp060L2.txt"
                         , sep = "\t", header = TRUE)

# Metadata
metDatDF <- read.xlsx("../metadata/PercentageReadsMapping.xlsx", 1)
pctDupDF <- read.xlsx("../metadata/PercentDuplicates.xlsx", 1)
libBatch1DF <- read.xlsx("../metadata/LibraryBatch_Well_ID.xlsx", 1, header = FALSE)
libBatch2DF <- read.xlsx("../metadata/LibraryBatch_Well_ID.xlsx", 2, header = FALSE)

# GC content and average length for each sample
load("../analysis/tables/Gene_Length_GC_Bias.rda")
lenBiasDF <- data.frame(GeneLengthScore = lenBias)
gcBiasDF <- data.frame(GCcontent = gcBias)
################################################################################

### Format data and add to metadata

##
# And add percent duplication data to metadata table

# Edit CellIDs to be same format
pctDupDF <- data.frame(PERCENT_DUPLICATION = t(pctDupDF)[-1, 1])
row.names(pctDupDF) <- gsub("\\..*\\.txt", "", row.names(pctDupDF))
# And add percent duplication data to metadata table
metDatDF <- merge(x = metDatDF, y = pctDupDF, by.x = "CellID", by.y = "row.names")

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
# Add length bias to metadata

row.names(lenBiasDF) <- gsub("_e.*", "", row.names(lenBiasDF))
lenBiasDF$CellID <- gsub("_e.*", "", row.names(lenBiasDF))
lenBiasDF <- lenBiasDF[2:(nrow(lenBiasDF)-1), ]
metDatDF <- merge(x = metDatDF, y = lenBiasDF, by.x = "CellID", by.y = "CellID")

##
# Add GC bias to metadata

row.names(gcBiasDF) <- gsub("_e.*", "", row.names(gcBiasDF))
gcBiasDF$CellID <- gsub("_e.*", "", row.names(gcBiasDF))
gcBiasDF <- gcBiasDF[2:(nrow(gcBiasDF)-1), ]
metDatDF <- merge(x = metDatDF, y = gcBiasDF, by.x = "CellID", by.y = "CellID")
