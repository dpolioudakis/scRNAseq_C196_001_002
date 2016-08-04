## Code to run data normalization for single cell RNA-seq data
## Inputs: Matrix of raw counts for each gene for each cell
##         Metadata file containing experimental covariates  

## Load libraries
library(xlsx);
library(genefilter);
library(statmod);
library(scLVM);

## Data directory
datadir = "/geschwindlabshares/RNAseq_singlecellfetal/C196-001_002_GR/normalization/"
## Data from 192 single cells
Run1data = paste0(datadir,"Exprs_HTSCexon_Run1.csv")
## Output directory for figures
outputdir = "/geschwindlabshares/RNAseq_singlecellfetal/C196-001_002_GR/normalization/separateERCC_CapPlate/figs/"
## Output directory for data
outputdirData = "/geschwindlabshares/RNAseq_singlecellfetal/C196-001_002_GR/normalization/separateERCC_CapPlate/"
## Metadata
metadata = "/geschwindlabshares/RNAseq_singlecellfetal/C196-001_JS/Metadata/CellDatabase_wBarcodes.xlsx"
## Minimum sequencing coverage to remove unusable cells
mincov = 500000

#######################################################################

## Read metadata
sampleInfo = read.xlsx(metadata,1)

## Filter metadata for cells meeting minimum coverage
sampleInfoFiltered = sampleInfo[sampleInfo$NumReads>mincov,]

## Read count data
run1Counts = read.table(Run1data, sep=",", header=TRUE, row.names=1)

## Make sure that cell names are identical between metadata and count data
names(run1Counts) <- sub("_.+","",names(run1Counts),perl=TRUE)
sampleInfoFiltered$CellID <- paste0("Cell",sampleInfoFiltered$CellID)
sampleInfoFiltered_CP1 = sampleInfoFiltered[sampleInfoFiltered$CapturePlateID=="CP1",]
sampleInfoFiltered_CP2 = sampleInfoFiltered[sampleInfoFiltered$CapturePlateID=="CP2",]

## Filter count data to only cells meeting minimum coverage
run1CountsFiltered <- run1Counts[,sampleInfoFiltered$CellID]

## TPM normalize
run1CountsFilteredTPM <- sweep(run1CountsFiltered,2,sampleInfoFiltered$NumReads/1000000,'/')
run1CountsFilteredTPM_CP1 <- run1CountsFilteredTPM[,sampleInfoFiltered_CP1$CellID]
run1CountsFilteredTPM_CP2 <- run1CountsFilteredTPM[,sampleInfoFiltered_CP2$CellID]

## ERCC normalize (separate by Capture Plate)
run1CountsFilteredTPMERCC <- run1CountsFilteredTPM[grepl("^ERCC",row.names(run1CountsFilteredTPM),perl=TRUE),]
run1CountsFilteredTPMERCC_CP1 <- run1CountsFilteredTPMERCC[,sampleInfoFiltered_CP1$CellID]
run1CountsFilteredTPMERCC_CP2 <- run1CountsFilteredTPMERCC[,sampleInfoFiltered_CP2$CellID]
# Calculate Scaling Factors
ERCCScalingFactor_CP1 <- estimateSizeFactorsForMatrix(run1CountsFilteredTPMERCC_CP1);
ERCCScalingFactor_CP2 <- estimateSizeFactorsForMatrix(run1CountsFilteredTPMERCC_CP2);
ERCCScalingFactor <- c(ERCCScalingFactor_CP1,ERCCScalingFactor_CP2)
run1CountsFilteredTPMNorm_CP1 <- sweep(run1CountsFilteredTPM_CP1, 2, ERCCScalingFactor_CP1, '/')
run1CountsFilteredTPMNorm_CP2 <- sweep(run1CountsFilteredTPM_CP2, 2, ERCCScalingFactor_CP2, '/')
run1CountsFilteredTPMNorm <- cbind(run1CountsFilteredTPMNorm_CP1,run1CountsFilteredTPMNorm_CP2)