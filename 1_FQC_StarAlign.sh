#!/bin/bash

# Damon Polioudakis
# 2015-01-11
# Code to perform FastQC and align a sample with STAR
# This script is meant to be fed into 2_FQC_StarAlign_QSUB.sh

### Add paths and variables (if needed)

export PATH=$PATH:/share/apps/STAR_2.4.0j/bin/Linux_x86_64

# Star index is Jason Stein's
STAR_INDEX=/geschwindlabshares/eQTL/GenomeBuild/IndexedGenome

# Assign arguments to variables

outDir=$1
inFqR1=$2
inFqR2=$3
inDir=$4
sampleID=$5

### Run FASTQC

echo ""
echo $sampleID
echo "Script is beginning..."

# Make a "fastqc" directory for all fastqc output - fastqc output for every
# sample will be output here
mkdir -p $inDir"/fastqc"

if [ ! -f "${inDir}/fastqc/${sampleID}_R2_fastqc.zip" ]; then
	/share/apps/FastQC-v0.11.2/fastqc --noextract --outdir $inDir"/fastqc" $inFqR1
	/share/apps/FastQC-v0.11.2/fastqc --noextract --outdir $inDir"/fastqc" $inFqR2
	echo "FastQC complete"
else
	echo "FastQC already performed for this sample lane"
fi

### STAR Alignment

## STAR and .SAM to .BAM

mkdir -p $outDir

if [ ! -f Aligned.sortedByCoord.out.bam ]; then
	# Jason's settings:
	STAR --runMode alignReads --genomeDir $STAR_INDEX --outFileNamePrefix $outDir/ --readFilesCommand zcat --readFilesIn $inFqR1 $inFqR2 --runThreadN 8 --outSAMtype BAM SortedByCoordinate
	# Jill's settings:
	# STAR --runMode alignReads --genomeDir $STAR_INDEX --readFilesCommand zcat --readFilesIn $inFqR1 $inFqR2 --runThreadN 8 --alignMatesGapMax 199  --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical
		# --outFilterIntronMotifs RemoveNoncanonical: filter out alignments that contain non-canonical junctions
		# outSAMstrandField intronMotif: strand derived from the intron motif. Reads with inconsistent and/or non-canonical introns are filtered out.
	echo "STAR alignment complete"
	# Jill's settings
	# samtools view -bS $outDir/Aligned.out.sam > $outDir/Aligned.out.bam
	# echo ".SAM now converted to .BAM"
	# rm -f $outDir/Aligned.out.sam
else
	echo "STAR alignment already complete for this sample lane"
fi
