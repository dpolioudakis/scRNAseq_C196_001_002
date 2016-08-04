#!/bin/bash

# Damon Polioudakis
# 2016-02-09
# Run HTSeq Counts
################################################################################
echo ""
echo "Starting 4_HTSC_Expression.sh... "$(date)
################################################################################

inBam=$1
sampleID=$2
outDir=$3
samDir=$4

gtf=../source/gencode.v19.annotation.gtf
outHTSC=${outDir}/${sampleID}_exon_union_count.txt
# Samtools sort output path, samtools adds .bam to ouput so leave suffix off
# file name
bamNameSrtdPfx=${samDir}/PEmatched_deduplicated_name_sorted
# Samtools sort output path
bamNameSrtd=${samDir}/PEmatched_deduplicated_name_sorted.bam
# If Samtools view conversion of .bam to .sam pipe into htseq-count does not
# work, can write .sam to disk and then input into htseq-count
# samPath=${samDir}/PEmatched_deduplicated_name_sorted.sam
################################################################################

echo ""
echo "In bam:"
echo ${inBam}
echo "SampleID:"
echo ${sampleID}
echo "Out dir:"
echo ${outDir}

if [ ! -s outHTSC ]; then
  # Samtools sort by bam by name (required by HTSeq)
  echo "Start samtools sorting by name on...:"
  echo ${inBam}
  # samtools adds .bam to ouput so left off .bam suffix from variable
  samtools sort -n ${inBam} ${bamNameSrtdPfx}
  echo "Done samtools sorting by name, output:"
  echo ${bamNameSrtd}

  # Convert name sorted bam to sam and pipe to HTseq
  echo "Started HTSC exon on...:"
  echo ${bamNameSrtd}
  # "-" tells HTseq to input from pipe
  samtools view -h ${bamNameSrtd} | /share/apps/anaconda/bin/htseq-count --stranded=no --mode=union --type=exon - $gtf >> ${outHTSC}
  rm -f ${bamNameSrtd}
  echo "Saving HTSC exon results to...: "${outHTSC}
else
  echo "${outHTSC} already exists"
fi
################################################################################

echo ""
echo "End of 4_HTSC_Expression.sh... "$(date)
