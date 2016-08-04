#!/bin/bash

# Damon Polioudakis
# 2016-03-16
# Count and compile mapped reads from VZ and CP from Jason RNAseq data for ATAC # project

# To submit this script:
# qsub -cwd -o logs/Compile_Mapped_Reads_Bulk_$(date +%Y%m%d).log -e logs/Compile_Mapped_Reads_Bulk_$(date +%Y%m%d).error -S /bin/bash -V -N Count_Map -q geschwind.q -l h_data=4G,h_rt=3:00:00 Compile_Mapped_Reads_Bulk.sh
################################################################################
echo "Starting Compile_Mapped_Reads_Bulk.sh... "$(date)

################################################################################

# Define Input Variables and Functions

outDir=../data/metadata
outFile=${outDir}/Bulk_Mapped_Reads.txt

################################################################################

mkdir -p $outDir

echo -e "SampleID\tNumber_Mapped_Reads" > ${outFile}

for inBam in $(awk 'FS="," {print $20}' < ../metadata/VZCP_sampleinfo.csv); do

  echo "In bam: "${inBam}
  # sampleID is the directory named by sample ID (removing
  # "../SxaQSEQsXbp060L2_bam/") from the path
  sampleID=$(basename ${inBam} .merge.remarkdup.bam)
	echo "In bam sample name: "${sampleID}

  echo "Starting counting mapped reads..."
  count=$(samtools view -F 0x904 ${inBam} | cut -f 1 | wc -l)
  echo "Done counting mapped reads..."

  echo -e ${sampleID}"\t"${count}
  echo -e ${sampleID}"\t"${count} >> ${outFile}

done
################################################################################

date +%F%t%T
echo "End of Compile_Mapped_Reads_Bulk.sh..."
