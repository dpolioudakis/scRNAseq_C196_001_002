#!/bin/bash

# Damon Polioudakis
# 2016-04-21
# Select fastqs for human only capture sites and merge lane 1 and 2 to prepare
# for variant calling
# Suggested script calls:
#   qsub -cwd -o logs/Merge_FASTQs_Variant_Calling_QSUB_$(date +%Y%m%d).log -e logs/Merge_FASTQs_Variant_Calling_QSUB_$(date +%Y%m%d).error -S /bin/bash -V -N MergeFQ -q geschwind.q -l h_data=4G,h_rt=12:00:00 Merge_FASTQs_Variant_Calling.sh
################################################################################

echo ""
echo "Starting Merge_FASTQs_Variant_Calling.sh"
echo ""
################################################################################

# Define Input Variables and Functions

inSampleID=../analysis/tables/Human_Only_Capture_Sites_10^5Hs_10^5Mm.txt
outDir=../data/fastq/Merged_For_Variant_Calling
mkdir -p ${outDir}
################################################################################

# Merge

ls ../data/fastq/SxaQSEQsXap089L2/*fastq.gz | while read pathFastq; do
  cat "${pathFastq}" ../data/fastq/SxaQSEQsXbp060L2/$(basename ${pathFastq}) > ${outDir}/$(basename ${pathFastq})
done
################################################################################

echo ""
echo "End of Merge_FASTQs_Variant_Calling.sh... "$(date)
