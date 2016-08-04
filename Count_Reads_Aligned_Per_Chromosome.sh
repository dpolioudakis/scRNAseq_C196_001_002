#!/bin/bash

#####TODO: Format the output of uniq -c, and remove the | grep

# Damon Polioudakis
# 2016-05-26
# Count reads aligned per chromosome

# Suggest script call:
# sh Count_Reads_Aligned_Per_Chromosome.sh 2>&1 | tee logs/Count_Reads_Aligned_Per_Chromosome_$(date +%Y%m%d).log
# Reminder: make /logs directory in code directory

# Outputs counts from all lanes into 1 file
################################################################################
echo ""
echo "Starting Count_Reads_Aligned_Per_Chromosome.sh..."$(date)
echo ""
################################################################################

# Define Input Variables and Functions

# Path to directory that contains bams
inParentDir=../data/bam/merged
# Path for parent directory to output processing files, will create
# subdirectories for each sample ID
outCounts=../analysis/tables/Compiled_Number_Reads_Aligned_Per_Chromosome.txt
################################################################################

> ${outCounts}

# Convert bams to sams, extract chromosome, count occurances

# -q 254 option selects reads with map quality score > 254
# The mapping quality MAPQ (column 5) is 255 for uniquely mapping reads
for inBam in ${inParentDir}/Cell*/Aligned.sortedByCoord.out.bam; do
  samtools view -q 254 ${inBam} | cut -f3 | sort -n | uniq -c | grep '.*Y.*'
  >> ${outCounts}
done




#
#
# mkdir -p ${outParentDir}
#
# # Header for output compiled statistics file
# echo -e "SampleID\tReads_In_Human\tReads_In_Mouse\tAlign_To_Both\tAlign_To_Human_Only\tAlign_To_Mouse_Only" > ${outStats}
#
# # Convert bams to sams, extract read names, and sort unique
# for inBamDir in ${inParentDir}/*/SxaQSEQsXap096L*/*; do
#
#   echo ""
#   echo "Converting bam to sam, extracting read names, and sort unique..."
#
#   inBam=${inBamDir}/Aligned.sortedByCoord.out.bam
#   echo "In Bam:"
#   echo ${inBam}
#
#   outDir=${outParentDir}/${inBamDir##${inParentDir}/}
#   echo "Out Dir:"
#   echo ${outDir}
#   mkdir -p ${outDir}
#
#   # -q 254 option selects reads with map quality score > 254
#   # The mapping quality MAPQ (column 5) is 255 for uniquely mapping reads
#   samtools view -q 254 ${inBam} | cut -f1 | sort -u > ${outDir}/Reads.txt
# done
#
# # Compare sorted read name files line by line to identify reads common to both
# # mouse and human alignments
# for outDirHuman in ${outParentDir}/NoERCC/Sxa*/*; do
#
#   echo ""
#   echo "Identifying reads common to both bams..."
#
#   echo "Out Dir Human:"
#   echo ${outDirHuman}
#
#   outDirMouse=$(echo ${outDirHuman} | sed s/NoERCC/mouse/)
#   echo "Out Dir Mouse:"
#   echo ${outDirMouse}
#
#   outDir=${outParentDir}/${outDirHuman##../data/compare_human_to_mouse/NoERCC/}
#   echo "Out Dir:"
#   echo ${outDir}
#   mkdir -p ${outDir}
#
#   sampleID=$(basename ${outDirHuman})
#   echo "Sample ID:"
#   echo ${sampleID}
#
#   # Reads Aligned to Both
#   comm -12 ${outDirHuman}/Reads.txt ${outDirMouse}/Reads.txt > ${outDir}/Align_To_Both.txt
#
#   # Reads Aligned to Human Only
#   comm -23 ${outDirHuman}/Reads.txt ${outDirMouse}/Reads.txt > ${outDir}/Align_To_Human_Only.txt
#
#   # Reads Aligned to Mouse Only
#   comm -13 ${outDirHuman}/Reads.txt ${outDirMouse}/Reads.txt > ${outDir}/Align_To_Mouse_Only.txt
#
#   # Compile and write out numbers aligning to both, human only, or mouse only
#   echo -ne ${sampleID}"\t" >> ${outStats}
#   echo -ne $(wc -l ${outDirHuman}/Reads.txt | cut -f1 -d" ")"\t" >> ${outStats}
#   echo -ne $(wc -l ${outDirMouse}/Reads.txt | cut -f1 -d" ")"\t" >> ${outStats}
#   echo -ne $(wc -l ${outDir}/Align_To_Both.txt | cut -f1 -d" ")"\t" >> ${outStats}
#   echo -ne $(wc -l ${outDir}/Align_To_Human_Only.txt | cut -f1 -d" ")"\t" >> ${outStats}
#   echo -e $(wc -l ${outDir}/Align_To_Mouse_Only.txt | cut -f1 -d" ") >> ${outStats}
#
# done
#
# # Remove intermediate files
# rm -r ${outParentDir}/mouse
# rm -r ${outParentDir}/NoERCC
# rm -r ${outParentDir}/SxaQSEQsXap096L*
################################################################################

echo ""
echo "End of Count_Reads_Aligned_Per_Chromosome.sh... "$(date)
