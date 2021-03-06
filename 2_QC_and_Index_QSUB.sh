#!/bin/bash

# Damon Polioudakis
# 2016-01-28
# Code to submit a group of samples - run in either one lane OR multiple lanes - for Stage 2 of RNASeq Pipeline (QC)
# To submit this script: qsub -cwd -o logs/2_QC_and_Index_QSUB_$(date +%Y%m%d).log -e logs/2_QC_and_Index_QSUB_$(date +%Y%m%d).error -S /bin/bash -V -N Submit_2 -q geschwind.q -l h_data=4G,h_rt=1:00:00 2_QC_and_Index_QSUB.sh

inParentDir=../data/bam
inLaneDir=${inParentDir}/SxaQSEQsXap089L2
outParentDir=../data/bam/merged

# Loop through sample directories inside each lane directory
for sampleDir in ${inLaneDir}/*; do

	echo ""
	echo "Sample dir:"
	echo ${sampleDir}

	# inParentDir are all lanes for each sample
	bnSampleDir=$(basename ${sampleDir})
	echo "In parent dir"
	echo ${inParentDir}

	outDir=${outParentDir}/${bnSampleDir}
	echo "Out dir:"
	echo ${outDir}

	# To check if final file output from QC_and_Index.sh script exists
	outFpath=${outDir}/PEmatched_markdup_sorted.bai

	# If final file output from QC_and_Index.sh script does not exist for that
	# sample, then qsub QC_and_Index.sh script for that sample
	if [ ! -f ${outFpath} ]; then

		echo ${outFpath}" does not exist, qsubbing ${inParentDir} to QC_and_Index.sh"

		# Log names for QC_and_Index.sh script
		outLog=logs/2_QC_and_Index_$(date +%Y%m%d)_${bnSampleDir}.log
		outErr=logs/2_QC_and_Index_$(date +%Y%m%d)_${bnSampleDir}.error

		# Note name of job (-N) will only display 1st 10 characters
		qsub -cwd -o ${outLog} -e ${outErr} -V -S /bin/bash -q geschwind.q -N QC_${bnSampleDir} -l h_data=36G,h_rt=24:00:00 2_QC_and_Index.sh ${inParentDir} ${outDir} ${bnSampleDir}

	else
		echo "${outFpath} already exists..."
	fi

done
