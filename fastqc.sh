#!/bin/bash

# Check if the destination folder argument is provided
if [ $# -lt 2 ]; then
    echo "Usage: $0 <fastq_files> <destination_folder>"
    exit 1
fi

# Loop through all fastq.gz files passed as arguments
for fastq_file in "${@:1:$#-1}"
do
    sample="${fastq_file%.fastq.gz}"
    echo "Processing $sample"
    # Run FastQC on each fastq.gz file
    fastqc "$fastq_file"
done

# Run MultiQC to generate a report for all fastq.gz files and save it in the specified destination folder
destination_folder="${@:$#}"
multiqc . -o "multiqc_reports/$destination_folder"


