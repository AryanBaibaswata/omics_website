#!/bin/bash


# Loop through all fastq.gz files passed as arguments
for fastq_file in "$@"
do
    sample="${fastq_file%.fastq.gz}"
    echo "Processing $sample"
    # Run FastQC on each fastq.gz file
    fastqc "$fastq_file"
done

# Run MultiQC to generate a report for all fastq.gz files and save it in the multiqc_reports folder
multiqc . -o multiqc_reports

# Generate the HTML report
