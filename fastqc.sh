#!/bin/bash

# Loop through all fastq.gz files passed as arguments
for fastq_file in "$@"
do
    sample="${fastq_file%.fastq.gz}"
    echo "Processing $sample"
    # Run FastQC on each fastq.gz file
    fastqc "$fastq_file"
done

# Run MultiQC to generate a report for all fastq.gz files
multiqc .

# Generate the HTML report
echo '<!DOCTYPE html>
<html>
<head>
    <title>Fastq Report</title>
</head>
<body>
    <h1>Report for MultiQC</h1>
    <p>Processing completed successfully.</p>
</body>
</html>' > report.html
