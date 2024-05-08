for i in *.fastq.gz
do
{
sample="${i%.fastq.gz}"
echo $sample

fastqc ${sample}

}
done

multiqc .

#add html generation steps here

echo "<!DOCTYPE html>
<html>
<head>
    <title>Fastq Report</title>
</head>
<body>
    <h1>Report for ${multiqc_report.html}</h1>
    <p>Processing completed successfully.</p>
</body>
</html>" > report.html