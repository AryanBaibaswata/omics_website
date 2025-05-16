#!/bin/bash
set -e
basedir="C:\omics_website\uploads"
samples=(
    "SRR15991306"
)
GENOMEIDX1="/mnt/c/Users/asus/genome_assembly/genome_sars_cov2/NC_045512.2.fasta"

for sample_name in "${samples[@]}"; do
    echo "Step-1.0: FastQC Quality Control Report for ${sample_name}"
    fastqc -o "${basedir}/fastqc_output/" "${basedir}/${sample_name}_R1.fastq.gz" "${basedir}/${sample_name}_R2.fastq.gz"

    echo "Step-1.1: Fastp Quality Control for ${sample_name}"
    fastp -i "${basedir}/${sample_name}_R1.fastq.gz" -o "${basedir}/${sample_name}_P1.fastq"           -I "${basedir}/${sample_name}_R2.fastq.gz" -O "${basedir}/${sample_name}_P2.fastq"           --thread 4 -h "${basedir}/fastp-${sample_name}.html" 2> "${basedir}/fastp-${sample_name}.log"

    echo "Step-2: Read Alignment for ${sample_name}"
    bowtie2 -p 64 -x "${GENOMEIDX1}" -1 "${basedir}/${sample_name}_P1.fastq" -2 "${basedir}/${sample_name}_P2.fastq" -S "${basedir}/${sample_name}.sam"

    echo "Step-3: Conversion Of Sam To BAM File for ${sample_name}"
    samtools view -b "${basedir}/${sample_name}.sam" -o "${basedir}/${sample_name}.bam"

    echo "Step-4: Alignment Metrics for ${sample_name}"
    samtools flagstat "${basedir}/${sample_name}.bam" > "${basedir}/${sample_name}.flagstat.txt"

    echo "Step-5: Conversion of BAM To Sorted BAM for ${sample_name}"
    samtools sort "${basedir}/${sample_name}.bam" -o "${basedir}/${sample_name}.sorted.bam"

    echo "Step-9: Removing duplicate reads from Sorted Bam Files for ${sample_name}"
    samtools rmdup -S "${basedir}/${sample_name}.sorted.bam" "${basedir}/${sample_name}.duprem.bam"

    echo "Step-6: Deriving Low Coverage Bed File for ${sample_name}"
    samtools depth "${basedir}/${sample_name}.sorted.bam" | awk '$3 < 5 {print $1"	"$2"	"$3}' > "${basedir}/coverage_${sample_name}.txt"

    echo "Step-7: Extracting start end coordinates of missing read segments for ${sample_name}"
    python3 "${basedir}/extract_intervalsforBED.py" "${basedir}/coverage_${sample_name}.txt" > "${basedir}/output_${sample_name}.bed"

    echo "Step-8: Performing N-masking for ${sample_name}"
    bedtools maskfasta -fi "${GENOMEIDX1}" -bed "${basedir}/output_${sample_name}.bed" -mc N -fo "${basedir}/${sample_name}_masked.fasta"

 

    echo "Step-10: Generation of VCF for ${sample_name}"
    bcftools mpileup -d 1000000 -f "${GENOMEIDX1}" "${basedir}/${sample_name}.duprem.bam" | bcftools call -cv --ploidy 1 -Oz -o "${basedir}/${sample_name}.vcf.gz"

    echo "Step-11: Generation of VCF Index for ${sample_name}"
    bcftools index "${basedir}/${sample_name}.vcf.gz"

    echo "Step-12: Generation of viral genome fasta for ${sample_name}"
    cat "${basedir}/${sample_name}_masked.fasta" | bcftools consensus "${basedir}/${sample_name}.vcf.gz" > "${basedir}/${sample_name}_genome.fa"
done
