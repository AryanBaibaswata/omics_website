#!/bin/bash

# Check for required arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input.bam> <output.bed>"
    exit 1
fi

INPUT_BAM="$1"
OUTPUT_BED="$2"

# Validate that the input BAM file exists
if [ ! -f "$INPUT_BAM" ]; then
    echo "Error: Input BAM file '$INPUT_BAM' not found."
    exit 1
fi

# Create a temporary file to store genome size information
GENOME_FILE=$(mktemp)

# Get the reference genome sizes from the BAM file header
samtools view -H "$INPUT_BAM" | grep "@SQ" | awk -F'\t' '{
    for(i=2; i<=NF; i++) {
        if($i ~ /^SN:/) { chr=substr($i, 4) }
        if($i ~ /^LN:/) { len=substr($i, 4) }
    }
    if(chr && len) { print chr"\t"len }
}' > "$GENOME_FILE"

# Extract aligned regions and convert them to BED format
samtools view -F 4 "$INPUT_BAM" | awk '{
    chr=$3; start=$4-1; end=start+length($10);
    if(chr != "*" && start >= 0) {
        print chr"\t"start"\t"end
    }
}' | sort -k1,1 -k2,2n | bedtools merge > aligned_regions.bed

# Identify unaligned regions
bedtools complement -i aligned_regions.bed -g "$GENOME_FILE" > "$OUTPUT_BED"

# Clean up temporary files
rm -f aligned_regions.bed "$GENOME_FILE"

echo "Unaligned regions have been saved to '$OUTPUT_BED'."

