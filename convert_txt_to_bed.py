import argparse
import subprocess

def convert_bam_to_bed(sample_id, input_txt, output_bed):
    # Step 1: Generate coverage information
    coverage_file = input_txt
    #with open(coverage_file, 'w') as coverage_out:
        #subprocess.run(["samtools", "depth", input_txt], stdout=coverage_out)

    # Step 2: Convert coverage information to BED format
    with open(coverage_file, 'r') as infile, open(output_bed, 'w') as outfile:
        for line in infile:
            fields = line.strip().split()
            chrom = fields[0]
            pos = int(fields[1]) - 1  # BED format is 0-based
            coverage = fields[2]
            # Create a BED interval for each position with the coverage as the score
            outfile.write(f"{chrom}\t{pos}\t{pos + 1}\t{coverage}\n")
    print(f"${coverage_file} ${outfile}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert BAM to BED using samtools depth")
    parser.add_argument("sample_id", help="Sample ID")
    parser.add_argument("input_bam", help="Input BAM file")
    parser.add_argument("output_bed", help="Output BED file")
    args = parser.parse_args()

convert_bam_to_bed(args.sample_id, args.input_bam, args.output_bed)