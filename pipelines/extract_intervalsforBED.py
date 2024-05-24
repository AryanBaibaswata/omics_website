input_file = "coverage.txt"
output_file = "output.bed"

with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        fields = line.strip().split()
        chrom = fields[0]
        pos = int(fields[1]) - 1  # BED format is 0-based
        coverage = fields[2]
        # Create a BED interval for each position with the coverage as the score
        outfile.write(f"{chrom}\t{pos}\t{pos + 1}\t{coverage}\n")
