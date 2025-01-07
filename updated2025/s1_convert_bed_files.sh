#!/bin/bash

# Define the reference file and output directory
reference=~/Desktop/BNKA01_contigs_conversion_final.txt
output_dir=converted

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Create a temporary file with the mapping for `awk`
awk 'NR > 1 {print $2, $1}' "$reference" > chrom_map.tmp

# Process each BED file in the current directory
for bedfile in *.bed; do
    # Skip if no BED files are found
    [ -e "$bedfile" ] || continue

    # Generate the output file path
    output_file="$output_dir/$bedfile"

    # Convert the chromosome names in the BED file
    awk -v OFS="\t" '
    BEGIN {
        while ((getline < "chrom_map.tmp") > 0) map[$1] = $2;
        close("chrom_map.tmp");
    }
    {
        if ($1 in map) $1 = map[$1];
        print
    }' "$bedfile" > "$output_file"

    echo "Processed $bedfile -> $output_file"
done

# Clean up the temporary file
rm -f chrom_map.tmp

echo "All files have been converted and saved in $output_dir."

