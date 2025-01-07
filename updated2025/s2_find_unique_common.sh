#!/bin/bash

# Define input BED files
files=("112cell_30_31.bed" "64cell_28_29.bed" "ef1_common_cinte.bed" "lateG_26_27.bed" "midN1_25.bed")
#files=("ef1_common_general.bed" "midN1_25_general.bed" "64cell_28_29_general.bed" "lateG_26_27_general.bed" "112cell_30_31_general.bed")

# Ensure BEDTools is installed
if ! command -v bedtools &> /dev/null; then
    echo "Error: BEDTools is not installed."
    exit 1
fi

# Sort BED files
for file in "${files[@]}"; do
    sort -k1,1 -k2,2n "$file" > "${file%.bed}_sorted.bed"
done

# Unique peaks for each file
for file in "${files[@]}"; do
    sorted_file="${file%.bed}_sorted.bed"
    others=$(ls *_sorted.bed | grep -v "$sorted_file" | tr '\n' ' ')
    bedtools intersect -v -a "$sorted_file" -b $others > "${file%.bed}_unique.bed"
done

# Find common peaks across all files
bedtools intersect -a "${files[0]%.bed}_sorted.bed" \
    -b "${files[1]%.bed}_sorted.bed" \
    -b "${files[2]%.bed}_sorted.bed" \
    -b "${files[3]%.bed}_sorted.bed" \
    -b "${files[4]%.bed}_sorted.bed" > common_peaks.bed

# Cleanup sorted files (optional)
rm *_sorted.bed

cat *unique.bed > unique_peaks.bed
