#!/bin/bash

# Set target and reference genome names
#TARGET_GENOME="Ciona_intestinalis"  # Edit this as needed
#REFERENCE_GENOME="Ciona_robusta"    # Edit this as needed
TARGET_GENOME="Ciona_robusta"
REFERENCE_GENOME="Ciona_intestinalis"
TARGET_FASTA="/Users/saudat/Desktop/HT/HT.Ref.fasta"
TARGET_SIZE="./Ros.genome.size"
# I need to also adjust getfasta to automatically chose the right reference file

#/Users/saudat/Desktop/Ros/Ros.fasta
#/Users/saudat/Desktop/HT/
#/Users/saudat/Desktop/HT/HT.Ref.fasta

#halAlignmentDepth --outWiggle alignmentDepthFile /Users/saudat/Desktop/svm_all/rep_crobu/ancestor_ciona_newick_final_ros ${REFERENCE_GENOME} 


# Check for required inputs
if [[ -z "$TARGET_GENOME" || -z "$REFERENCE_GENOME" ]]; then
  echo "Error: Please define both TARGET_GENOME and REFERENCE_GENOME variables."
  exit 1
fi

# Create an output directory if it doesn't exist
OUTPUT_DIR="processed_files"
mkdir -p "$OUTPUT_DIR"

# Loop through all input BED files
for FILE in *.bed; do
  # Extract the base filename (without extension)
  BASENAME=$(basename "$FILE" .bed)

  # Define intermediate and output filenames
  HIGHSCORES_FILE="${BASENAME}_HScores.bed"
  LIFTED_HIGHSCORES_FILE="${BASENAME}_HScores_lifted.bed"
  LIFTED_BED_FILE="${BASENAME}_${TARGET_GENOME}.bed"
  ORTHOLOG_VALIDATED_FILE="${BASENAME}_${TARGET_GENOME}_validated.bed"
  FINAL_BED_FILE="${BASENAME}_${TARGET_GENOME}_final.bed"
  FINAL_FASTA_FILE="${OUTPUT_DIR}/${BASENAME}_${TARGET_GENOME}.fa"

  # Step 1: Get max score positions
  python3 /Users/saudat/halLiftover-postprocessing/getMaxScorePositionFromWig.py \
    --bedFileName "$FILE" \
    --wigFileName alignmentDepthFile \
    --chromSizesFileName "$TARGET_SIZE" \
    --highestScoreLocationFileName "$HIGHSCORES_FILE" \
    --tempDir=/tmp

  # Step 2: Lift max positions
  halLiftover /Users/saudat/Desktop/svm_all/rep_crobu/ancestor_ciona_newick_final_ros \
    "$REFERENCE_GENOME" "$HIGHSCORES_FILE" \
    "$TARGET_GENOME" "$LIFTED_HIGHSCORES_FILE"

  # Step 3: Lift the full BED file
  halLiftover /Users/saudat/Desktop/svm_all/rep_crobu/ancestor_ciona_newick_final_ros \
    "$REFERENCE_GENOME" "$FILE" \
    "$TARGET_GENOME" "$LIFTED_BED_FILE"

  # Step 4: Verify orthologous positions
  python /Users/saudat/halLiftover-postprocessing/orthologFind.py \
    -min_len 50 -max_len 3000 -protect_dist 5 \
    -qFile "$FILE" \
    -tFile "$LIFTED_BED_FILE" \
    -sFile "$LIFTED_HIGHSCORES_FILE" \
    -oFile "$ORTHOLOG_VALIDATED_FILE" \
    -mult_keepone

  # Step 5: Format to 6-column BED file
  cut -f 1,2,3,5,6 "$ORTHOLOG_VALIDATED_FILE" > "$FINAL_BED_FILE"

  # Step 7: Extract FASTA sequences
  bedtools getfasta \
    -fi "$TARGET_FASTA" \
    -bed "$FINAL_BED_FILE" \
    -fo "$FINAL_FASTA_FILE" \
    -nameOnly

  # Clean up intermediate files
  rm -f "$HIGHSCORES_FILE" "$LIFTED_HIGHSCORES_FILE" "$LIFTED_BED_FILE" "$ORTHOLOG_VALIDATED_FILE"
done

echo "Processing completed. Outputs saved in $OUTPUT_DIR"

rm *png
rm *failed
