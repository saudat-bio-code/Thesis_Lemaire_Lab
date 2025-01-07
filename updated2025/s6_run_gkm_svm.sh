#!/bin/bash

# Set the directory containing your files
input_dir="."
model_dir="$input_dir/models"
random_kmers="../../random10mers.fasta"

# Create the models directory if it doesn't exist
mkdir -p "$model_dir"

# Iterate over all paired files
# Define the input directory if needed, or set it to the current directory
input_dir="."

# Loop over all files matching 'posSet_*crobu_FiltedEnhancer.fa'
for focal_file in posSet_*cinte_FiltedEnhancer.fa; do
  base_name="${focal_file#posSet_}"  # Remove the 'posSet_' prefix
  # Extract the base name by removing the 'posSet_' prefix and 'cinte_FiltedEnhancer.fa' suffix
  base_name="${base_name%cinte_FiltedEnhancer.fa}"  # Remove the suffix

  # Construct the ancestor file name
  ancestor_file="${input_dir}/${base_name}anc_FiltedEnhancer.fa"
  orig_file="${input_dir}/${base_name}cinte_FiltedEnhancer.fa"
  echo "$ancestor_file"
  # Ensure the ancestor file exists
  if [[ -f "$ancestor_file" ]]; then
    # Define output file names
    model="${model_dir}/${base_name}cinte_MODEL"
    prediction="${base_name}_cinte.random10mers.txt"
    echo "negSet_${base_name}cinte_FiltedEnhancer.fa"
    # Run gkmtrain
    gkmtrain -l 10  "$focal_file" "negSet_${base_name}cinte_FiltedEnhancer.fa" "$model"

    # Wait for training to complete before predicting
    #wait

    # Run gkmpredict
    gkmpredict "$random_kmers" "${model}.model.txt" "$prediction"

    echo "$prediction"
    # Run testPosSelecPrint.pl
    perl /Users/saudat/Desktop/bed_files/testinggkm/testPosSelec.pl "$orig_file" "$ancestor_file" "$prediction" 5000 highertail "$model"

    echo "Completed analysis for $base_name"
  else
    echo "Ancestor file not found for $focal_file. Skipping."
  fi
done



# Loop over all files matching 'posSet_*crobu_FiltedEnhancer.fa'
for focal_file in *crobu_FiltedEnhancer.fa; do
  base_name="${focal_file#posSet_}"  # Remove the 'posSet_' prefix
  # Extract the base name by removing the 'posSet_' prefix and 'crobu_FiltedEnhancer.fa' suffix
  base_name="${base_name%crobu_FiltedEnhancer.fa}"  # Remove the 'crobu_FiltedEnhancer.fa' suffix

  # Construct the ancestor file name
  ancestor_file="${input_dir}/${base_name}anc_FiltedEnhancer.fa"
  orig_file="${input_dir}/${base_name}crobu_FiltedEnhancer.fa"

  # Ensure the ancestor file exists
  if [[ -f "$ancestor_file" ]]; then
    # Define output file names
    model="${model_dir}/${base_name}crobu_MODEL"
    prediction="${base_name}_crobu.random10mers.txt"

    # Run gkmtrain
    gkmtrain -l 10  "$focal_file" "negSet_${base_name}crobu_FiltedEnhancer.fa" "$model"

    # Wait for training to complete before predicting
    #wait

    # Run gkmpredict
    gkmpredict "$random_kmers" "${model}.model.txt" "$prediction"

    echo "$prediction"
    # Run testPosSelecPrint.pl
    perl /Users/saudat/Desktop/bed_files/testinggkm/testPosSelec.pl "$orig_file" "$ancestor_file" "$prediction" 5000 highertail "$model"

    echo "Completed analysis for $base_name"
  else
    echo "Ancestor file not found for $focal_file. Skipping."
  fi
done
