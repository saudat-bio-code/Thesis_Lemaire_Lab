library(BSgenome.Crobusta.HT.KY) 
library(BSgenome.Cintestinalis)
source("~/Desktop/scripts/Thesis_Lemaire_Lab/scripts/get_nulll_seq_edited_2025.R")

#setwd("/Users/saudat/Desktop/bed_files/final_bed_files/cinte_files/converted/general_peaks/genic_regions/")
setwd("/Users/saudat/Desktop/bed_files/final_bed_files/cinte_files/converted/general_peaks/genic_regions//")

trf = import("~/Desktop/genomes/Ros.fasta.2.5.7.80.10.50.500.bed")
genome =  BSgenome.Cintestinalis
#bed_files <- list.files(pattern = "\\cinte_FilteredGenicRegions.bed$")
bed_files <- list.files(pattern = "\\cinte_FiltedEnhancer.bed$")


# Loop through each BED file and apply the function
for (bed_file in bed_files) {
  # Extract ID (filename without extension)
  ID <- tools::file_path_sans_ext(basename(bed_file))
  # Construct output filenames
  #outputBedFN <- paste0("negSet_", ID, ".bed")
  outputPosFastaFN <- paste0("posSet_", ID, ".fa")
  outputNegFastaFN <- paste0("negSet_", ID, ".fa")
  #bed_file_gr = import(bed_file)
  # Apply the function
  genNullSeqs_ranges(inputBedFN = bed_file,
                     genomeVersion = 'no',
                     #outputBedFN = outputBedFN,
                     outputPosFastaFN = outputPosFastaFN,
                     outputNegFastaFN = outputNegFastaFN,
                     xfold = 1,
                     repeat_match_tol = 0.02,
                     GC_match_tol = 0.02,
                     length_match_tol = 0.02,
                     batchsize = 5000,
                     nMaxTrials = 20,
                     genome = genome
  )
  
  # Print progress
  cat("Processed:", bed_file, "\n")
}
# Print completion message
cat("All BED files processed!\n")



trf = import("~/Desktop/genomes/HT.Ref.fasta.2.5.7.80.10.50.500.bed")
genome =  BSgenome.Crobusta.HT.KY
#bed_files <- list.files(pattern = "\\crobu_FilteredGenicRegions.bed$")
bed_files <- list.files(pattern = "\\crobu_FiltedEnhancer.bed$")

# Loop through each BED file and apply the function
for (bed_file in bed_files) {
  # Extract ID (filename without extension)
  ID <- tools::file_path_sans_ext(basename(bed_file))
  # Construct output filenames
  #outputBedFN <- paste0("negSet_", ID, ".bed")
  outputPosFastaFN <- paste0("posSet_", ID, ".fa")
  outputNegFastaFN <- paste0("negSet_", ID, ".fa")
  #bed_file_gr = import(bed_file)
  # Apply the function
  genNullSeqs_ranges(inputBedFN = bed_file,
                     genomeVersion = 'no',
                     #outputBedFN = outputBedFN,
                     outputPosFastaFN = outputPosFastaFN,
                     outputNegFastaFN = outputNegFastaFN,
                     xfold = 1,
                     repeat_match_tol = 0.02,
                     GC_match_tol = 0.02,
                     length_match_tol = 0.02,
                     batchsize = 5000,
                     nMaxTrials = 20,
                     genome = genome
  )
  
  # Print progress
  cat("Processed:", bed_file, "\n")
}
# Print completion message
cat("All BED files processed!\n")
