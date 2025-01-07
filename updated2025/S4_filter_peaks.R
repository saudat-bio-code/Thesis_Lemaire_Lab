library(data.table)
library(BSgenome)
library(Biostrings)
library(peakToGene)
library(ChIPpeakAnno)
library(rtracklayer)
library(IRanges)

# Define the directory paths
input_dir <- "/Users/saudat/Desktop/bed_files/final_bed_files/cinte_files/converted/unique_peaks/"
output_dir <- "/Users/saudat/Desktop/bed_files/final_bed_files/cinte_files/converted/unique_peaks/enhancers/"
#snps <- import( "/Users/saudat/Desktop/bed_files/final_bed_files/cinte_files/converted/s4_mutations/combined_sorted.cinte.bed")
snps <- import( "/Users/saudat/Desktop/bed_files/final_bed_files/cinte_files/converted/s4_mutations/combined_sorted.crobu.bed")

dir.create(output_dir, showWarnings = FALSE)
ky <- getFeatures('~/Downloads/HT.Gene.gff3') 
# Set the working directory containing the `.bed` files
setwd("/Users/saudat/Desktop/bed_files/final_bed_files/cinte_files/converted/unique_peaks/")
# Helper functions


cleanPeaksBasedOnSnps<-function (input_peaks, snps) {
  overlaps =IRanges::findOverlaps(input_peaks,snps)
  qc = table(table(overlaps@from)>2)
  cat(paste("number of peaks with more than or equal to 2 snps", qc[2],"\n"))
  cat(paste("number of peaks with less than or equal to 2 snps", qc[1],"\n"))
  sub_na = ifelse(table(overlaps@from)>2, table(overlaps@from),NA)
  x<-sub_na[!is.na(sub_na)]
  output = input_peaks[as.integer(names(x))]
  return(output)
}


cleanPeaksBasedOnAncestry<-function (crobu_peaks,cinte_peaks, ancestor_peaks) {
  qc = table(ancestor_peaks@elementMetadata$name  %in% cinte_peaks@elementMetadata$name & ancestor_peaks@elementMetadata$name  %in% crobu_peaks@elementMetadata$name)
  cat(paste("number of shared peaks in all files", qc[2],"\n"))
  cat(paste("number of filtered peaks that are not found in all files", qc[1],"\n"))
  common_names = Reduce(intersect, list(
    ancestor_peaks@elementMetadata$name,
    cinte_peaks@elementMetadata$name,
    crobu_peaks@elementMetadata$name
  ))
  common_names = common_names[!duplicated(common_names)]
  filtered_cinte_peaks = cinte_peaks[cinte_peaks@elementMetadata$name %in% common_names]
  filtered_crobu_peaks = crobu_peaks[crobu_peaks@elementMetadata$name %in% common_names]
  filtered_ancestor_peaks = ancestor_peaks[ancestor_peaks@elementMetadata$name %in% common_names]
  assign("common_peaks", filtered_crobu_peaks, envir = .GlobalEnv)
  assign("cinte_peaks", filtered_cinte_peaks, envir = .GlobalEnv)
  assign("crobu_peaks", filtered_crobu_peaks, envir = .GlobalEnv)
  assign("ancestor_peaks", filtered_ancestor_peaks, envir = .GlobalEnv)
  
}

# Initialize result table
results_table <- data.frame(
  File = character(),
  Step1 = integer(),
  Step2 = integer(),
  Step3 = integer(),
  Step4 = integer(),
  Step5 = integer(),
  stringsAsFactors = FALSE
)

# Process each `.bed` file in the directory
#bed_files <- list.files(pattern = "\\_cinte.bed$")
bed_files <- list.files(pattern = "\\_Ciona_robusta_final.bed$")

for (file in bed_files) {
  cat("Processing:", file, "\n")
  #core_name <- gsub("_cinte.bed", "", file)
  core_name <- gsub("_Ciona_robusta_final.bed", "", file)

  # Step 1: Import and filter for upstream enhancers and promoters
  peaks <- import(file)
  names(peaks) =peaks$name
  overlaps <- lapply(ky, getOverlaps, peaks)
  featpeaks <- sapply(overlaps, function(x) unique(x[, 1]))
  peaks <- peaks[peaks@elementMetadata$name %in% featpeaks$upstream |
                   peaks@elementMetadata$name %in% featpeaks$promoter]
  step1_count <- length(peaks)
  
  # Step 2: Discard peaks smaller than 75 bp
  peaks <- peaks[end(peaks) - start(peaks) > 75]
  step2_count <- length(peaks)
  
  # Step 3: Filter peaks with >= 2 SNPs
  # (Assumes `snps` is a GRanges object containing SNPs)
  # Replace with actual SNP data file
  peaks <- cleanPeaksBasedOnSnps(peaks, snps)
  step3_count <- length(peaks)
  
  
  # Step 4: Clean based on ancestry
  # (Assumes `ancestor_peaks`, `cinte_peaks`, and `crobu_peaks` are loaded)
  ancestor_peaks <- import(paste0(core_name, "_Anc0_final.bed"))
  #cinte_peaks <- peaks
  #crobu_peaks <- import(paste0(core_name,"_Ciona_robusta_final.bed")
  crobu_peaks <- peaks
  cinte_peaks <- import(paste0(core_name,"_cinte.bed"))
  cleanPeaksBasedOnAncestry(crobu_peaks, cinte_peaks, ancestor_peaks)
  step4_count <- length(common_peaks)
  
  # Step 5: Remove duplicated peaks
  peaks <- common_peaks[!duplicated(common_peaks@ranges) & !duplicated(common_peaks@elementMetadata$name)]
  step5_count <- length(peaks)
  
  # Export the final peaks
  #export.bed(peaks, con = paste0(output_dir, core_name, "common_FilteredGenicRegions.bed")) 
  
  peaks <- cinte_peaks[!duplicated(cinte_peaks@elementMetadata$name)] #cinte
  export.bed(peaks, con = paste0(output_dir, core_name, "cinte_FiltedEnhancer.bed"))
  
  peaks <- crobu_peaks[!duplicated(crobu_peaks@elementMetadata$name)] #crobu
  export.bed(peaks, con = paste0(output_dir, core_name, "crobu_FiltedEnhancer.bed"))
  
  peaks <- ancestor_peaks[!duplicated(ancestor_peaks@elementMetadata$name)] #crobu
  export.bed(peaks, con = paste0(output_dir, core_name, "anc_FiltedEnhancer.bed"))
  
  # Append results to the table
  results_table <- rbind(results_table, data.frame(
    File = file,
    Step1 = step1_count, #filter enhancers
    Step2 = step2_count, #75 bp enhancers
    Step3 = step3_count, #at least 2 snps
    Step4 = step4_count, #remove duplicated
    Step5 = step5_count, #ancestry
    stringsAsFactors = FALSE
  ))
}

results_table= results_table[!duplicated(results_table$File),]
# Save the results table to a CSV file
write.csv(results_table, "filtering_summary.csv", row.names = FALSE)

cat("Processing complete. Summary saved to filtering_summary.csv.\n")
