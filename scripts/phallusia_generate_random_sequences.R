library(BSgenome)
library(gkmSVM)
library(Biostrings)
library(peakToGene)
#library(BSgenome.Ciinte.JDB.Roscoff)
library(BSgenome.Phmammilata)
library(data.table)

#forgeBSgenomeDataPkg("~/seed") 
#inputBedFN = peaks in bed format, different for each stage?
#try with an old genome?
genome =  BSgenome.Phmammilata

setwd("~/Downloads/bed_files/Pmamilata/ssh_imported_final/")
trf = import("Phmamm_MTP2014_genome_2018.fasta.2.7.7.80.10.50.500.bed")
cell112_30_31_uni.bed =import("Pm02Pm09_WT_112_peaks2_unique.bed")
cell64_28_29_uni.bed<-import("Pm12Pm13_WT_64_peaks2_unique.bed")
midN1_25_uni.bed<-import("Pm14-MidNeurula_peaks_unique.bed")#these two got problems
lateG_26_27_uni.bed<-import("Pm01-LateGastrula_peaks_unique.bed")

#import negative set
cell112_30_31_uni.bed =import("negSet_112cell.bed")


hist(width(lateG_26_27_uni.bed), breaks = 1000,xlim=c(0,600), 
     main = "Positive set", xlab = "Lengths of the peaks/regulatory element")

hist(width(trf), breaks = 1000,xlim=c(0,1000), 
     main = "Positive set", xlab = "Lengths of the peaks/regulatory element")

setwd("~")
trf = import("HT.Ref.fasta.2.7.7.80.10.50.500.bed")
trf = import("2021-Roscoff.fna.2.7.7.80.10.50.500.bed")

#cell112_30_31_uni.bed= cell112_30_31_uni.bed[end(cell112_30_31_uni.bed) -start(cell112_30_31_uni.bed) > 50] 
#cell64_28_29_uni.bed= cell64_28_29_uni.bed[end(cell64_28_29_uni.bed) -start(cell64_28_29_uni.bed) > 50] 
#midN1_25_uni.bed= midN1_25_uni.bed[end(midN1_25_uni.bed) -start(midN1_25_uni.bed) > 50] 
#lateG_26_27_uni.bed= lateG_26_27_uni.bed[end(lateG_26_27_uni.bed) -start(lateG_26_27_uni.bed) > 50] 

#export.bed(cell112_30_31_uni.bed,con='112cell_30_31_uni.bed')
#export.bed(cell64_28_29_uni.bed,con='64cell_28_29_uni.bed')
#export.bed(midN1_25_uni.bed,con='midN1_25_uni.bed')
#export.bed(lateG_26_27_uni.bed,con='lateG_26_27_uni.bed')


genNullSeqs_ranges(
  "Pm02Pm09_WT_112_peaks2_unique.bed", 
  genomeVersion='no', 
  outputBedFN = 'negSet_112cell.bed', 
  outputPosFastaFN = 'posSet_112cell.fa',
  outputNegFastaFN = 'negSet_112cell.fa', 
  xfold = 1, 
  repeat_match_tol = 0.02, 
  GC_match_tol = 0.02, 
  length_match_tol = 0.02, 
  batchsize = 5000, 
  nMaxTrials = 20, 
  genome =  genome)

genNullSeqs_ranges(
  "Pm12Pm13_WT_64_peaks2_unique.bed", 
  genomeVersion='no', 
  outputBedFN = 'negSet_64cell.bed', 
  outputPosFastaFN = 'posSet_64cell.fa',
  outputNegFastaFN = 'negSet_64cell.fa', 
  xfold = 1, 
  repeat_match_tol = 0.02, 
  GC_match_tol = 0.02, 
  length_match_tol = 0.02, 
  batchsize = 5000, 
  nMaxTrials = 20, 
  genome =  genome)


genNullSeqs_ranges(
  "Pm14-MidNeurula_peaks_unique.bed", 
  genomeVersion='no', 
  outputBedFN = 'negSet_midN1.bed', 
  outputPosFastaFN = 'posSet_midN1.fa',
  outputNegFastaFN = 'negSet_midN1.fa', 
  xfold = 1, 
  repeat_match_tol = 0.02, 
  GC_match_tol = 0.02, 
  length_match_tol = 0.02, 
  batchsize = 5000, 
  nMaxTrials = 20, 
  genome =  genome)

genNullSeqs_ranges(
  "Pm01-LateGastrula_peaks_unique.bed", 
  genomeVersion='no', 
  outputBedFN = 'negSet_lateG.bed', 
  outputPosFastaFN = 'posSet_lateG.fa',
  outputNegFastaFN = 'negSet_lateG.fa', 
  xfold = 1, 
  repeat_match_tol = 0.02, 
  GC_match_tol = 0.02, 
  length_match_tol = 0.02, 
  batchsize = 5000, 
  nMaxTrials = 20, 
  genome =  genome)
#  "pot_Late_G_Pm_GGACTCCT_peaks_unique",




setwd("")

anc = import("112cell_30_31_uni_anc_halper2.bed")
original = import("112cell_30_31_uni.bed")
original= original[original@elementMetadata$name %in% anc@elementMetadata$name]
export.bed(original,"112cell_30_31_uni.bed")

anc = import("64cell_28_29_uni_anc_halper2.bed")
original = import("64cell_28_29_uni.bed")
original= original[original@elementMetadata$name %in% anc@elementMetadata$name]
export.bed(original,"64cell_28_29_uni.bed")

anc = import("midN1_25_uni_halper2.bed")
original = import("midN1_25_uni.bed")
original= original[original@elementMetadata$name %in% anc@elementMetadata$name]
export.bed(original,"midN1_25_uni.bed")

anc = import("lateG_26_27_uni_halper2.bed")
original = import("lateG_26_27_uni.bed")
original= original[original@elementMetadata$name %in% anc@elementMetadata$name]
export.bed(original,"lateG_26_27_uni.bed")

anc = import("ef1_common_prenet_uni_anc_halper2.bed")
original = import("ef1_common_prenet_uni.bed")
original= original[original@elementMetadata$name %in% anc@elementMetadata$name]
export.bed(original,"ef1_common_prenet_uni.bed")



setwd("~/Downloads/bed_files/ancestral_peaks/cinte")


anc = import("112cell_30_31_uni_cinte_halper2.bed")
original = import("112cell_30_31_uni_anc_halper2.bed")
original= original[original@elementMetadata$name %in% anc@elementMetadata$name]
export.bed(original,"112cell_30_31_uni_anc_halper2.bed")

anc = import("64cell_28_29_uni_cinte_halper2.bed")
original = import("64cell_28_29_uni_anc_halper2.bed")
original= original[original@elementMetadata$name %in% anc@elementMetadata$name]
export.bed(original,"64cell_28_29_uni_anc_halper2.bed")

anc = import("midN1_25_uni_cinte_halper2.bed")
original = import("midN1_25_uni_halper2.bed")
original= original[original@elementMetadata$name %in% anc@elementMetadata$name]
export.bed(original,"midN1_25_uni_halper2.bed")

anc = import("lateG_26_27_uni_cinte_halper2.bed")
original = import("lateG_26_27_uni_halper2.bed")
original= original[original@elementMetadata$name %in% anc@elementMetadata$name]
export.bed(original,"lateG_26_27_uni_halper2.bed")

anc = import("ef1_common_prenet_uni_cinte_halper2.bed")
original = import("ef1_common_prenet_uni_anc_halper2.bed")
original= original[original@elementMetadata$name %in% anc@elementMetadata$name]
export.bed(original,"ef1_common_prenet_uni_anc_halper2.bed")


it is 10 am now what else can I do?














