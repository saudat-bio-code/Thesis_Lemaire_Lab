#GCA_R <- system.file("extdata", "GCA_R.fa.gz", package="BSgenome")
#GCA_P <- system.file("extdata", "GCA_P.fa.gz", package="BSgenome") junk stuff
#seed_files <- system.file("extdata", "GCA_R.fa.gz", package="BSgenome")
#NAME
#BSgenome.Ciinte.JDB.Roscoff
#BSgenome.Ciinte.JDB.Plymouth
#title
#Cioina intestinalis 2021 GCA_018327825.1_Cint_typeB-Roscoff
#Cioina intestinalis 2021 GCA_018327805.1_Cint_typeB-Plymouth
#Description

library(BSgenome)
library(gkmSVM)
library(Biostrings)
library(peakToGene)
library(BSgenome.Ciinte.JDB.Roscoff)
library(BSgenome.Cintestinalis.KY)

#forgeBSgenomeDataPkg("~/seed") 
#inputBedFN = peaks in bed format, different for each stage?
#try with an old genome?


genNullSeqs(
    "dmel2_3_unique.bed", 
    genomeVersion='dmelT2_3.bed', 
    outputBedFN = 'negSet_dmel2_3_unique.bed.bed', 
    outputPosFastaFN = 'posSet_dmel2_3_unique.bed.fa',
    outputNegFastaFN = 'negSet_dmel2_3_unique.bed.fa', 
    xfold = 1, 
    repeat_match_tol = 0.02, 
    GC_match_tol = 0.02, 
    length_match_tol = 0.02, 
    batchsize = 5000, 
    nMaxTrials = 20, 
    genome =  BSgenome.Ciinte.JDB.Roscoff)

genNullSeqs(
  "dmel4_5_unique.bed", 
  genomeVersion='dmelT2_3.bed', 
  outputBedFN = 'negSet_dmel4_5_unique.bed.bed', 
  outputPosFastaFN = 'posSet_dmel4_5_unique.bed.fa',
  outputNegFastaFN = ' ', 
  xfold = 1, 
  repeat_match_tol = 0.02, 
  GC_match_tol = 0.02, 
  length_match_tol = 0.02, 
  batchsize = 5000, 
  nMaxTrials = 20, 
  genome =  BSgenome.Ciinte.JDB.Roscoff)

genNullSeqs(
  "dmel6_7_unique.bed", 
  genomeVersion='dmelT2_3.bed', 
  outputBedFN = 'negSet_dmel6_7_unique.bed.bed', 
  outputPosFastaFN = 'posSet_dmel6_7_unique.bed.fa',
  outputNegFastaFN = 'negSet_dmel6_7_unique.bed.fa', 
  xfold = 1, 
  repeat_match_tol = 0.02, 
  GC_match_tol = 0.02, 
  length_match_tol = 0.02, 
  batchsize = 5000, 
  nMaxTrials = 20, 
  genome =  BSgenome.Ciinte.JDB.Roscoff)

genNullSeqs(
  "midN1_unique.bed", 
  genomeVersion='dmelT2_3.bed', 
  outputBedFN = 'negSet_midN1_unique.bed.bed', 
  outputPosFastaFN = 'posSet_midN1_unique.bed.fa',
  outputNegFastaFN = 'negSet_midN1_unique.bed.fa', 
  xfold = 1, 
  repeat_match_tol = 0.02, 
  GC_match_tol = 0.02, 
  length_match_tol = 0.02, 
  batchsize = 5000, 
  nMaxTrials = 20, 
  genome =  BSgenome.Ciinte.JDB.Roscoff)


gkmtrain posSet_midN1_unique.bed.fa negSet_midN1_unique.bed.fa -l=10 midN1_model
gkmtrain posSet_dmel6_7_unique.bed.fa negSet_dmel6_7_unique.bed.fa -l=10 dmel6_7_model
gkmtrain posSet_dmel4_5_unique.bed.fa negSet_dmel4_5_unique.bed.fa -l=10 dmel4_5_model
gkmtrain posSet_dmel2_3_unique.bed.fa negSet_dmel2_3_unique.bed.fa -l=10 dmel2_3_model



genNullSeqs(
  "2_june.bed", 
  genomeVersion='dmelT2_3.bed', 
  outputBedFN = 'negSet_2J.bed', 
  outputPosFastaFN = 'posSet2J.fa',
  outputNegFastaFN = 'negSet2J.fa', 
  xfold = 1, 
  repeat_match_tol = 0.02, 
  GC_match_tol = 0.02, 
  length_match_tol = 0.02, 
  batchsize = 5000, 
  nMaxTrials = 20, 
  genome = BSgenome.Cintestinalis.KH.JoinedScaffoldOld)


trf = import("trf_crobu19.bed")

genNullSeqs(
  "peaks_merged_last_stage.bed", 
  genomeVersion='dmelT2_3.bed', 
  outputBedFN = 'negSet_last_stage.bed', 
  outputPosFastaFN = 'posSet_last_stage.fa',
  outputNegFastaFN = 'negSet_last_stage.fa', 
  xfold = 1, 
  repeat_match_tol = 0.02, 
  GC_match_tol = 0.02, 
  length_match_tol = 0.02, 
  batchsize = 5000, 
  nMaxTrials = 20, 
  genome = BSgenome.Cintestinalis.KY)

library('peakToGene')
peaks_merged_last_stage<-import("~/Desktop/peaks_merged_last_stage.bed")

dmel2_3_unique.bed = import("~/Downloads/bed_files/new/dmel2_3_unique.bed")
dmel6_7_unique.bed = import("~/Downloads/bed_files/new/dmel6_7_unique.bed")
dmel6_7_unique.bed = import("merged_edited_clean.bed")

features2020 = getFeatures("~/Downloads/liftoff2019ref_2021regions.gff3")
features2019 = getFeatures("~/Downloads/HT.Gene.gff3")
features2020_gr = import("~/Downloads/liftoff2019ref_2021regions.gff3")
features2020_gr = import("~/Downloads/liftoff2019ref_2021regions.gff3")

annoData <- toGRanges(features2020, feature="gene")
features2019_gr = import("~/Downloads/HT.Gene.gff3")



overlaps.anno <- annotatePeakInBatch(dmel2_3_unique.bed, AnnotationData= features2020_gr, 
                                     output="overlapping", maxgap=5000L)
overlaps.anno <- annotatePeakInBatch(dmel6_7_unique.bed, AnnotationData= features2020_gr, 
                                     output="overlapping", maxgap=5000L,  multiple = F)

overlaps_dmel2_3_unique.bed <- lapply(features2020,findOverlaps,dmel2_3_unique.bed)

overlaps_peaks_merged_last_stage <- lapply(features2,getOverlaps,dmel6_7_unique.bed)

# YAHOO I made it through~~~


# Tomorrow I will be generating peak files
#or I can go home and do it using Galaxy
#I need to align all files so reproducibility
#yahoo I got more time!
gr1.copy <- gr1
gr1.copy$score <- 1
binOverFeature(gr1, gr1.copy, annotationData=annoData,
               radius=5000, nbins=10, FUN=c(sum, length),
               ylab=c("score", "count"), 
               main=c("Distribution of aggregated peak scores around TSS", 
                      "Distribution of aggregated peak numbers around TSS"))

