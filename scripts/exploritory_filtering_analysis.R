

#so I just realised there might have been issues with a lot of repeats but this would have affected not the right files
#One cool way to find out if I am doing the right thing with my analysis is to repeat the steps but to align it on the same genome
#chieck with lifoff
IRanges::findOverlaps( ef1_common_prenet_uni.bed,features2020)
txdb<-GenomicFeatures::makeTxDbFromGFF("~/la",organism="Gallus gallus",taxonomyId=9031)


setwd("~/")
trf = import("2021-Roscoff.fna.2.7.7.80.10.50.500.bed") 
setwd("~/Downloads/bed_files/test_repeats_neg_set/")
cell112_30_31_uni.bed =import("112cell_30_31_uni.bed")
cell64_28_29_uni.bed<-import("64cell_28_29_uni.bed")
ef1_common_prenet_uni.bed<-import("ef1_common_prenet_uni.bed") #these two got problems
midN1_25_uni.bed<-import("midN1_25_uni.bed")#these two got problems
lateG_26_27_uni.bed<-import("lateG_26_27_uni.bed")
IRanges::findOverlaps( cell112_30_31_uni.bed,trf) #1254/5514
IRanges::findOverlaps( cell64_28_29_uni.bed,trf) #1026/4509
IRanges::findOverlaps( midN1_25_uni.bed,trf)# 1612/14225
IRanges::findOverlaps( lateG_26_27_uni.bed,trf) #1092/9228
IRanges::findOverlaps( ef1_common_prenet_uni.bed,trf)#297/10473




#### duplications
setwd("~/Downloads/bed_files/")
lateG_26_27_uni.bed <- read.table("lateG_26_27_uni.bed",sep = "\t" )
midN1_25_uni.bed <- read.table("midN1_25_uni.bed",sep = "\t" )
cell112_30_31_uni.bed <- read.table("112cell_30_31_uni.bed",sep = "\t" )
cell64_28_29_uni.bed <- read.table("64cell_28_29_uni.bed",sep = "\t" )
ef1_common_prenet_uni.bed <- read.table("ef1_common_prenet_uni.bed",sep = "\t" )

colnames(ef1_common_prenet_uni.bed)<-c("chr", "start", "end", "name","score","strand","s1","s2","s3","s4")
colnames(cell64_28_29_uni.bed)<-c("chr", "start", "end", "name","score","strand","s1","s2","s3","s4")
colnames(cell112_30_31_uni.bed)<-c("chr", "start", "end", "name","score","strand","s1","s2","s3","s4")
colnames(midN1_25_uni.bed)<-c("chr", "start", "end", "name","score","strand","s1","s2","s3","s4")
colnames(lateG_26_27_uni.bed)<-c("chr", "start", "end", "name","score","strand","s1","s2","s3","s4")

ef1_common_prenet_uni.bed=ef1_common_prenet_uni.bed[!duplicated(ef1_common_prenet_uni.bed$name),]
table(duplicated(ef1_common_prenet_uni.bed$start))
table(duplicated(midN1_25_uni.bed$start))
table(duplicated(cell64_28_29_uni.bed$start))
table(duplicated(cell112_30_31_uni.bed$start))
table(duplicated(lateG_26_27_uni.bed$start))



trf = import("2021-Roscoff.fna.2.7.7.80.10.50.500.bed") 

setwd("~")
trf = import("HT.Ref.fasta.2.7.7.80.10.50.500.bed")
setwd("~/Downloads/bed_files/test_repeats_neg_set/")
cell112_30_31_uni.bed =import("112cell_30_31_uni_cut.bed")
cell64_28_29_uni.bed<-import("64cell_28_29_uni_cut.bed")
ef1_common_prenet_uni.bed<-import("ef1_common_prenet_uni_cut.bed") #these two got problems
midN1_25_uni.bed<-import("midN1_25_uni_cut.bed")#these two got problems
lateG_26_27_uni.bed<-import("lateG_26_27_uni_cut.bed")
IRanges::findOverlaps( cell112_30_31_uni.bed,trf) # 382/3681
IRanges::findOverlaps( cell64_28_29_uni.bed,trf) #2143/3146
IRanges::findOverlaps( midN1_25_uni.bed,trf) #2142/8540
IRanges::findOverlaps( lateG_26_27_uni.bed,trf) # 581/6351
IRanges::findOverlaps( ef1_common_prenet_uni.bed,trf) #151/10000
#numbers too good? is this repeated dataset or not?
enNullSeqs_ranges(
  "112cell_30_31_uni_cut.bed", 
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

negset_cell112_30_31_uni.bed =import("negSet_112cell.bed")
IRanges::findOverlaps( negset_cell112_30_31_uni.bed,trf) # 251 of 382??
IRanges::findOverlaps( negset_cell112_30_31_uni.bed@ranges,trf) # 247 of 382??
IRanges::findOverlaps( negset_cell112_30_31_uni.bed@ranges,trf) #  I get 404 at randoom lol! how does this work, now I get 444
#it is #394 in enNulSeqs()
#will this affect classification in positive or negatve way?


setwd("~")
trf = import("HT.Ref.fasta.2.7.7.80.10.50.500.bed")
setwd("~/Downloads/bed_files/test_repeats_neg_set/")

negSet_112cell.bed =import("negSet_112cell.bed")
negSet_64cell.bed =import("negSet_64cell.bed")
negSet_ef1.bed =import("negSet_ef1.bed")
negSet_lateG.bed =import("negSet_lateG.bed")
negSet_midN1.bed =import("negSet_midN1.bed")
#random/trf/actual/total sqs
IRanges::findOverlaps( negSet_112cell.bed,trf) # 419 / 236 / 382 / 2602
IRanges::findOverlaps( negSet_64cell.bed,trf) # 2129 / 1477 / 2143 / 13725
IRanges::findOverlaps( negSet_ef1.bed,trf) # 1035 / 59 / 151 / 6321 
IRanges::findOverlaps( negSet_lateG.bed,trf) # 526 / 388 / 581 / 3722
IRanges::findOverlaps( negSet_midN1.bed,trf) #2789 /1068 / 2172 / 12660

library(caret)
heatmap = read.csv("heatmap_performance - Sheet1(3).csv")
heatmap_m = mutate(heatmap)
#heatmap_mc = as.matrix(heatmap_m)
heatmap_mx = confusionMatrix(heatmap_m) # conf_mat_tab
heatmap(heatmap_mx)

heatmap_cl = heatmap[,-1]
heatmap_cl_m = as.matrix(heatmap_cl)

library("ggplot2")
ggplot(hm, aes(x=heatmap_cl, y=new, fill=value)) + geom_tile() 

+ scale_x_discrete(limits = rev(levels(heatmap_cl)))

hm <- hm %>%
  mutate(x = factor(x), # alphabetical order by default
         heatmap_cl = factor(heatmap_cl, levels = rev(unique(heatmap_cl)))) # force reverse alphabetical order

################some more filgering 
setwd("~/Downloads/bed_files/ancestral_peaks//")

cell112_30_31_uni.bed =import("112cell_30_31_uni.bed")
cell64_28_29_uni.bed<-import("64cell_28_29_uni.bed")
ef1_common_prenet_uni.bed<-import("ef1_common_prenet_uni.bed") #these two got problems
midN1_25_uni.bed<-import("midN1_25_uni.bed")#these two got problems
lateG_26_27_uni.bed<-import("lateG_26_27_uni.bed")

cell112_30_31_uni_anc.bed =import("112cell_30_31_uni_anc.bed")
cell64_28_29_uni_anc.bed<-import("64cell_28_29_uni_anc.bed")
ef1_common_prenet_uni_anc.bed<-import("ef1_common_prenet_uni_anc.bed") #these two got problems
midN1_25_uni_anc.bed<-import("midN1_25_uni_anc.bed")#these two got problems
lateG_26_27_uni_anc.bed<-import("lateG_26_27_uni_anc.bed")

table(cell112_30_31_uni.bed@elementMetadata$name %in% cell112_30_31_uni_anc.bed@elementMetadata$name)

table(cell112_30_31_uni_anc.bed@elementMetadata$name %in% cell112_30_31_uni.bed@elementMetadata$name) 


table( cell64_28_29_uni_anc.bed@elementMetadata$name  %in% cell64_28_29_uni.bed@elementMetadata$name   )
table( ef1_common_prenet_uni_anc.bed@elementMetadata$name  %in%  ef1_common_prenet_uni.bed@elementMetadata$name  )
table( midN1_25_uni_anc.bed@elementMetadata$name  %in%  midN1_25_uni.bed@elementMetadata$name  )
table( lateG_26_27_uni_anc.bed@elementMetadata$name  %in%  lateG_26_27_uni.bed@elementMetadata$name  )

table(duplicated(cell112_30_31_uni_anc.bed@elementMetadata$name))
table(duplicated(  cell64_28_29_uni_anc.bed@elementMetadata$name))
table(duplicated(  ef1_common_prenet_uni_anc.bed@elementMetadata$name))
table(duplicated(  midN1_25_uni_anc.bed@elementMetadata$name))
table(duplicated(  lateG_26_27_uni_anc.bed@elementMetadata$name))



################## RNA - ATAC correlations

setwd("~/Downloads/bed_files/rna_atac_correlations")
sorted_all_conditions.bed = fread("sorted_all_conditions.bed")
colnames(sorted_all_conditions.bed)<-c("chr", "start", "end", "name","score","strand","source","type", "note", "ID")
sorted_all_conditions.bed = toGRanges(sorted_all_conditions.bed)



setwd("~/Downloads")
all_conditions_FPKM_results.csv = read.csv("all_conditions_FPKM_results.csv")
KH2012 = import("KHNCBI.2018.Gene.gff3")
ranges_all_conditions = KH2012[KH2012@elementMetadata$ID %in% all_conditions_FPKM_results.csv$X]
#this file is perfectly what I needed, the only problem I have with it is that I need to have RNA results saved in $score column
head(ranges_all_conditions) #maybe I can fix it if I need that info
#now this is ready and I can save it
export.bed(ranges_all_conditions , "all_condtions.bed")
export.gff3(ranges_all_conditions , "all_condtions.bed")

#importing convereted files yay
cell112_30_31_uni.bed = fread("112cell_30_31_uni_genes.bed")
cell64_28_29_uni.bed = fread("64cell_28_29_uni_genes.bed")
midN1_25_uni.bed = fread("midN1_25_uni_genes.bed")
lateG_26_27_uni.bed =fread("lateG_26_27_uni_genes.bed") 
ef1_common_prenet_uni.bed =fread("ef1_common_prenet_uni_genes.bed")

all_conditions_FPKM_results.csv$label = "other"

all_conditions_FPKM_results.csv$label = ifelse(all_conditions_FPKM_results.csv$X %in% common_deve_genes, "Filtered", all_conditions_FPKM_results.csv$label)

#cell112
peak_clean = cell112_30_31_uni.bed[!duplicated(cell112_30_31_uni.bed$V4),V4]
peak_clean = peak_clean[!duplicated(peak_clean)]
all_conditions_FPKM_results.csv$label = ifelse(all_conditions_FPKM_results.csv$X %in% peak_clean, "B", all_conditions_FPKM_results.csv$label)
#lateg
peak_clean = lateG_26_27_uni.bed[!duplicated(lateG_26_27_uni.bed$V4),V4]
peak_clean = peak_clean[!duplicated(peak_clean)]
all_conditions_FPKM_results.csv$label = ifelse(all_conditions_FPKM_results.csv$X %in% peak_clean , "C", all_conditions_FPKM_results.csv$label)
#midn1
peak_clean = midN1_25_uni.bed[!duplicated(midN1_25_uni.bed$V4),V4]
peak_clean = peak_clean[!duplicated(peak_clean)]
all_conditions_FPKM_results.csv$label = ifelse(all_conditions_FPKM_results.csv$X %in% peak_clean, "D", all_conditions_FPKM_results.csv$label)
#cell64
peak_clean = cell64_28_29_uni.bed[!duplicated(cell64_28_29_uni.bed$V4),V4]
peak_clean = peak_clean[!duplicated(peak_clean)]
all_conditions_FPKM_results.csv$label = ifelse(all_conditions_FPKM_results.csv$X %in% peak_clean , "A", all_conditions_FPKM_results.csv$label)
#ef1
peak_clean = ef1_common_prenet_uni.bed[!duplicated(ef1_common_prenet_uni.bed$V4),V4]
peak_clean = peak_clean[!duplicated(peak_clean)]
all_conditions_FPKM_results.csv$label = ifelse(all_conditions_FPKM_results.csv$X %in% peak_clean, "E", all_conditions_FPKM_results.csv$label)
all_conditions_FPKM_results.csv$label = ifelse(all_conditions_FPKM_results.csv$X %in% common_deve_genes, "Filtered", all_conditions_FPKM_results.csv$label)


table(all_conditions_FPKM_results.csv$label)


#cell64
cell64_28_29_uni.bed = cell64_28_29_uni.bed[!duplicated(cell64_28_29_uni.bed$V4),V4]
cell64_28_29_uni.bed = cell64_28_29_uni.bed[!duplicated(cell64_28_29_uni.bed)]
#cell112
cell112_30_31_uni.bed = cell112_30_31_uni.bed[!duplicated(cell112_30_31_uni.bed$V4),V4]
cell112_30_31_uni.bed = cell112_30_31_uni.bed[!duplicated(cell112_30_31_uni.bed)]
#midn1
midN1_25_uni.bed = midN1_25_uni.bed[!duplicated(midN1_25_uni.bed$V4),V4]
midN1_25_uni.bed = midN1_25_uni.bed[!duplicated(midN1_25_uni.bed)]
#ef1
ef1_common_prenet_uni.bed = ef1_common_prenet_uni.bed[!duplicated(ef1_common_prenet_uni.bed$V4),V4]
ef1_common_prenet_uni.bed = ef1_common_prenet_uni.bed[!duplicated(ef1_common_prenet_uni.bed)]
#lateg
lateG_26_27_uni.bed = lateG_26_27_uni.bed[!duplicated(lateG_26_27_uni.bed$V4),V4]
lateG_26_27_uni.bed = lateG_26_27_uni.bed[!duplicated(lateG_26_27_uni.bed)]

gene_list = list(cell64 =cell64_28_29_uni.bed, cell112 = cell112_30_31_uni.bed, midn1 = midN1_25_uni.bed, lateG = lateG_26_27_uni.bed, ef1 = ef1_common_prenet_uni.bed)
library(ggvenn)
ggvenn(
  gene_list, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

#gene_list = list(cell64 =cell64_28_29_uni.bed, cell112 = cell112_30_31_uni.bed, midn1 = midN1_25_uni.bed, lateG = lateG_26_27_uni.bed, ef1 = ef1_common_prenet_uni.bed)
#gene_list2 = list(cell64 =cell64_28_29_uni.bed, cell112 = cell112_30_31_uni.bed, midn1 = midN1_25_uni.bed)

#common_deve_genes = Reduce(intersect, gene_list)
#common_deve_genes2 = Reduce(intersect, gene_list2)

gene_list = c(cell64_28_29_uni.bed,cell112_30_31_uni.bed,midN1_25_uni.bed,ef1_common_prenet_uni.bed,lateG_26_27_uni.bed, ef1_common_prenet_uni.bed)
common_deve_genes = gene_list[duplicated(gene_list)]

#first we sort
all_conditions_FPKM_results.csv = all_conditions_FPKM_results.csv[order(all_conditions_FPKM_results.csv$label),]
head(all_conditions_FPKM_results.csv)
#rownames
rownames(all_conditions_FPKM_results.csv) = all_conditions_FPKM_results.csv$X
#then we get rid of extra columns
all_conditions_FPKM_filtered = all_conditions_FPKM_results.csv[,-c(1,2,3,6,7,14,15,16,17)]
head(all_conditions_FPKM_filtered)
pheatmap(all_conditions_FPKM_filtered)

library(cluster)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(csaw)
library(pheatmap)
library(limma)
library(tidyr)

count.matrix = as.matrix(all_conditions_FPKM_filtered)
#dds = DESeqDataSetFromMatrix(all_conditions_FPKM_filtered_m, colData = colnames(all_conditions_FPKM_filtered)) 
#vsd = vst(all_conditions_FPKM_filtered, blind = FALSE) # first transform the data
#df <- as.data.frame(colData(ddsM1)[,c("Condition","sex", "age", "SampleID", "Batch")], stringsAsFactors = FALSE)
#rownames(df) = paste0(df$SampleID, "_", df$age)
#df = df[,c("Batch", "age", "Condition")]
#count.matrix <- removeBatchEffect(assay(vsd), batch = colData(ddsM1)$Batch, design = model.matrix( ~ colData(ddsM1)$Condition))
# generate Z-scores from counts

count.matrix = as.matrix(all_conditions_FPKM_filtered)
Means.vsd = rowMeans(count.matrix)
SD.vsd = rowSds(count.matrix)
norm.vsd = (count.matrix - Means.vsd) / SD.vsd
#colnames(norm.vsd) = rownames(all_conditions_FPKM_filtered)
x = norm.vsd
clusters = pam(x, 2) # use pam clustering here to obtain consistent clusters, you canHEA choose the cluster number 
#x = x[names(sort(clusters$clustering)),]
pheatmap(x, cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames=T,
         show_colnames = T)

