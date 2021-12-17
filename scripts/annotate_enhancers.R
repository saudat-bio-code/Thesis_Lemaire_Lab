setwd("~/Downloads/bed_files/merging/")
cell112_30_31_uni.bed =import("112cell_30_31_uni.bed")
cell64_28_29_uni.bed<-import("64cell_28_29_uni.bed")
ef1_common_prenet_uni.bed<-import("ef1_common_prenet_uni.bed") #these two got problems
midN1_25_uni.bed<-import("midN1_25_uni.bed")#these two got problems
lateG_26_27_uni.bed<-import("lateG_26_27_uni.bed")
#####enhancers only
annotated_peak<- annotatePeakInBatch(ef1_common_prenet_uni.bed, AnnotationData= features2020_gr, 
                                     output="overlapping", maxgap=5000L,  multiple = F)
upstream_annotated_peak = annotated_peak[(elementMetadata(annotated_peak)[,"insideFeature"] %in% "upstream" ),]
upstream_annotated_peak = upstream_annotated_peak[!duplicated(elementMetadata(upstream_annotated_peak)["name"]),]
upstream_annotated_peak #ef1_common_prenet_uni.bed
export.bed(upstream_annotated_peak,con='enhancers/ef1_common_prenet_uni.bed')


annotated_peak<- annotatePeakInBatch(cell64_28_29_uni.bed, AnnotationData= features2020_gr, 
                                     output="overlapping", maxgap=5000L,  multiple = F)
upstream_annotated_peak = annotated_peak[(elementMetadata(annotated_peak)[,"insideFeature"] %in% "upstream" ),]
upstream_annotated_peak = upstream_annotated_peak[!duplicated(elementMetadata(upstream_annotated_peak)["name"]),]
upstream_annotated_peak #cell64_28_29_uni.bed
export.bed(upstream_annotated_peak,con='enhancers/64cell_28_29_uni.bed')

annotated_peak<- annotatePeakInBatch(cell112_30_31_uni.bed, AnnotationData= features2020_gr, 
                                     output="overlapping", maxgap=5000L,  multiple = F)
upstream_annotated_peak = annotated_peak[(elementMetadata(annotated_peak)[,"insideFeature"] %in% "upstream" ),]
upstream_annotated_peak = upstream_annotated_peak[!duplicated(elementMetadata(upstream_annotated_peak)["name"]),]
upstream_annotated_peak #cell112_30_31_uni.bed
export.bed(upstream_annotated_peak,con='enhancers/112cell_30_31_uni.bed')

annotated_peak<- annotatePeakInBatch(midN1_25_uni.bed, AnnotationData= features2020_gr, 
                                     output="overlapping", maxgap=5000L,  multiple = F)
upstream_annotated_peak = annotated_peak[(elementMetadata(annotated_peak)[,"insideFeature"] %in% "upstream" ),]
upstream_annotated_peak = upstream_annotated_peak[!duplicated(elementMetadata(upstream_annotated_peak)["name"]),]
upstream_annotated_peak #midN1_25_uni.bed
export.bed(upstream_annotated_peak,con='enhancers/midN1_25_uni.bed')

annotated_peak<- annotatePeakInBatch(lateG_26_27_uni.bed, AnnotationData= features2020_gr, 
                                     output="overlapping", maxgap=5000L,  multiple = F)
upstream_annotated_peak = annotated_peak[(elementMetadata(annotated_peak)[,"insideFeature"] %in% "upstream" ),]
upstream_annotated_peak = upstream_annotated_peak[!duplicated(elementMetadata(upstream_annotated_peak)["name"]),]
upstream_annotated_peak #lateG_26_27_uni.bed
export.bed(upstream_annotated_peak,con='enhancers/lateG_26_27_uni.bed')


setwd("~/Downloads/bed_files/merging/converted/")
cell112_30_31_uni.bed =import("112cell_30_31_uni.bed")
cell64_28_29_uni.bed<-import("64cell_28_29_uni.bed")
ef1_common_prenet_uni.bed<-import("ef1_common_prenet_uni.bed") #these two got problems
midN1_25_uni.bed<-import("midN1_25_uni.bed")#these two got problems
lateG_26_27_uni.bed<-import("lateG_26_27_uni.bed")

annotated_peak<- annotatePeakInBatch(ef1_common_prenet_uni.bed, AnnotationData= features2019_gr, 
                                     output="overlapping", maxgap=5000L,  multiple = F)
upstream_annotated_peak = annotated_peak[(elementMetadata(annotated_peak)[,"insideFeature"] %in% "upstream" ),]
upstream_annotated_peak = upstream_annotated_peak[!duplicated(elementMetadata(upstream_annotated_peak)["name"]),]
upstream_annotated_peak #ef1_common_prenet_uni.bed
export.bed(upstream_annotated_peak,con='enhancers/ef1_common_prenet_uni.bed')

annotated_peak<- annotatePeakInBatch(cell64_28_29_uni.bed, AnnotationData= features2019_gr, 
                                     output="overlapping", maxgap=5000L,  multiple = F)
upstream_annotated_peak = annotated_peak[(elementMetadata(annotated_peak)[,"insideFeature"] %in% "upstream" ),]
upstream_annotated_peak = upstream_annotated_peak[!duplicated(elementMetadata(upstream_annotated_peak)["name"]),]
upstream_annotated_peak #cell64_28_29_uni.bed
export.bed(upstream_annotated_peak,con='enhancers/64cell_28_29_uni.bed')

annotated_peak<- annotatePeakInBatch(cell112_30_31_uni.bed, AnnotationData= features2019_gr, 
                                     output="overlapping", maxgap=5000L,  multiple = F)
upstream_annotated_peak = annotated_peak[(elementMetadata(annotated_peak)[,"insideFeature"] %in% "upstream" ),]
upstream_annotated_peak = upstream_annotated_peak[!duplicated(elementMetadata(upstream_annotated_peak)["name"]),]
upstream_annotated_peak #cell112_30_31_uni.bed
export.bed(upstream_annotated_peak,con='enhancers/112cell_30_31_uni.bed')

annotated_peak<- annotatePeakInBatch(midN1_25_uni.bed, AnnotationData= features2019_gr, 
                                     output="overlapping", maxgap=5000L,  multiple = F)
upstream_annotated_peak = annotated_peak[(elementMetadata(annotated_peak)[,"insideFeature"] %in% "upstream" ),]
upstream_annotated_peak = upstream_annotated_peak[!duplicated(elementMetadata(upstream_annotated_peak)["name"]),]
upstream_annotated_peak #midN1_25_uni.bed
export.bed(upstream_annotated_peak,con='enhancers/midN1_25_uni.bed')

annotated_peak<- annotatePeakInBatch(lateG_26_27_uni.bed, AnnotationData= features2019_gr, 
                                     output="overlapping", maxgap=5000L,  multiple = F)
upstream_annotated_peak = annotated_peak[(elementMetadata(annotated_peak)[,"insideFeature"] %in% "upstream" ),]
upstream_annotated_peak = upstream_annotated_peak[!duplicated(elementMetadata(upstream_annotated_peak)["name"]),]
upstream_annotated_peak #lateG_26_27_uni.bed
export.bed(upstream_annotated_peak,con='enhancers/lateG_26_27_uni.bed')
##############################################################################################################################



