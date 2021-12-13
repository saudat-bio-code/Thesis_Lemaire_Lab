##INPUTS 
data<-read.table('/srv/scratch/annashch/stemchell/het/output/remapped_merged_peaks/all.idr.bed.merged',header=F) 
pr_files=c("~/ef1_alpha_16_june.bw",
"~/ef1_alpha_2_june.bw")
names(pr_files)=c("ef1_16", "ef1_2")
outputfname="crobu_counts.regression.csv"
readlength=51 

#############################################################################

library('GenomicRanges') 
library('rtracklayer') 

colnames(data)<-c("chr","start","end")
bed<-with(data,GRanges(chr,IRanges(start+1,end)))
bw=BigWigFileList(pr_files) 
counts <- matrix(NA, nrow = length(bed), ncol = length(bw))
colnames(counts) <- names(bw)
chromnames=levels(seqnames(bed))

for(i in seq_len(length(bw))) 
{
    last=1 	 
    coverage <- import(bw[[i]], as = 'RleList')
    for(j in chromnames)
    {
    range_vals=ranges(bed[seqnames(bed)==j])
    cur_coverage=coverage[[j]] 
    if(is.null(cur_coverage))
    {
    counts[last:(last+length(range_vals)-1),i]=matrix(0,nrow=length(range_vals),ncol=1)
    }
    else
    {
    newvals=sum(Views(cur_coverage, ranges(bed[seqnames(bed)==j])))
    counts[last:(last+length(newvals)-1), i] <-newvals 
    }   
    last=last+length(range_vals) 
    }

}
## Divide by read length and round to integer numbers
counts <- round(counts / readlength, 0)
write.csv(counts,file=outputfname)

