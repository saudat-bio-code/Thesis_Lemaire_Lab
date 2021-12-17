xsetwd("~/svm_all/enhancers_crobu/")
setwd("~/svm_all/all_crobu/")
setwd("~/svm_all/all_ciinte/")
setwd("~/svm_all/all_ciinte/moodel")
setwd("~/svm_all/all_ciinte/long_without_trf/")
setwd("~/svm_all/all_ciinte/no_trf")
setwd("~/svm_all/all_ciinte/no_trf_mixed")
setwd("~/svm_all/repeated_macs_crobu_ancestry/")
setwd("~/svm_all/repeated_macs_crobu")
setwd("~/Downloads/bed_files/test_repeats_neg_set/without_trf")
setwd("~/Downloads/bed_files/changed_t_param_plus_ancestr")

library("data.table")
library("RColorBrewer") 
library("plyr")
library("qvalue")
library("topGO")
library("ROCR")
library("tseries")

pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3")) 
rocPrArea<-function(positive, negative) {

  positive$state<-rep(1,tiems=nrow(positive))
  negative$state<-rep(0,tiems=nrow(negative))
  aucResult<-c()
  pr_plotValueX<-list()
  pr_plotValueY<-list()

  for (i in 2:5) {
    pos<-cbind(positive[[i]],positive[,6])
    neg<-cbind(negative[[i]],negative[,6])
    tempData<-rbind(pos,neg)
    pred <- prediction(tempData$V1,tempData$state) 
    perf <- performance( pred, "tpr", "fpr" )
    
    auc_temp <-performance( pred, measure = "auc")
    aucResult[i-1]<-round(unlist(slot(auc_temp, "y.values")),digits = 3)

    pr_temp<-data.frame(unlist(perf@x.values),unlist(perf@y.values))
    names(pr_temp)<-c("x","y")
    pr_temp<-na.omit(pr_temp)
    pr_plotValueX[[i-1]]<-pr_temp$x
    pr_plotValueY[[i-1]]<-pr_temp$y
  } 
  return(list(pr_plotValueX,pr_plotValueY,aucResult))
}


dev.off()
pdf("plots_cmam_all.pdf")
#pdf("plots_ciinte_long.pdf")
##### take T1 enhancers as an example #####
## 1_stage 64_cell cross validation 
cv<-fread("64cell_model.cvpred.txt")
crossV<-cv[,c(2,3)]
## T1 enhancers prediction based on 2_stage 112_cell 
pos2<-fread("64cellfasta_bypr112cell.txt")
neg2<-fread("NS_64cellfasta_bypr112cell.txt")
## T1 enhancers prediction based on 3_stage late_gastrula 
pos3<-fread("64cellfasta_byprlateG.txt") 
neg3<-fread("NS_64cellfasta_byprlateG.txt")
## T1 enhancers prediction based on 4_stage mid_neurula 
pos4<-fread("64cellfasta_byprmidN.txt")
neg4<-fread("NS_64cellfasta_byprmidN.txt")
## T1 enhancers prediction based on 5_stage mid_tailbud 
pos5<-fread("64cellfasta_bypref1.txt")
neg5<-fread("NS_64cellfasta_bypref1.txt")

## calculate AUC
# 1_stage 64_cell cross validation 
pred <- prediction(crossV$V2, crossV$V3) 
perf <- performance( pred, "tpr", "fpr" )
auc_temp <-performance( pred, measure = "auc")
aucT1<-unlist(slot(auc_temp, "y.values"))
# T1 enhancers prediction based on other models
pos<-cbind(pos2,pos3$V2,pos4$V2,pos5$V2)
neg<-cbind(neg2,neg3$V2,neg4$V2,neg5$V2)
stage1<-rocPrArea(pos,neg)
aucOthers<-list(stage1[[1]],stage1[[2]],stage1[[3]])

## plot
plotColors<-c(pal[3],pal[7],pal[4],pal[5],pal[2])
plotLegend<-c("1_stage 64_cell", "2_stage 112_cell", "3_stage late_gastrula", "4_stage mid_neurula", "5_stage mid_tailbud")

par(mar=c(5,5,2,2))
plot(unlist(perf@x.values),unlist(perf@y.values),lwd = 3,cex=1.5,col=plotColors[1],main="Ciona enhancers",xlab="False positive rate",ylab="Ture positive rate",
     cex.axis=1.5,cex.lab=1.5,cex.main=1.5,type="l")
legend(0.2,0.5,paste(plotLegend[1]," ","(AUC = ",round(aucT1,digits = 3),")",sep=""),lwd = 2,bty="n",cex=1.4,col = plotColors[1]) 

for (i in 1:4) {
  lines(unlist(aucOthers[[1]][i]),unlist(aucOthers[[2]][i]),col=plotColors[i+1],lwd=3)
}

legendY<-c(0.4,0.3,0.2,0.1)
for (i in 1:4) {
  legendName<-paste(plotLegend[i+1]," ","(AUC = ",round(aucOthers[[3]][i],digits = 3),")",sep="")
  legend(0.2,legendY[i],legendName,lwd = 3,bty="n",cex=1.4,col = plotColors[[i+1]]) 
  
}

####stage 2

##### take T1 enhancers as an example #####
## 1_stage 64_cell cross validation 
cv<-fread("112cell_model.cvpred.txt")
crossV<-cv[,c(2,3)]
## T1 enhancers prediction based on 2_stage 112_cell 
pos2<-fread("112cellfasta_bypr64cell.txt")
neg2<-fread("NS_112cellfasta_bypr64cell.txt")
## T1 enhancers prediction based on 3_stage late_gastrula 
pos3<-fread("112cellfasta_byprlateG.txt")
neg3<-fread("NS_112cellfasta_byprlateG.txt")
## T1 enhancers prediction based on 4_stage mid_neurula 
pos4<-fread("112cellfasta_byprmidN.txt")
neg4<-fread("NS_112cellfasta_byprmidN.txt")
## T1 enhancers prediction based on 5_stage mid_tailbud 
pos5<-fread("112cellfasta_bypref1.txt")
neg5<-fread("NS_112cellfasta_bypref1.txt")

## calculate AUC
# 1_stage 64_cell cross validation 
pred <- prediction(crossV$V2, crossV$V3) 
perf <- performance( pred, "tpr", "fpr" )
auc_temp <-performance( pred, measure = "auc")
aucT1<-unlist(slot(auc_temp, "y.values"))
# T1 enhancers prediction based on other models
pos<-cbind(pos2,pos3$V2,pos4$V2,pos5$V2)
neg<-cbind(neg2,neg3$V2,neg4$V2,neg5$V2)
stage1<-rocPrArea(pos,neg)
aucOthers<-list(stage1[[1]],stage1[[2]],stage1[[3]])

## plot
plotColors<-c(pal[3],pal[7],pal[4],pal[5],pal[2])
plotLegend<-c("2_stage 112_cell","1_stage 64_cell", "3_stage late_gastrula", "4_stage mid_neurula", "5_stage mid_tailbud")

par(mar=c(5,5,2,2))
plot(unlist(perf@x.values),unlist(perf@y.values),lwd = 3,cex=1.5,col=plotColors[1],main="Ciona enhancers",xlab="False positive rate",ylab="Ture positive rate",
     cex.axis=1.5,cex.lab=1.5,cex.main=1.5,type="l")
legend(0.2,0.5,paste(plotLegend[1]," ","(AUC = ",round(aucT1,digits = 3),")",sep=""),lwd = 3,bty="n",cex=1.4,col = plotColors[1]) 

for (i in 1:4) {
  lines(unlist(aucOthers[[1]][i]),unlist(aucOthers[[2]][i]),col=plotColors[i+1],lwd=3)
}

legendY<-c(0.4,0.3,0.2,0.1)
for (i in 1:4) {
  legendName<-paste(plotLegend[i+1]," ","(AUC = ",round(aucOthers[[3]][i],digits = 3),")",sep="")
  legend(0.2,legendY[i],legendName,lwd = 3,bty="n",cex=1.4,col = plotColors[[i+1]]) 
  
}

####stage 3


##### take T1 enhancers as an example #####
## 1_stage 64_cell cross validation 
cv<-fread("lateG_model.cvpred.txt")
crossV<-cv[,c(2,3)]
## T1 enhancers prediction based on 2_stage 112_cell 
pos2<-fread("lateGfasta_bypr64cell.txt")
neg2<-fread("NS_lateGfasta_bypr64cell.txt")
## T1 enhancers prediction based on 3_stage late_gastrula 
pos3<-fread("lateGfasta_byprmidN.txt")
neg3<-fread("NS_lateGfasta_byprmidN.txt")
## T1 enhancers prediction based on 4_stage mid_neurula 
pos4<-fread("lateGfasta_bypr112cell.txt")
neg4<-fread("NS_lateGfasta_bypr112cell.txt")
## T1 enhancers prediction based on 5_stage mid_tailbud 
pos5<-fread("lateGfasta_bypref1.txt")
neg5<-fread("NS_lateGfasta_bypref1.txt")

## calculate AUC
# 1_stage 64_cell cross validation 
pred <- prediction(crossV$V2, crossV$V3) 
perf <- performance( pred, "tpr", "fpr" )
auc_temp <-performance( pred, measure = "auc")
aucT1<-unlist(slot(auc_temp, "y.values"))
# T1 enhancers prediction based on other models
pos<-cbind(pos2,pos3$V2,pos4$V2,pos5$V2)
neg<-cbind(neg2,neg3$V2,neg4$V2,neg5$V2)
stage1<-rocPrArea(pos,neg)
aucOthers<-list(stage1[[1]],stage1[[2]],stage1[[3]])

## plot
plotColors<-c(pal[3],pal[7],pal[4],pal[5],pal[2])
plotLegend<-c("3_stage late_gastrula", "1_stage 64_cell", "2_stage 112_cell", "4_stage mid_neurula", "5_stage mid_tailbud")

par(mar=c(5,5,2,2))
plot(unlist(perf@x.values),unlist(perf@y.values),lwd = 3,cex=1.5,col=plotColors[1],main="Ciona enhancers",xlab="False positive rate",ylab="Ture positive rate",
cex.axis=1.5,cex.lab=1.5,cex.main=1.5,type="l")
legend(0.2,0.5,paste(plotLegend[1]," ","(AUC = ",round(aucT1,digits = 3),")",sep=""),lwd = 3,bty="n",cex=1.4,col = plotColors[1]) 

for (i in 1:4) {
lines(unlist(aucOthers[[1]][i]),unlist(aucOthers[[2]][i]),col=plotColors[i+1],lwd=3)
}

legendY<-c(0.4,0.3,0.2,0.1)
for (i in 1:4) {
legendName<-paste(plotLegend[i+1]," ","(AUC = ",round(aucOthers[[3]][i],digits = 3),")",sep="")
legend(0.2,legendY[i],legendName,lwd = 3,bty="n",cex=1.4,col = plotColors[[i+1]]) 

}



####stage 4


##### take T1 enhancers as an example #####
## 1_stage 64_cell cross validation 
cv<-fread("midN1_model.cvpred.txt")
crossV<-cv[,c(2,3)]
## T1 enhancers prediction based on 2_stage 112_cell 
pos2<-fread("midN1fasta_bypr64cell.txt")
neg2<-fread("NS_midN1fasta_bypr64cell.txt")
## T1 enhancers prediction based on 3_stage late_gastrula 
pos3<-fread("midN1fasta_bypr112cell.txt")
neg3<-fread("NS_midN1fasta_bypr112cell.txt")
## T1 enhancers prediction based on 4_stage mid_neurula 
pos4<-fread("midN1fasta_byprlateG.txt")
neg4<-fread("NS_midN1fasta_byprlateG.txt")
## T1 enhancers prediction based on 5_stage mid_tailbud 
pos5<-fread("midN1fasta_bypref1.txt")
neg5<-fread("NS_midN1fasta_bypref1.txt")

## calculate AUC
# 1_stage 64_cell cross validation 
pred <- prediction(crossV$V2, crossV$V3) 
perf <- performance( pred, "tpr", "fpr" )
auc_temp <-performance( pred, measure = "auc")
aucT1<-unlist(slot(auc_temp, "y.values"))
# T1 enhancers prediction based on other models
pos<-cbind(pos2,pos3$V2,pos4$V2,pos5$V2)
neg<-cbind(neg2,neg3$V2,neg4$V2,neg5$V2)
stage1<-rocPrArea(pos,neg)
aucOthers<-list(stage1[[1]],stage1[[2]],stage1[[3]])

## plot
plotColors<-c(pal[3],pal[7],pal[4],pal[5],pal[2])
plotLegend<-c( "4_stage mid_neurula", "1_stage 64_cell", "2_stage 112_cell", "3_stage late_gastrula","5_stage mid_tailbud")

par(mar=c(5,5,2,2))
plot(unlist(perf@x.values),unlist(perf@y.values),lwd = 3,cex=1.5,col=plotColors[1],main="Ciona enhancers",xlab="False positive rate",ylab="Ture positive rate",
     cex.axis=1.5,cex.lab=1.5,cex.main=1.5,type="l")
legend(0.2,0.5,paste(plotLegend[1]," ","(AUC = ",round(aucT1,digits = 3),")",sep=""),lwd = 3,bty="n",cex=1.4,col = plotColors[1]) 

for (i in 1:4) {
  lines(unlist(aucOthers[[1]][i]),unlist(aucOthers[[2]][i]),col=plotColors[i+1],lwd=3)
}

legendY<-c(0.4,0.3,0.2,0.1)
for (i in 1:4) {
  legendName<-paste(plotLegend[i+1]," ","(AUC = ",round(aucOthers[[3]][i],digits = 3),")",sep="")
  legend(0.2,legendY[i],legendName,lwd = 3,bty="n",cex=1.4,col = plotColors[[i+1]]) 
  
}


####stage 5

##### take T1 enhancers as an example #####
## 1_stage 64_cell cross validation 
cv<-fread("ef1_model.cvpred.txt")
crossV<-cv[,c(2,3)]
## T1 enhancers prediction based on 2_stage 112_cell 
pos2<-fread("ef1fasta_bypr64cell.txt")
neg2<-fread("NS_ef1fasta_bypr64cell.txt")
## T1 enhancers prediction based on 3_stage late_gastrula 
pos3<-fread("ef1fasta_bypr112cell.txt")
neg3<-fread("NS_ef1fasta_bypr112cell.txt")
## T1 enhancers prediction based on 4_stage mid_neurula 
pos4<-fread("ef1fasta_byprlateG.txt")
neg4<-fread("NS_ef1fasta_byprlateG.txt")
## T1 enhancers prediction based on 5_stage mid_tailbud 
pos5<-fread("ef1fasta_byprmidN.txt")
neg5<-fread("NS_ef1fasta_byprmidN.txt")

## calculate AUC
# 1_stage 64_cell cross validation 
pred <- prediction(crossV$V2, crossV$V3) 
perf <- performance( pred, "tpr", "fpr" )
auc_temp <-performance( pred, measure = "auc")
aucT1<-unlist(slot(auc_temp, "y.values"))
# T1 enhancers prediction based on other models
pos<-cbind(pos2,pos3$V2,pos4$V2,pos5$V2)
neg<-cbind(neg2,neg3$V2,neg4$V2,neg5$V2)
stage1<-rocPrArea(pos,neg)
aucOthers<-list(stage1[[1]],stage1[[2]],stage1[[3]])

## plot
plotColors<-c(pal[3],pal[7],pal[4],pal[5],pal[2])
plotLegend<-c( "5_stage mid_tailbud", "1_stage 64_cell", "2_stage 112_cell", "3_stage late_gastrula", "4_stage mid_neurula")

par(mar=c(5,5,2,2))
plot(unlist(perf@x.values),unlist(perf@y.values),lwd = 3,cex=1.5,col=plotColors[1],main="Ciona enhancers",xlab="False positive rate",ylab="Ture positive rate",
     cex.axis=1.5,cex.lab=1.5,cex.main=1.5,type="l")
legend(0.2,0.5,paste(plotLegend[1]," ","(AUC = ",round(aucT1,digits = 3),")",sep=""),lwd = 3,bty="n",cex=1.4,col = plotColors[1]) 

for (i in 1:4) {
  lines(unlist(aucOthers[[1]][i]),unlist(aucOthers[[2]][i]),col=plotColors[i+1],lwd=3)
}

legendY<-c(0.4,0.3,0.2,0.1)
for (i in 1:4) {
  legendName<-paste(plotLegend[i+1]," ","(AUC = ",round(aucOthers[[3]][i],digits = 3),")",sep="")
  legend(0.2,legendY[i],legendName,lwd = 3,bty="n",cex=1.4,col = plotColors[[i+1]]) 
  
}

dev.off()

##########################################################################################################################
#####* deltaSVM summary *#####
##### function used for deltaSVM analysis #####
setwd("deltaSVM/")


dataMod<-function(deltaSVM) {
## change the first column into bed format
deltaSVM<-as.data.frame(deltaSVM)
names(deltaSVM)<-c("ID","deltaSVM","varNumb","pValue")
#splString<-strsplit(deltaSVM$ID,"_",fixed=TRUE)
#splString<-data.frame(unlist(splString))
#ID.bed<-matrix(splString$unlist.splString., ncol=5, byrow=TRUE)
#ID.bed<-ID.bed[,c(1:3)]
#ID.bed[,2]<-as.numeric(ID.bed[,2])-1
#deltaSVM<-cbind(ID.bed, deltaSVM[,c(2:4)])
#colnames(deltaSVM)[c(1:3)]<-c("chr","start","end")
#deltaSVM$start<-as.numeric(as.character(deltaSVM$start))
#deltaSVM$end<-as.numeric(as.character(deltaSVM$end))

## qvalue
deltaSVM<-deltaSVM[order(deltaSVM$pValue),]
deltaSVM$qValue<-qvalue(deltaSVM$pValue,pi0=1)$qvalues
return(deltaSVM)
}

pos_sel_64<-fread("deltaSVM/64_pos_sel_cinte_crobu1k")
pos_sel_64<-dataMod(pos_sel_64)

pos_sel_112<-fread("deltaSVM/112_pos_sel_cinte")
pos_sel_112<-dataMod(pos_sel_112)

lateG_pos_sel<-fread("deltaSVM/lateG_pos_sel_cinte")
lateG_pos_sel<-dataMod(lateG_pos_sel)

midN1_pos_sel<-fread("deltaSVM/midN1_pos_sel_cinte")
midN1_pos_sel<-dataMod(midN1_pos_sel)

ef1_pos_sel<-fread("deltaSVM/ef1_pos_sel")
ef1_pos_sel<-dataMod(ef1_pos_sel)

#pos_sel_112  pos_sel_64  ef1_pos_sel  lateG_pos_sel  midN1_pos_sel

## plot deltaSVM and pvalue
par(mfrow=c(1,2))
par(mar=c(7,5,4,2))
hist(pos_sel_64$deltaSVM,breaks = 50,main="C. intestinalis 64 cell stage",xlab="deltaSVM",xlim=c(-500,500),
cex.lab=2,cex.axis=2,cex.main=2,col=pal[3])
hist(pos_sel_64$pValue,breaks = 50,main="C. intestinalis  64 cell stage",xlab="Pvalue",xlim=c(0,1),
cex.lab=2,cex.axis=2,cex.main=2,col=pal[3])

hist(pos_sel_112$deltaSVM,breaks = 50,main="Cion robusta 112 cell stage",xlab="deltaSVM",xlim=c(-50,50),
cex.lab=2,cex.axis=2,cex.main=2,col=pal[7])
hist(pos_sel_112$pValue,breaks = 50,main="Cion robusta 112 cell stage",xlab="Pvalue",xlim=c(0,1),
cex.lab=2,cex.axis=2,cex.main=2,col=pal[7])

hist(lateG_pos_sel$deltaSVM,breaks = 50,main="Cion robusta Late Gastrula",xlab="deltaSVM",xlim=c(-50,50),
cex.lab=2,cex.axis=2,cex.main=2,col=pal[4])
hist(lateG_pos_sel$pValue,breaks = 50,main="Cion robusta Late Gastrula",xlab="Pvalue",xlim=c(0,1),
cex.lab=2,cex.axis=2,cex.main=2,col=pal[4])

hist(midN1_pos_sel$deltaSVM,breaks = 50,main="Cion robusta Mid Neurula",xlab="deltaSVM",xlim=c(-50,50),
cex.lab=2,cex.axis=2,cex.main=2,col=pal[5])
hist(midN1_pos_sel$pValue,breaks = 50,main="Cion robusta Mid Neurula",xlab="Pvalue",xlim=c(0,1),
cex.lab=2,cex.axis=2,cex.main=2,col=pal[5])

hist(ef1_pos_sel$deltaSVM,breaks = 50,main="C. intestinalis Mid Tailbud",xlab="deltaSVM",xlim=c(-50,50),
cex.lab=2,cex.axis=2,cex.main=2,col=pal[2])
hist(ef1_pos_sel$pValue,breaks = 50,main="C. intestinalis Mid Tailbud",xlab="Pvalue",xlim=c(0,1),
cex.lab=2,cex.axis=2,cex.main=2,col=pal[2])

####################################################################################################