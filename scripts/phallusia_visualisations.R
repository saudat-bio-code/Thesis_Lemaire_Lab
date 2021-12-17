
library("data.table")
library("RColorBrewer") 
library("plyr")
library("qvalue")
library("topGO")
library("ROCR")
library("tseries")
setwd("/home/salishayeva/Downloads/bed_files/Pmamilata/ssh_imported_final/predictions_with_phami_neg_set/")
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3")) 



rocPrArea2<-function(positive, negative) {
positive$state<-rep(1,tiems=nrow(positive))
negative$state<-rep(0,tiems=nrow(negative))
aucResult<-c()
pr_plotValueX<-list()
pr_plotValueY<-list()
for (i in 2:4) {
pos<-cbind(positive[[i]],positive[,5])
neg<-cbind(negative[[i]],negative[,5])
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
## T1 enhancers prediction based on 4_stage mid_neurula_uncertain 
#pos5<-fread("64cellfasta_bypref1.txt")
#neg5<-fread("NS_64cellfasta_bypref1.txt")

## calculate AUC
# 1_stage 64_cell cross validation 
pred <- prediction(crossV$V2, crossV$V3) 
perf <- performance( pred, "tpr", "fpr" )
auc_temp <-performance( pred, measure = "auc")
aucT1<-unlist(slot(auc_temp, "y.values"))
# T1 enhancers prediction based on other models
pos<-cbind(pos2,pos3$V2,pos4$V2)
neg<-cbind(neg2,neg3$V2,neg4$V2)
stage1<-rocPrArea2(pos,neg)
aucOthers<-list(stage1[[1]],stage1[[2]],stage1[[3]])

## plot
plotColors<-c(pal[3],pal[7],pal[4],pal[5])
plotLegend<-c("1_stage 64_cell", "2_stage 112_cell", "3_stage late_gastrula", "4_stage mid_neurula")

par(mar=c(5,5,2,2))
plot(unlist(perf@x.values),unlist(perf@y.values),lwd = 3,cex=1.5,col=plotColors[1],main="Phallusia enhancers",xlab="False positive rate",ylab="Ture positive rate",
     cex.axis=1.5,cex.lab=1.5,cex.main=1.5,type="l")
legend(0.2,0.5,paste(plotLegend[1]," ","(AUC = ",round(aucT1,digits = 3),")",sep=""),lwd = 2,bty="n",cex=1.4,col = plotColors[1]) 

for (i in 1:3) {
  lines(unlist(aucOthers[[1]][i]),unlist(aucOthers[[2]][i]),col=plotColors[i+1],lwd=3)
}

legendY<-c(0.4,0.3,0.2,0.1)
for (i in 1:3) {
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
## T1 enhancers prediction based on 4_stage mid_neurula_uncertain 
#pos5<-fread("112cellfasta_bypref1.txt")
#neg5<-fread("NS_112cellfasta_bypref1.txt")

## calculate AUC
# 1_stage 64_cell cross validation 
pred <- prediction(crossV$V2, crossV$V3) 
perf <- performance( pred, "tpr", "fpr" )
auc_temp <-performance( pred, measure = "auc")
aucT1<-unlist(slot(auc_temp, "y.values"))
# T1 enhancers prediction based on other models
pos<-cbind(pos2,pos3$V2,pos4$V2)
neg<-cbind(neg2,neg3$V2,neg4$V2)
stage1<-rocPrArea2(pos,neg)
aucOthers<-list(stage1[[1]],stage1[[2]],stage1[[3]])

## plot
plotColors<-c(pal[3],pal[7],pal[4],pal[5])
plotLegend<-c("2_stage 112_cell","1_stage 64_cell", "3_stage late_gastrula", "4_stage mid_neurula")

par(mar=c(5,5,2,2))
plot(unlist(perf@x.values),unlist(perf@y.values),lwd = 3,cex=1.5,col=plotColors[1],main="Phallusia enhancers",xlab="False positive rate",ylab="Ture positive rate",
     cex.axis=1.5,cex.lab=1.5,cex.main=1.5,type="l")
legend(0.2,0.5,paste(plotLegend[1]," ","(AUC = ",round(aucT1,digits = 3),")",sep=""),lwd = 3,bty="n",cex=1.4,col = plotColors[1]) 

for (i in 1:3) {
  lines(unlist(aucOthers[[1]][i]),unlist(aucOthers[[2]][i]),col=plotColors[i+1],lwd=3)
}

legendY<-c(0.4,0.3,0.2,0.1)
for (i in 1:3) {
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
## T1 enhancers prediction based on 4_stage mid_neurula_uncertain 
#pos5<-fread("lateGfasta_bypref1.txt")
#neg5<-fread("NS_lateGfasta_bypref1.txt")

## calculate AUC
# 1_stage 64_cell cross validation 
pred <- prediction(crossV$V2, crossV$V3) 
perf <- performance( pred, "tpr", "fpr" )
auc_temp <-performance( pred, measure = "auc")
aucT1<-unlist(slot(auc_temp, "y.values"))
# T1 enhancers prediction based on other models
pos<-cbind(pos2,pos3$V2,pos4$V2)
neg<-cbind(neg2,neg3$V2,neg4$V2)
stage1<-rocPrArea2(pos,neg)
aucOthers<-list(stage1[[1]],stage1[[2]],stage1[[3]])

## plot
plotColors<-c(pal[3],pal[7],pal[4],pal[5])
plotLegend<-c("3_stage late_gastrula", "1_stage 64_cell", "2_stage 112_cell", "4_stage mid_neurula")

par(mar=c(5,5,2,2))
plot(unlist(perf@x.values),unlist(perf@y.values),lwd = 3,cex=1.5,col=plotColors[1],main="Phallusia enhancers",xlab="False positive rate",ylab="Ture positive rate",
     cex.axis=1.5,cex.lab=1.5,cex.main=1.5,type="l")
legend(0.2,0.5,paste(plotLegend[1]," ","(AUC = ",round(aucT1,digits = 3),")",sep=""),lwd = 3,bty="n",cex=1.4,col = plotColors[1]) 

for (i in 1:3) {
  lines(unlist(aucOthers[[1]][i]),unlist(aucOthers[[2]][i]),col=plotColors[i+1],lwd=3)
}

legendY<-c(0.4,0.3,0.2,0.1)
for (i in 1:3) {
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
## T1 enhancers prediction based on 4_stage mid_neurula_uncertain 
#pos5<-fread("midN1fasta_bypref1.txt")
#neg5<-fread("NS_midN1fasta_bypref1.txt")

## calculate AUC
# 1_stage 64_cell cross validation 
pred <- prediction(crossV$V2, crossV$V3) 
perf <- performance( pred, "tpr", "fpr" )
auc_temp <-performance( pred, measure = "auc")
aucT1<-unlist(slot(auc_temp, "y.values"))
# T1 enhancers prediction based on other models
pos<-cbind(pos2,pos3$V2,pos4$V2)
neg<-cbind(neg2,neg3$V2,neg4$V2)
stage1<-rocPrArea2(pos,neg)
aucOthers<-list(stage1[[1]],stage1[[2]],stage1[[3]])

## plot
plotColors<-c(pal[3],pal[7],pal[4],pal[5])
plotLegend<-c( "4_stage mid_neurula", "1_stage 64_cell", "2_stage 112_cell", "3_stage late_gastrula")

par(mar=c(5,5,2,2))
plot(unlist(perf@x.values),unlist(perf@y.values),lwd = 3,cex=1.5,col=plotColors[1],main="Phallusia enhancers",xlab="False positive rate",ylab="Ture positive rate",
     cex.axis=1.5,cex.lab=1.5,cex.main=1.5,type="l")
legend(0.2,0.5,paste(plotLegend[1]," ","(AUC = ",round(aucT1,digits = 3),")",sep=""),lwd = 3,bty="n",cex=1.4,col = plotColors[1]) 

for (i in 1:3) {
  lines(unlist(aucOthers[[1]][i]),unlist(aucOthers[[2]][i]),col=plotColors[i+1],lwd=3)
}

legendY<-c(0.4,0.3,0.2,0.1)
for (i in 1:3) {
  legendName<-paste(plotLegend[i+1]," ","(AUC = ",round(aucOthers[[3]][i],digits = 3),")",sep="")
  legend(0.2,legendY[i],legendName,lwd = 3,bty="n",cex=1.4,col = plotColors[[i+1]]) 
  
}


dev.off()

