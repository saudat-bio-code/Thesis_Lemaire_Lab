
##### take T1 enhancers as an example #####
## T1 model cross validation 
cv<-fread("dmel4_5_model.cvpred.txt")
crossV<-cv[,c(2,3)]
## T4 enhancers prediction based on T1 model 
pos2<-fread("4_5fasta_bypr6_7.txt")
neg2<-fread("NS_4_5fasta_bypr6_7.txt")
## T4 enhancers prediction based on T2 model 
pos3<-fread("4_5fasta_byprmidN.txt")
neg3<-fread("NS_4_5fasta_byprmidN.txt")
## T4 enhancers prediction based on T3 model 
pos4<-fread("4_5fasta_bypr2_3.txt")
neg4<-fread("NS_4_5fasta_bypr2_3.txt")

## calculate AUC
# T1 model cross validation 
pred <- prediction(crossV$V2, crossV$V3) 
perf <- performance( pred, "tpr", "fpr" )
auc_temp <-performance( pred, measure = "auc")
aucT1<-unlist(slot(auc_temp, "y.values"))
# T1 enhancers prediction based on other models

# mode 1 2 3
pos<-cbind(pos2,pos3$V2,pos4$V2)
neg<-cbind(neg2,neg3$V2,neg4$V2)

stage1<-rocPrArea2(pos,neg)
aucOthers<-list(stage1[[1]],stage1[[2]],stage1[[3]])

## plot
plotColors<-c(pal[3],pal[4],pal[5],pal[2])
plotLegend<-c( "64_cell model", " 112_cell model", " Mid Neurula", "Late Gastrula")


par(mar=c(5,5,2,2))
plot(unlist(perf@x.values),unlist(perf@y.values),lwd = 3,cex=1.5,col=plotColors[1],main="First stage enhancers",xlab="False positive rate",ylab="Ture positive rate",
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


#############################################################################################################################
#  6-7 enhancers


##### take T1 enhancers as an example #####
## T1 model cross validation 
cv<-fread("dmel6_7_model.cvpred.txt")
crossV<-cv[,c(2,3)]
## T4 enhancers prediction based on T1 model 
pos2<-fread("6_7fasta_bypr4_5.txt")
neg2<-fread("NS_6_7fasta_bypr4_5.txt")
## T4 enhancers prediction based on T2 model 
pos3<-fread("6_7fasta_byprmidN.txt")
neg3<-fread("NS_6_7fasta_byprmidN.txt")
## T4 enhancers prediction based on T3 model 
pos4<-fread("6_7fasta_bypr2_3.txt")
neg4<-fread("NS_6_7fasta_bypr2_3.txt")

## calculate AUC
# T1 model cross validation 
pred <- prediction(crossV$V2, crossV$V3) 
perf <- performance( pred, "tpr", "fpr" )
auc_temp <-performance( pred, measure = "auc")
aucT1<-unlist(slot(auc_temp, "y.values"))
# T1 enhancers prediction based on other models

# mode 1 2 3
pos<-cbind(pos2,pos3$V2,pos4$V2)
neg<-cbind(neg2,neg3$V2,neg4$V2)

stage1<-rocPrArea2(pos,neg)
aucOthers<-list(stage1[[1]],stage1[[2]],stage1[[3]])

## plot
plotColors<-c(pal[3],pal[4],pal[5],pal[2])
plotLegend<-c( " 112_cell model","64_cell model", " Mid Neurula", "Late Gastrula")


par(mar=c(5,5,2,2))
plot(unlist(perf@x.values),unlist(perf@y.values),lwd = 3,cex=1.5,col=plotColors[1],main="First stage enhancers",xlab="False positive rate",ylab="Ture positive rate",
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

#############################################################################################################################
# Midneurula nehancer


##### take T1 enhancers as an example #####
## T1 model cross validation 
cv<-fread("midN1_model.cvpred.txt")
crossV<-cv[,c(2,3)]
## T4 enhancers prediction based on T1 model 
pos2<-fread("midN1fasta_bypr4_5.txt")
neg2<-fread("NS_midN1fasta_bypr4_5.txt")
## T4 enhancers prediction based on T2 model 
pos3<-fread("midN1fasta_bypr6_7.txt")
neg3<-fread("NS_midN1fasta_bypr6_7.txt")
## T4 enhancers prediction based on T3 model 
pos4<-fread("midN1fasta_bypr2_3.txt")
neg4<-fread("NS_midN1fasta_bypr6_7.txt")

## calculate AUC
# T1 model cross validation 
pred <- prediction(crossV$V2, crossV$V3) 
perf <- performance( pred, "tpr", "fpr" )
auc_temp <-performance( pred, measure = "auc")
aucT1<-unlist(slot(auc_temp, "y.values"))
# T1 enhancers prediction based on other models

# mode 1 2 3
pos<-cbind(pos2,pos3$V2,pos4$V2)
neg<-cbind(neg2,neg3$V2,neg4$V2)

stage1<-rocPrArea2(pos,neg)
aucOthers<-list(stage1[[1]],stage1[[2]],stage1[[3]])

## plot
plotColors<-c(pal[3],pal[4],pal[5],pal[2])
plotLegend<-c( "Mid Neurula model", " 64_cell model"," 112_cell model",  "Late Gastrula")


par(mar=c(5,5,2,2))
plot(unlist(perf@x.values),unlist(perf@y.values),lwd = 3,cex=1.5,col=plotColors[1],main="First stage enhancers",xlab="False positive rate",ylab="Ture positive rate",
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



