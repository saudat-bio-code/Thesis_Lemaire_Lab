

#####** positive selection **#####
#####* ROC analysis for SVM models *#####
##### function used for ROC analysis #####
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

################################################################################################
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




##### take T1 enhancers as an example #####
## T1 model cross validation 
cv<-fread("dmel2_3_model.cvpred.txt")
crossV<-cv[,c(2,3)]
## T4 enhancers prediction based on T1 model 
pos2<-fread("2_3fasta_bypr4_5.txt")
neg2<-fread("NS_2_3fasta_bypr4_5.txt")
## T4 enhancers prediction based on T2 model 
pos3<-fread("2_3fasta_bypr6_7.txt")
neg3<-fread("NS_2_3fasta_bypr6_7.txt")
## T4 enhancers prediction based on T3 model 
pos4<-fread("2_3fasta_byprmidN.txt")
neg4<-fread("NS_2_3fasta_byprmidN.txt")

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

stage1<-rocPrArea(pos,neg)
aucOthers<-list(stage1[[1]],stage1[[2]],stage1[[3]])

## plot
plotColors<-c(pal[3],pal[4],pal[5],pal[2])
plotLegend<-c("Late Gastrula",  "64_cell model", " 112_cell model", " Mid Neurula")


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




























#####* deltaSVM summary *#####
##### function used for deltaSVM analysis #####
dataMod<-function(deltaSVM) {
  ## change the first column into bed format
  deltaSVM<-as.data.frame(deltaSVM)
  names(deltaSVM)<-c("ID","deltaSVM","varNumb","pValue")
  splString<-strsplit(deltaSVM$ID,"_",fixed=TRUE)
  splString<-data.frame(unlist(splString))
  ID.bed<-matrix(splString$unlist.splString., ncol=5, byrow=TRUE)
  ID.bed<-ID.bed[,c(1:3)]
  ID.bed[,2]<-as.numeric(ID.bed[,2])-1
  deltaSVM<-cbind(ID.bed, deltaSVM[,c(2:4)])
  colnames(deltaSVM)[c(1:3)]<-c("chr","start","end")
  deltaSVM$start<-as.numeric(as.character(deltaSVM$start))
  deltaSVM$end<-as.numeric(as.character(deltaSVM$end))
  
  ## qvalue
  deltaSVM<-deltaSVM[order(deltaSVM$pValue),]
  deltaSVM$qValue<-qvalue(deltaSVM$pValue,pi0=1)$qvalues
  return(deltaSVM)
}

dmelT1Enh_deltaSVM<-fread("deltaSVM/dmel24_specific_enhancer_dm3_deltaSVM_highertailTest.txt")
dmelT1Enh_deltaSVM<-dataMod(dmelT1Enh_deltaSVM)

## plot deltaSVM and pvalue
par(mfrow=c(1,2))
par(mar=c(7,5,4,2))
hist(dmelT1Enh_deltaSVM$deltaSVM,breaks = 50,main="D.melanogaster T1",xlab="deltaSVM",xlim=c(-50,50),
     cex.lab=2,cex.axis=2,cex.main=2,col=pal[3])
hist(dmelT1Enh_deltaSVM$pValue,breaks = 50,main="D.melanogaster T1",xlab="Pvalue",xlim=c(0,1),
     cex.lab=2,cex.axis=2,cex.main=2,col=pal[3])
