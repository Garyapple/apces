############## INIT SYSTEM ################
rm(list = ls(all = TRUE)) 
setwd("F:/Classification")
workpath="F:/Classification/"

######### Ensemble Feature Selection #########
source("EGSG.CV.R")

datapath<-paste(workpath,"Datasets/",sep="")
datasets=list.files(datapath)

outputpath<-paste(workpath,"EGSG",sep="")
if(file.exists(outputpath)==FALSE)dir.create(outputpath)

sn.res=matrix(0,10,7)
sp.res=matrix(0,10,7)
acc.res=matrix(0,10,7)
gmean.res=matrix(0,10,7)

for(k in 1:length(datasets)){
  gc(FALSE,TRUE) 
  msg<-paste("Processing the dataset:  ",datasets[k],sep="")
  print(msg)
  
  filespath<-paste(datapath,datasets[k],sep="")
  expfile<-paste(filespath,"/dataset(normalized).csv",sep="") # When the train set is provided directly
  dataset<-read.csv(expfile,TRUE)
  classL<-dataset[,ncol(dataset)]
  data<-dataset[,-ncol(dataset)]
  data=as.matrix(data)
  rm(dataset)
  tabL=table(classL)
  classL<- factor(ifelse(classL==names(tabL[1]), 'control', 'treatment'))
  
  expfile<-paste(filespath,"/dataset(discretized).csv",sep="") # When the train set is provided directly
  dataset<-read.csv(expfile,TRUE)
  disc.data<-dataset[,-ncol(dataset)]
  disc.data=as.matrix(disc.data)
  rm(dataset)
  
  for(i in 1:10){
    tmpres=EGSG.CV(data=data,disc.data=disc.data,folds=5,
                   classlabels=classL,nFit=21)
    sn.res[i,k]=tmpres$sn
    sp.res[i,k]=tmpres$sp
    acc.res[i,k]=tmpres$acc
    gmean.res[i,k]=tmpres$gmean
  }
  
}

means=apply(sn.res,2,mean)
sn.res=rbind(sn.res,means)
rownames(sn.res)=c(1:10,"mean")
expfile<-paste(outputpath,"/","5CV_sn(svm_thr).csv",sep="")
write.csv(sn.res,expfile)

means=apply(sp.res,2,mean)
sp.res=rbind(sp.res,means)
rownames(sp.res)=c(1:10,"mean")
expfile<-paste(outputpath,"/","5CV_sp(svm_thr).csv",sep="")
write.csv(sp.res,expfile)

means=apply(acc.res,2,mean)
acc.res=rbind(acc.res,means)
rownames(acc.res)=c(1:10,"mean")
expfile<-paste(outputpath,"/","5CV_acc(svm_thr).csv",sep="")
write.csv(acc.res,expfile)

means=apply(gmean.res,2,mean)
gmean.res=rbind(gmean.res,means)
rownames(gmean.res)=c(1:10,"mean")
expfile<-paste(outputpath,"/","5CV_gmean(svm_thr).csv",sep="")
write.csv(gmean.res,expfile)
