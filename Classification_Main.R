############## INIT SYSTEM ################
rm(list = ls(all = TRUE)) 
setwd("F:/Classification")
workpath="F:/Classification/"

######### Ensemble Feature Selection #########
source("APBEFS.CV.R")

datapath<-paste(workpath,"Datasets/",sep="")
datasets=list.files(datapath)

outputpath<-paste(workpath,"APBEFS",sep="")
if(file.exists(outputpath)==FALSE)dir.create(outputpath)

for(k in 1:7){
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
  
  range=seq(1,50,by=4)
  sn.res=matrix(0,10,length(range))
  sp.res=matrix(0,10,length(range))
  acc.res=matrix(0,10,length(range))
  gmean.res=matrix(0,10,length(range))
  
  count=1
  for(nFit in range){
    sn.tmpres=array(0,10)
    sp.tmpres=array(0,10)
    acc.tmpres=array(0,10)
    gmean.tmpres=array(0,10)
    for(i in 1:10){
      tmpres=APBEFS.CV(data=data,classlabels=classL,folds=5,nFit=nFit,topK=300)
      sn.tmpres[i]=tmpres$sn
      sp.tmpres[i]=tmpres$sp
      acc.tmpres[i] = tmpres$acc
      gmean.tmpres[i] = tmpres$gmean
    }
    sn.res[,count]=sn.tmpres
    sp.res[,count]=sp.tmpres
    acc.res[,count]=acc.tmpres
    gmean.res[,count]=gmean.tmpres
    count=count+1
  }
  means=apply(sn.res,2,mean)
  sn.res=rbind(sn.res,means)
  rownames(sn.res)=c(1:10,"mean")
  colnames(sn.res)=c(range)
  expfile<-paste(outputpath,"/",substr(datasets[k],1,3),"_5CV_sn(svm_thr).csv",sep="")
  write.csv(sn.res,expfile)
  
  means=apply(sp.res,2,mean)
  sp.res=rbind(sp.res,means)
  rownames(sp.res)=c(1:10,"mean")
  colnames(sp.res)=c(range)
  expfile<-paste(outputpath,"/",substr(datasets[k],1,3),"_5CV_sp(svm_thr).csv",sep="")
  write.csv(sp.res,expfile)
  
  means=apply(acc.res,2,mean)
  acc.res=rbind(acc.res,means)
  rownames(acc.res)=c(1:10,"mean")
  colnames(acc.res)=c(range)
  expfile<-paste(outputpath,"/",substr(datasets[k],1,3),"_5CV_acc(svm_thr).csv",sep="")
  write.csv(acc.res,expfile)
  
  means=apply(gmean.res,2,mean)
  gmean.res=rbind(gmean.res,means)
  rownames(gmean.res)=c(1:10,"mean")
  colnames(gmean.res)=c(range)
  expfile<-paste(outputpath,"/",substr(datasets[k],1,3),"_5CV_gmean(svm_thr).csv",sep="")
  write.csv(gmean.res,expfile)
  
}
