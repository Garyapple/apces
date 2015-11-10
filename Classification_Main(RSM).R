############## INIT SYSTEM ################
rm(list = ls(all = TRUE)) 
setwd("F:/Classification")
workpath="F:/Classification/"

require(caret,quietly=TRUE)
require(e1071)

datapath<-paste(workpath,"Datasets/",sep="")
datasets=list.files(datapath)

outputpath<-paste(workpath,"RSM",sep="")
if(file.exists(outputpath)==FALSE)dir.create(outputpath)

feature.num =256
nFit=21 

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
  
  sn.res=matrix(0,10,folds)
  sp.res=matrix(0,10,folds)
  acc.res=matrix(0,10,folds)
  gmean.res=matrix(0,10,folds)
  
  for(times in 1:10){
    cv.index=createFolds(classL,k=folds)
    
    count =1
    for (fold in 1:folds) {
      print(paste("Fold:",fold))
      test <- cv.index[[fold]]
      traindata=data[-test,]
      trainlabels=classL[-test]
      testdata = data[test,]
      testlabels = classL[test]
      
      N=length(trainlabels)
      control.num=sum(as.numeric(trainlabels=="control"))
      treatment.num=N-control.num
      
      # Random select features
      len=ncol(traindata)
      feature.smp<-matrix(0,nFit,feature.num)
      for(j in 1:nFit){
        feature.smp[j,]=sample(1:len,feature.num,replace=FALSE) # random select feature sets  
      }
      
      # Train base classifies
      threshold = (control.num - treatment.num)/(N+2) #a=1
      pred.allres = matrix(0,nFit,length(testlabels))
      
      for(idx in 1:nFit){
        print(idx)
        features=feature.smp[idx,]
        fit <- svm(classL~.,data.frame(data= traindata[,features],classL=trainlabels)) # save the tunne model
        fit.predict<-predict(fit,data.frame(data=testdata[,features]),decision.values = TRUE)
        pred.allres[idx,]=ifelse(attr(fit.predict,"decision.values") > threshold,1,0)
      }
      
      predict.votes = apply(pred.allres,2,sum)
      fit.predict<-factor(ifelse(predict.votes > nFit/2,"control","treatment"))
      
      if(nlevels(fit.predict)>1){
        tmpres = confusionMatrix(fit.predict,testlabels)
        sp.res[fold,count]=tmpres$byClass["Specificity"]
        sn.res[fold,count]=tmpres$byClass["Sensitivity"]
        gmean.res[fold,count] <- sqrt(sp.res[fold,count]*sn.res[fold,count])
        acc.res[fold,count] <- tmpres$overall["Accuracy"]
      }
      else{
        if(fit.predict[1]=="control"){
          sn.res[fold,count]=1
          sp.res[fold,count]=0
        }
        else{
          sn.res[fold,count]=0
          sp.res[fold,count]=1    
        }
        gmean.res[fold,count] <- 0
        acc.res[fold,count] <- length(which(as.character(fit.predict)==as.character(testlabels)))/length(test)
      }
      count=count+1
    }
  }
  
  means=apply(sn.res,1,mean)
  sn.res=cbind(sn.res,means)
  rownames(sn.res)=c(1:10)
  colnames(sn.res)=c(1:5,"mean")
  expfile<-paste(outputpath,"/",substr(datasets[k],1,3),"_5CV_sn(svm_thr).csv",sep="")
  write.csv(sn.res,expfile)
  
  means=apply(sp.res,1,mean)
  sp.res=rbind(sp.res,means)
  rownames(sp.res)=c(1:10)
  colnames(sn.res)=c(1:5,"mean")
  expfile<-paste(outputpath,"/",substr(datasets[k],1,3),"_5CV_sp(svm_thr).csv",sep="")
  write.csv(sp.res,expfile)
  
  means=apply(acc.res,1,mean)
  acc.res=rbind(acc.res,means)
  rownames(acc.res)=c(1:10)
  colnames(sn.res)=c(1:5,"mean")
  expfile<-paste(outputpath,"/",substr(datasets[k],1,3),"_5CV_acc(svm_thr).csv",sep="")
  write.csv(acc.res,expfile)
  
  means=apply(gmean.res,1,mean)
  gmean.res=rbind(gmean.res,means)
  rownames(gmean.res)=c(1:10)
  colnames(sn.res)=c(1:5,"mean")
  expfile<-paste(outputpath,"/",substr(datasets[k],1,3),"_5CV_gmean(svm_thr).csv",sep="")
  write.csv(gmean.res,expfile)
  
}
