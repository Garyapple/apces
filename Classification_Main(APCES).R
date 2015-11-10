############## INIT SYSTEM ################
rm(list = ls(all = TRUE)) 
setwd("F:/Classification")
workpath="F:/Classification/"

require(apcluster,quietly=TRUE)
require(WGCNA,quietly=TRUE)
require(caret,quietly=TRUE)
require(e1071)
require(caTools)
require(GeneSelector)

######### Ensemble Feature Selection #########
source("APCES.CV.R")

datapath<-paste(workpath,"Datasets/",sep="")
datasets=list.files(datapath)

outputpath<-paste(workpath,"APCES",sep="")
if(file.exists(outputpath)==FALSE)dir.create(outputpath)

folds=5
top_K=500
nFit=45
# ranges=seq(1,nFit,by=4)

for(k in 7:7){
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
  
  sn.res=matrix(0,10,1)
  sp.res=matrix(0,10,1)
  acc.res=matrix(0,10,1)
  gmean.res=matrix(0,10,1)
  
  for(times in 1:10){
    cv.index=createFolds(classL,k=folds)            #把类标签分成5等份
    sn.tmpres=matrix(0,folds,1)
    sp.tmpres=matrix(0,folds,1)
    acc.tmpres=matrix(0,folds,1)
    gmean.tmpres=matrix(0,folds,1)
    
    for (fold in 1:folds) {
      print(paste("Fold:",fold))
      test <- cv.index[[fold]]
      testdata = data[test,]
      testlabels = classL[test]
      if(fold==1)
      {
        check <- cv.index[[fold+1]]
        checkdata = data[check,]
        checklabels = classL[check]
        traindata = data[c(-test,-check),]
        trainlabels = classL[c(-test,-check)]
      }
      else
      {
        check <- cv.index[[fold-1]]
        checkdata = data[check,]
        checklabels = classL[check]
        traindata = data[c(-test,-check),]
        trainlabels = classL[c(-test,-check)]
      }
      
      # Aggregation 
      N=length(trainlabels)
      control.num=sum(as.numeric(trainlabels=="control"))
      treatment.num=N-control.num
      
      if(control.num<treatment.num){
        flag="treatment";
        basedata=traindata[which(trainlabels=="control"),]
        minnum=control.num
      }
      else{
        flag="control";
        basedata=traindata[which(trainlabels=="treatment"),]
        minnum=treatment.num
      }
      
      repeats=21 # rank times
      AggRanks=c()
      for(iter in 1:repeats){
        if(flag=="treatment"){
          index=which(trainlabels=="treatment")
          data.smp=sample(1:treatment.num,minnum)
          classLab=c(rep(0,control.num),rep(1,control.num))
        }
        else{
          index=which(trainlabels=="control")
          data.smp=sample(1:control.num,minnum)
          classLab=c(rep(1,treatment.num),rep(0,treatment.num))
        }
        data.add=data[index[data.smp],]
        data.tmp=rbind(basedata,data.add)
        rm(data.add,data.smp,index)
        
        # Aggregation 
        eBayes.res<-RankingWelchT(t(data.tmp),classLab)
        eBayes.res = sort(eBayes.res@ranking,index.return=TRUE)$ix[1:top_K]
        if(iter==1){
          AggRanks=eBayes.res
        }
        else{
          AggRanks=rbind(AggRanks,eBayes.res)
        }
        rm(data.tmp,eBayes.res)
      }
      AggRanks = as.vector(AggRanks)
      feature.topK=unique(AggRanks)
      
      rm(AggRanks)
      
      # AP Clustering
      control.samples = which(trainlabels=="control")
      treatment.samples =  which(trainlabels=="treatment")
      
      dis1=adjacency(traindata[control.samples,feature.topK],type = "signed",corFnc = "bicor",
                     power=2,corOptions = "maxPOutliers = 0.02")
      dis2=adjacency(traindata[treatment.samples,feature.topK],type = "signed",corFnc = "bicor",
                     power=2,corOptions = "maxPOutliers = 0.02")
      dis=(dis1+dis2)/2-1
      apres<-apcluster(s=dis)
      exemplarIndex<-apres@exemplars
      K=length(exemplarIndex)
      
      # Random select features
      feature.smp<-matrix(0,nFit,K)
      if(nFit==1){
        feature.smp[1,]=feature.topK[exemplarIndex]
      }
      else{
        feature.smp[1,]=feature.topK[exemplarIndex] # use  the center of feature group
        for(j in 1:K){
          cluster=apres[[j]]
          len=length(cluster)
          tmp=sample(1:len,(nFit-1),replace=TRUE) # random select feature sets  
          cluster=feature.topK[cluster]
          feature.smp[2:nFit,j]=cluster[tmp]
        }
      }
      
      # Train base classifies
      threshold = (control.num - treatment.num)/(N+2) #a=1
      pred.checkallres = matrix(0,nFit,length(checklabels))
      pred.allres = matrix(0,nFit,length(testlabels))
      j<-1
      for(idx in 1:nFit){
        print(idx)
        features=feature.smp[idx,]
        fit <- svm(classL~.,data.frame(data= traindata[,features],classL=trainlabels)) # save the tunne model
        fit.predict1<-predict(fit,data.frame(data=checkdata[,features]),decision.values = TRUE)
        pred.checkallres[idx,]=ifelse(attr(fit.predict1,"decision.values") > threshold,1,0)
        checklabels <- factor(ifelse(classL[check]=="control",1,0))
        Caa = sum(as.numeric((pred.checkallres[idx,] == checklabels) & (checklabels == 1))) #计算kappa系数的相关参数
        Can = sum(as.numeric((pred.checkallres[idx,] != checklabels) & (checklabels == 0)))
        Cna = sum(as.numeric((pred.checkallres[idx,] != checklabels) & (checklabels == 1)))
        Cnn = sum(as.numeric((pred.checkallres[idx,] == checklabels) & (checklabels == 0)))
        Pra = (Caa+Cnn)/length(checklabels)
        Pre = (Caa+Can)*(Caa+Cna)/length(checklabels)/length(checklabels)+(Cnn+Can)*(Cnn+Cna)/length(checklabels)/length(checklabels)
        kappa = (Pra-Pre)/(1-Pre)
        kappathreshold = 0.57
        if(kappa > kappathreshold)
        {
#           testlabels <- factor(ifelse((testlabels==1),'control','treatment'))
          fit.predict<-predict(fit,data.frame(data=testdata[,features]),decision.values = TRUE)
          pred.allres[j,]=ifelse(attr(fit.predict,"decision.values") > threshold,1,0)
          j=j+1
        }

      }
      
#       count=1
#       for(range in ranges){
        pred.tmpres = pred.allres[1:j-1,]
        if(j>2){
          predict.votes = apply(pred.tmpres,2,sum)
        }
        else{
          predict.votes = pred.tmpres
        }
        
        fit.predict<-factor(ifelse(predict.votes > (j-1)/2,"control","treatment"))
        
        if(nlevels(fit.predict)>1){
          tmpres = confusionMatrix(fit.predict,testlabels)
          sp.tmpres[fold,1]=tmpres$byClass["Specificity"]
          sn.tmpres[fold,1]=tmpres$byClass["Sensitivity"]
          gmean.tmpres[fold,1] <- sqrt(sp.tmpres[fold,1]*sn.tmpres[fold,1])
          acc.tmpres[fold,1] <- tmpres$overall["Accuracy"]
        }
        else{
          if(fit.predict[1]=="control"){
            sn.tmpres[fold,]=1
            sp.tmpres[fold,]=0
          }
          else{
            sn.tmpres[fold,]=0
            sp.tmpres[fold,]=1    
          }
          gmean.tmpres[fold,] <- 0
          acc.tmpres[fold,] <- length(which(as.character(fit.predict)==as.character(testlabels)))/length(test)
        }
#         count=count+1
#       }
    }
    
    sn.res[times,]=apply(sn.tmpres,2,mean)
    sp.res[times,]=apply(sp.tmpres,2,mean)
    acc.res[times,]=apply(acc.tmpres,2,mean)
    gmean.res[times,]=apply(gmean.tmpres,2,mean)
  }
  means=apply(sn.res,2,mean)
  sn.res=rbind(sn.res,means)
  rownames(sn.res)=c(1:10,"mean")
  colnames(sn.res)=c(1)
  expfile<-paste(outputpath,"/",substr(datasets[k],1,3),"_5CV_sn(svm_thr).csv",sep="")
  write.csv(sn.res,expfile)
  
  means=apply(sp.res,2,mean)
  sp.res=rbind(sp.res,means)
  rownames(sp.res)=c(1:10,"mean")
  colnames(sp.res)=c(1)
  expfile<-paste(outputpath,"/",substr(datasets[k],1,3),"_5CV_sp(svm_thr).csv",sep="")
  write.csv(sp.res,expfile)
  
  means=apply(acc.res,2,mean)
  acc.res=rbind(acc.res,means)
  rownames(acc.res)=c(1:10,"mean")
  colnames(acc.res)=c(1)
  expfile<-paste(outputpath,"/",substr(datasets[k],1,3),"_5CV_acc(svm_thr).csv",sep="")
  write.csv(acc.res,expfile)
  
  means=apply(gmean.res,2,mean)
  gmean.res=rbind(gmean.res,means)
  rownames(gmean.res)=c(1:10,"mean")
  colnames(gmean.res)=c(1)
  expfile<-paste(outputpath,"/",substr(datasets[k],1,3),"_5CV_gmean(svm_thr).csv",sep="")
  write.csv(gmean.res,expfile)
  
}
