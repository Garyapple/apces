# APBEFS: AP Based Ensemble feature selection methods

APBEFS <- function(traindata,classlabels,topK=300,nFit=45){
  require(apcluster,quietly=TRUE)
  require(WGCNA,quietly=TRUE)
  require(caret,quietly=TRUE)
  require(e1071)
  require(caTools)
  require(GeneSelector)
  
  # Aggregation 
  N=length(classlabels)
  control.num=sum(as.numeric(classlabels=="control"))
  treatment.num=N-control.num
  
  if(control.num<treatment.num){
    flag="treatment";
    basedata=traindata[which(classlabels=="control"),]
    minnum=control.num
  }
  else{
    flag="control";
    basedata=traindata[which(classlabels=="treatment"),]
    minnum=treatment.num
  }
  
  repeats=21 # rank times
  AggRanks=c()
  for(iter in 1:repeats){
    if(flag=="treatment"){
      index=which(classlabels=="treatment")
      data.smp=sample(1:treatment.num,minnum)
      classLab=c(rep(0,control.num),rep(1,control.num))
    }
    else{
      index=which(classlabels=="control")
      data.smp=sample(1:control.num,minnum)
      classLab=c(rep(1,treatment.num),rep(0,treatment.num))
    }
    data.add=data[index[data.smp],]
    data.tmp=rbind(basedata,data.add)
    rm(data.add,data.smp,index)
    
    # Aggregation 
    eBayes.res<-RankingLimma(t(data.tmp),classLab)
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
  control.samples = which(classlabels=="control")
  treatment.samples =  which(classlabels=="treatment")
  
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
  fit.ensembles=vector(mode="list",length=nFit) # models list
  fit.features=vector(mode="list",length=nFit) # models list

  threshold = (control.num - treatment.num)/(N+2) #a=1
               
  for(idx in 1:nFit){
    print(idx)
    features=feature.smp[idx,]
    modeldata = traindata[,features]
    fit.ensembles[[idx]] <- svm(classL~.,data.frame(data= modeldata,classL=classlabels)) # save the tunne model
    fit.features[[idx]]=c(features)
  }
  
  ensemble.res<- list(models=fit.ensembles,
                      models.features=fit.features,
                      models.thresholds=threshold,
                      models.num=nFit)
  return(ensemble.res)
}

predict.APBEFS<-function(EFSResult,newdata){
  fit.models = EFSResult$models
  feature.subsets = EFSResult$models.features
  fit.thresholds= EFSResult$models.thresholds
  nModel = EFSResult$models.num
  
  # Only one test object
  newdata=as.matrix(newdata)
  data.dim=dim(newdata)
  if(min(data.dim)==1){
    dim(newdata)=c(1,max(data.dim))
    newdata=rbind(newdata,newdata)
  }
  
  # Prediction
  pred.ensemble=c()
  for(ii in 1:nModel){
    features=feature.subsets[[ii]]
    modeldata=newdata[,features]
    current.fit=fit.models[[ii]]
    pred.tmp=predict(current.fit,data.frame(data=modeldata),decision.values =TRUE)
    current.class=ifelse(attr(pred.tmp,"decision.values") > fit.thresholds,"control","treatment")
    if(ii==1){
      pred.ensemble=current.class
    }
    else{
      pred.ensemble=data.frame(pred.ensemble,current.class)
    }
  }
  
  # Combine the predictions
  if(nModel>1){
    pred.class=apply(pred.ensemble,1,function(x){
      return(names(which.max(table(x))))
    })
  }
  else{
    pred.class <- pred.ensemble
  }
  
  # Only one object
  if(min(data.dim)==1){
    pred.class=pred.class[1]
  }
  return(pred.class)
}