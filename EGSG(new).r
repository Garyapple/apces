EGSG <- function(norm.data,disc.data,classlabels,nFit=31){
  require(caret)
  require(e1071)
  require(RWeka)
  source("myInfo.r")
  
  print("----------------------->")
  gc(FALSE,TRUE)
  info.cc=array(1,ncol(disc.data))
  attr.entropies = apply(disc.data,2,my.entropy)
  class.entropy = my.entropy(classlabels)
  ac.joint.entropies=apply(disc.data,2,function(t){
    return(my.joint.entropy(t,classlabels))
  })
  ac.mi=class.entropy+attr.entropies-ac.joint.entropies
  info.cc=ac.mi/ac.joint.entropies
  
  deta=quantile(info.cc,0.3)
  org.index=which(info.cc>=deta)
  norm.data=norm.data[,org.index]
  disc.data=disc.data[,org.index]
  attr.entropies=attr.entropies[org.index]
  info.cc=info.cc[org.index]
  Top_k=1500
  
  # sort
  if(ncol(disc.data)<=Top_k){
    sign.gene=sort(info.cc,decreasing=TRUE,index.return=TRUE)$ix
  }
  else{
    sign.gene=sort(info.cc,decreasing=TRUE,index.return=TRUE)$ix[1:Top_k]  
  }
  
  org.index=org.index[sign.gene]
  norm.data=norm.data[,sign.gene]
  disc.data=disc.data[,sign.gene]
  attr.entropies=attr.entropies[sign.gene]
  info.cc=info.cc[sign.gene]
  
  rm(ac.joint.entropies,ac.mi,class.entropy,sign.gene)
  
  gc(FALSE,TRUE)
  GS=list()
  GS.center=c()
  Grem=c(1:ncol(disc.data))
  Grel=c()
  i=1
  repeat{
    print(i)
    center=Grem[1]
    Grel=Grem[-1]
    GS.center=c(GS.center,center)
    
    t=disc.data[,center]
    mat=disc.data[,Grel]
    aa.joint.entropies=apply(data.matrix(mat),2,function(x){
      return(my.joint.entropy(x,t))
    })
    
    e.xy=attr.entropies[center]+attr.entropies[Grel]
    aa.mi=e.xy-aa.joint.entropies
    
    info.fc=matrix(1,1,length(Grel))
    info.fc=aa.mi/aa.joint.entropies
    
    index=c(info.fc>=info.cc[Grel])
    GS.i=Grel[index]
    GS[[i]]=c(center,GS.i)
    i=i+1
    
    index=c(info.fc<info.cc[Grel])
    Grem=Grel[index]
    if(length(Grem)==1){
      GS[[i]]=Grem
      GS.center=c(GS.center,Grem)
      break
    }
    if(length(Grem)<1){
      break
    }
  }
  rm(index,GS.i,Grem,info.cc,info.fc,attr.entropies,center,i,disc.data)
  
  EP=length(GS)
  smp<-matrix(0,nFit,EP)
  smp[1,]=GS.center # use  the best gene subsets
  t=15 # the maximize number of genes in each cluster could be selected 
  
  if(nFit>1){
    for(j in 1:EP){
      cluster=GS[[j]]
      len=length(cluster)
      if(len>=t){
        tmp=sample(1:t,(nFit-1),replace=TRUE) # random select feature sets
      }
      else{
        tmp=sample(1:len,(nFit-1),replace=TRUE) # random select feature sets
      }
      smp[2:nFit,j]=cluster[tmp]
    }
    rm(tmp,len,cluster)
  }
  
  ############### Classification seting #################
  gc(FALSE,TRUE) 
  fit.ensemble=vector(mode="list",length=nFit) # models list
  fit.features=vector(mode="list",length=nFit) # models list
  
  for(i in 1:21){
    features=smp[i,]
    modeldata=norm.data[,features]
    #     fit.ensemble[[i]] <- IBk(classL~.,data.frame(data=modeldata,classL=classlabels),
    #                              control=Weka_control(K=3)) # save the tunne model
    fit.ensemble[[i]] <- svm(classL~.,data.frame(data=modeldata,classL=classlabels)) # save the tunne model
    fit.features[[i]]<-org.index[features]
  }    
  # the final result of fit ensemble 
  ensemble.res<- list(models=fit.ensemble,
                      models.features=fit.features,
                      models.num=nFit)
  return(ensemble.res)
}  

predict.EGSG<-function(EFSResult,newdata,threshold){
  fit.models = EFSResult$models
  feature.subsets = EFSResult$models.features
  nModel = EFSResult$models.num
  
  newdata=as.matrix(newdata)
  data.dim=dim(newdata)
  if(min(data.dim)==1){
    dim(newdata)=c(1,max(data.dim))
    newdata=rbind(newdata,newdata)
  }
  
  pred.allres = matrix(0,nModel,nrow(newdata))
  for(ii in 1:nModel){
    features=feature.subsets[[ii]]
    modeldata=newdata[,features]
    current.fit=fit.models[[ii]]
    current.pred=predict(current.fit,data.frame(data=modeldata),decision.values = TRUE) # SVM_THR
    pred.allres[ii,]=ifelse(attr(current.pred,"decision.values") > threshold,1,0)
  }
  predict.votes = apply(pred.allres,2,sum)
  pred.class<-factor(ifelse(predict.votes > nModel/2,"control","treatment"))
  
  if(min(data.dim)==1){
    pred.class=pred.class[1]
  }
  return(pred.class)
}