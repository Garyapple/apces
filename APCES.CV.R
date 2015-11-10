# APBEFS.CV: cross validation of AP based ensemble feature selection methods

APBEFS.CV <- function(data,classlabels,nFit=45,topK=300,folds=5){
  library(caret)
  source("APCES.R")
  
  N=length(classlabels)
  
  res=c()
  sn.res = array(0,folds)
  sp.res = array(0,folds)
  acc.res = array(0,folds)
  gmean.res = array(0,folds)
  
  cv.index=createFolds(classlabels,k=folds)
  for (i in 1:folds) {
    print(paste("Fold:",i))
    test <- cv.index[[i]]
    
    fit <- APBEFS(traindata=data[-test,],classlabels=classlabels[-test],
                  topK=topK,nFit=nFit)
    fit.predict<-predict.APBEFS(fit, data[test,])
    
    if(nlevels(factor(fit.predict))>1){
      tmpres = confusionMatrix(fit.predict,classlabels[test])
      sp.res[i]=tmpres$byClass["Specificity"]
      sn.res[i]=tmpres$byClass["Sensitivity"]
      gmean.res[i] <- sqrt(sp.res[i]*sn.res[i])
      acc.res[i] <- tmpres$overall["Accuracy"]
    }
    else{
      if(fit.predict[1]=="control"){
        sn.res[i]=1
        sp.res[i]=0
      }
      else{
        sn.res[i]=0
        sp.res[i]=1    
      }
      gmean.res[i] <- 0
      acc.res[i] <- length(which(fit.predict==classlabels[test]))/length(test)
    }
  }
  res$sn=mean(sn.res)
  res$sp=mean(sp.res)
  res$acc=mean(acc.res)
  res$gmean=mean(gmean.res)
  return(res)
}