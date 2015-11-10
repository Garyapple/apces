
EGSG.CV <- function(data,disc.data,classlabels,nFit=31,folds=5){
  library(caret)
  source("EGSG(new).R")
  
  N=length(classlabels)
  
  sn.tmpres=array(0,folds)
  sp.tmpres=array(0,folds)
  acc.tmpres=array(0,folds)
  gmean.tmpres=array(0,folds)
  
  cv.index=createFolds(classlabels,k=folds)
  
  for (fold in 1:folds) {
    print(paste("Fold:",fold))
    test <- cv.index[[fold]]
    fit <- EGSG(norm.data=data[-test,],disc.data=disc.data[-test,],
                classlabels=classlabels[-test],nFit=nFit)
    
    N=length(classlabels[-test])
    control.num=sum(as.numeric(classlabels[-test]=="control"))
    treatment.num=N-control.num
    
    threshold = (control.num - treatment.num)/(N+2) #a=1
    fit.predict<-predict.EGSG(fit, data[test,],threshold)
    
    if(nlevels(fit.predict)>1){
      tmpres = confusionMatrix(fit.predict,classlabels[test])
      sp.tmpres[fold]=tmpres$byClass["Specificity"]
      sn.tmpres[fold]=tmpres$byClass["Sensitivity"]
      gmean.tmpres[fold] <- sqrt(sp.tmpres[fold]*sn.tmpres[fold])
      acc.tmpres[fold] <- tmpres$overall["Accuracy"]
    }
    else{
      if(fit.predict[1]=="control"){
        sn.tmpres[fold]=1
        sp.tmpres[fold]=0
      }
      else{
        sn.tmpres[fold]=0
        sp.tmpres[fold]=1    
      }
      gmean.tmpres[fold] <- 0
      acc.tmpres[fold] <- length(which(as.character(fit.predict)==as.character(classlabels[test])))/length(test)
    }
  }
  res=c()
  res$sn=mean(sn.tmpres)
  res$sp=mean(sp.tmpres)
  res$acc=mean(acc.tmpres)
  res$gmean=mean(gmean.tmpres)
  
  return(res)
}