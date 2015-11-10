inputpath = "F:/My Data/Datasets"
outputpath = "F:/Classification/Datasets"
if(file.exists(outputpath)==FALSE)dir.create(outputpath)
memory.limit(4000)

library(RWeka,quietly=TRUE)

######### Leukemia(ALL-AML) ########
gc(reset=TRUE)
expfile<-paste(inputpath,"/ALL-AML_Leukemia/ALL-AML_train.arff",sep="") # When the train set is provided directly
trainset=read.arff(expfile)
expfile<-paste(inputpath,"/ALL-AML_Leukemia/ALL-AML_test.arff",sep="") # When the train set is provided directly
testset=read.arff(expfile)
dataset=rbind(trainset,testset)
rm(trainset,testset)

labels=dataset[,ncol(dataset)]
print(labels)
index=c(labels=="ALL")
data1=dataset[index,]
data1[,ncol(dataset)]=0

index=c(labels=="AML")
data2=dataset[index,]
data2[,ncol(dataset)]=1

rm(dataset)
finaldataset = rbind(data1,data2)

expfile<-paste(outputpath,"/Leukemia(ALL-AML)",sep="")
if(file.exists(expfile)==FALSE)dir.create(expfile)
expfile<-paste(expfile,"/dataset.csv",sep="")
write.csv(finaldataset,expfile, row.names = FALSE)

rm(labels,data1,data2,finaldataset,expfile,index)

write.arff(dataset,file="dataset(normalized).arff")
# data=data.frame(id=c(1:nrow(dataset)),dmm=c(1:nrow(dataset)),
#                 y=dataset[,ncol(dataset)],data=dataset[,-ncol(dataset)])
# expfile="dataset(norm).csv"
# write.csv(data,expfile, row.names = FALSE)

######### BreastCancer ########
gc(reset=TRUE)
expfile<-paste(inputpath,"/BreastCancer/breastCancer-train.arff",sep="") # When the train set is provided directly
trainset=read.arff(expfile)
expfile<-paste(inputpath,"/BreastCancer/breastCancer-test.arff",sep="") # When the train set is provided directly
testset=read.arff(expfile)
dataset=rbind(trainset,testset)
rm(trainset,testset)

labels=dataset[,ncol(dataset)]
print(labels)

index=c(labels=="non-relapse")
data1=dataset[index,]
data1[,ncol(dataset)]=0

index=c(labels=="relapse")
data2=dataset[index,]
data2[,ncol(dataset)]=1

rm(dataset)
finaldataset = rbind(data1,data2)

expfile<-paste(outputpath,"/BreastCancer",sep="")
if(file.exists(expfile)==FALSE)dir.create(expfile)
expfile<-paste(expfile,"/dataset.csv",sep="")
write.csv(finaldataset,expfile, row.names = FALSE)

rm(labels,data1,data2,finaldataset,expfile,index)

######### CNS ########
gc(reset=TRUE)
expfile<-paste(inputpath,"/CNS/centralNervousSystem-outcome.arff",sep="") # When the train set is provided directly
dataset=read.arff(expfile)

labels=dataset[,ncol(dataset)]
print(labels)

index=c(labels=="Class1")
data1=dataset[index,]
data1[,ncol(dataset)]=0

index=c(labels=="Class0")
data2=dataset[index,]
data2[,ncol(dataset)]=1

rm(dataset)
finaldataset = rbind(data1,data2)

expfile<-paste(outputpath,"/CNS",sep="")
if(file.exists(expfile)==FALSE)dir.create(expfile)
expfile<-paste(expfile,"/dataset.csv",sep="")
write.csv(finaldataset,expfile, row.names = FALSE)

rm(labels,data1,data2,finaldataset,expfile,index)


######### ColonTumor ########
gc(reset=TRUE)
expfile<-paste(inputpath,"/ColonTumor/colonTumor.arff",sep="") # When the train set is provided directly
dataset=read.arff(expfile)

labels=dataset[,ncol(dataset)]
print(labels)

index=c(labels=="positive")
data1=dataset[index,]
data1[,ncol(dataset)]=0

index=c(labels=="negative")
data2=dataset[index,]
data2[,ncol(dataset)]=1

rm(dataset)
finaldataset = rbind(data1,data2)

expfile<-paste(outputpath,"/ColonTumor",sep="")
if(file.exists(expfile)==FALSE)dir.create(expfile)
expfile<-paste(expfile,"/dataset.csv",sep="")
write.csv(finaldataset,expfile, row.names = FALSE)

rm(labels,data1,data2,finaldataset,expfile,index)

######### DLBCL-Stanford ########
gc(reset=TRUE)
expfile<-paste(inputpath,"/DLBCL-Stanford/DLBCL-Stanford.arff",sep="") # When the train set is provided directly
dataset=read.arff(expfile)

labels=dataset[,ncol(dataset)]
print(labels)

index=c(labels=="germinal")
data1=dataset[index,]
data1[,ncol(dataset)]=0

index=c(labels=="activated")
data2=dataset[index,]
data2[,ncol(dataset)]=1

rm(dataset)
finaldataset = rbind(data1,data2)

library(pcaMethods)
data=finaldataset[,-ncol(finaldataset)]
labels=finaldataset[,ncol(finaldataset)]
rm(data1,data2,finaldataset)
res=llsImpute(data,allVariables=TRUE)
filled.data <- completeObs(res)
finaldataset=data.frame(filled.data,labels)

expfile<-paste(outputpath,"/DLBCL(Stanford)",sep="")
if(file.exists(expfile)==FALSE)dir.create(expfile)
expfile<-paste(expfile,"/dataset.csv",sep="")
write.csv(finaldataset,expfile, row.names = FALSE)

rm(finaldataset,filled.data,data,labels,res)

######### LungCancer ########
gc(reset=TRUE)
expfile<-paste(inputpath,"/LungCancer-Harvard1/lung-harvard.arff",sep="") # When the train set is provided directly
dataset=read.arff(expfile)

labels=dataset[,ncol(dataset)]
print(labels)

index=c(labels=="NORMAL")
data1=dataset[index,]
data1[,ncol(dataset)]=0

index=c(labels!="NORMAL")
data2=dataset[index,]
data2[,ncol(dataset)]=1
rm(dataset)

finaldataset = rbind(data1,data2)
expfile<-paste(outputpath,"/LungCancer",sep="")
if(file.exists(expfile)==FALSE)dir.create(expfile)
expfile<-paste(expfile,"/dataset.csv",sep="")
write.csv(finaldataset,expfile, row.names = FALSE)

rm(labels,data1,data2,finaldataset,expfile,index)

######### ProststeCancer ########
gc(reset=TRUE)
expfile<-paste(inputpath,"/ProstateCancer/prostate_tumorVSNormal_train.arff",sep="") # When the train set is provided directly
trainset=read.arff(expfile)
expfile<-paste(inputpath,"/ProstateCancer/prostate_tumorVSNormal_test.arff",sep="") # When the train set is provided directly
testset=read.arff(expfile)
dataset=rbind(trainset,testset)
rm(trainset,testset)

labels=dataset[,ncol(dataset)]
print(labels)

index=c(labels=="Normal")
data1=dataset[index,]
data1[,ncol(dataset)]=0

index=c(labels=="Tumor")
data2=dataset[index,]
data2[,ncol(dataset)]=1

rm(dataset)
finaldataset = rbind(data1,data2)

expfile<-paste(outputpath,"/ProstateCancer",sep="")
if(file.exists(expfile)==FALSE)dir.create(expfile)
expfile<-paste(expfile,"/dataset.csv",sep="")
write.csv(finaldataset,expfile, row.names = FALSE)

rm(labels,data1,data2,finaldataset,expfile,index)

