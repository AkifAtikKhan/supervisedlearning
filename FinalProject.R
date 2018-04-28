##Packages to run all the functions. 

library(rando)

## Misclass Function 

misclass = function(fit,y) 
{
  temp <- table(fit,y)
  cat("Table of Misclassification\n")
  cat("(row = predicted, col = actual)\n")
  print(temp)
  cat("\n\n")
  numcor <- sum(diag(temp))
  numinc <- length(y) - numcor
  mcr <- numinc/length(y)
  cat(paste("Misclassification Rate = ",format(mcr,digits=3)))
  cat("\n")
}
####

## Final Project
WaterTrain = read.csv(file="http://course1.winona.edu/bdeppa/Stat%20425/Data/WaterSol%20(Train).csv")
WaterTest = read.csv(file="http://course1.winona.edu/bdeppa/Stat%20425/Data/WaterSol%20(Test).csv")
WaterTrain$y = as.factor(WaterTrain$y)

#### Principal Components. 
WaterTrain.pc = read.csv(file="http://course1.winona.edu/bdeppa/Stat%20425/Data/WaterTrainPC.csv")
WaterTest.pc = read.csv(file="http://course1.winona.edu/bdeppa/Stat%20425/Data/WaterTestPC.csv")
WaterTrain.pc$y = as.factor(WaterTrain.pc$y)
water.trainpc = WaterTrain.pc[sam,]
water.validpc = WaterTrain.pc[-sam,]

####
Train.X = WaterTrain[,1:71]
Test.X = WaterTest[,-1]  # Remove idnum column, the rest are predictors
Water.X = rbind(Train.X,Test.X)
Water.X = scale(Water.X)
trainX = Water.X[1:4000,]
testX = Water.X[-c(1:4000),]
Wtrain = data.frame(Solubility=as.factor(WaterTrain$y),trainX)
Wtest = data.frame(testX)
water.train = Wtrain[sam,]
water.valid = Wtrain[-sam,]
####

library(xgboost)
library(Matrix)

set.seed(1)
sam = sample(1:4000,floor(4000*.6666),replace=F)

water.train = Wtrain[sam,]
water.valid = Wtrain[-sam,]

##### XGBOOST and MATRIX!!#### DEPPA'S NOTES. 
Train.X = Matrix::sparse.model.matrix(Solubility~.-1,data=water.train)
Valid.X = Matrix::sparse.model.matrix(Solubility~.-1,data=water.valid)
Test.X = Matrix::sparse.model.matrix(sample(c(0,1),1631,replace=TRUE)~.-1,data=Wtest)
Train.Y = as.list(as.numeric(water.train$Solubility)-1)
Valid.Y = as.list(as.numeric(water.valid$Solubility)-1)
dtrain = xgboost::xgb.DMatrix(data=Train.X,label=Train.Y)
dvalid = xgboost::xgb.DMatrix(data=Valid.X,label=Valid.Y)
watchlist1 = list(train=dtrain,test=dvalid)
### 
mod1 = xgboost::xgb.train(data=dtrain,max.depth=4,eta=.1,nthread=6,nround=300,watchlist=watchlist1,objective="binary:logistic",eval_metric="auc",colsample_bytree=.5)
yhat.gb = predict(mod1,type="prob", newdata = dvalid)

## * need to calculate the '1's from a vector matrix for ROC Curve.* 

### ROC CURVES FOR XGBOOST MODEL ####
pred. = prediction(yhat1,Watersol.valid$y)
perf = performance(pred,'tpr',"fpr")
plot(perf)  draws the ROC curve
perf.auc = performance(pred,”auc”)

##########
## Measuring Model Performance if greater than 0.5 THEN YES else NO. for  both 0's and 1'S.  
err <- mean(as.numeric(yhat.xg > 0.5) != water.valid$Solubility)
print(paste("classification error=", err))
###

####  $$$ My codes for correlation plots + looking for  best PC'S. $$$$ #####
Acor = cor(water.train[,-1])
corrplot::corrplot(Acor,order="hclust",tl.cex=.01)
water.pca = prcomp(Wtrain[,-1],scale=TRUE)
summary(water.pca)
#12 Principal Components seems best. 
water12 = water.pca$x[,1:12]
W.trainpc = data.frame(Solubility=water.train$Solubility,water12[sam])
W.testpc = data.frame(W.trainpc[-sam])
#######



library(xgboost)
library(Matrix)
install.packages("ROCR")
library(ROCR)


##### SVN ####
require(kernlab)
## Loading required package: kernlab
require(caret)
library(kernlab)
library(caret) 
require(pROC)
set.seed(333)

### Models for Support Vector Machines. 
ctrl = trainControl(method="repeatedcv",repeats=5,summaryFunction=twoClassSummary,classProbs=TRUE)
svm.tune = train(x = water.train[, -1], y=make.names(water.train$Solubility),method="svmRadial",tuneLength=9,preProc=c("center","scale"),metric="ROC",trControl=ctrl)
svm.tune

ksvm.tuned = ksvm(Solubility~.,data=water.train,kpar=list(sigma=.005),C=2.5,cross=5)
ksvm.tuned
yhat.ksvm = predict(ksvm.tuned,newdata = water.train)
misclass(yhat.ksvm,water.train$Solubility)

#######

##Multi cores Trial. 
#This is a package from the applied predictive modeling book:
#here's the link about what it does: https://www.r-bloggers.com/simple-parallel-randomforest-with-domc-package/

##Random Forest using DoMC. You could Multicores for ANY MODEL Techniques. 
install.packages("doMC")
library(doMC)
library(foreach)

## Different PC Combinations to see how many is optimal. 
water.pc15 = water.trainpc[,1:16]
water.pc15val = water.validpc[,1:16]

water.pc10 = water.trainpc[,1:11]
water.pc10val = water.validpc[,1:11]

water.pc12 = water.trainpc[,1:12]
water.pc12val = water.validpc[,1:12]

registerDoMC(4) #number of cores on the machine
darkAndScaryForest <- foreach(y= seq(50), .combine=combine) %dopar%{
  rf = randomForest::randomForest(Solubility~., data = water.train, ntree=50, 
  norm.votes=FALSE) }


############ Based on the Trial and Error it seems that not using PC's gives a more optimal 
##****classification error****  but it could be overfitting!! 



parallel.ypred = predict(darkAndScaryForest, newdata = water.valid, type = "prob")

water.pc12$y = as.factor(water.pc12$y)
water.pc12val$y = as.factor(water.pc12val$y)


## Measuring Model Performance if greater than 0.5 THEN YES else NO. for  both 0's and 1'S.  
err <- mean(as.numeric(parallel.ypred > 0.5) != water.valid$Solubility)
print(paste("classification error=", err))

###

