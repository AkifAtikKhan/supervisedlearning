# (1). 
#

conductors.train = read.csv(file.choose(),header=T,sep=',')
conductors.test = read.csv(file.choose(),header=T,sep=',')

##Log transform The two energy variables. 

conductors.train$formation_energy_ev_natom = log(conductors.train$formation_energy_ev_natom + 1) 
conductors.train$bandgap_energy_ev = log(conductors.train$bandgap_energy_ev + 1) 

conductors.test$bandgap_energy_ev=  log(conductors.test$bandgap_energy_ev + 1) 
conductors.test$formation_energy_ev_natom = log(conductors.test$formation_energy_ev_natom + 1) 

              ###           Basic simple linear  regression                  ########
lm.formenergy = lm(formation_energy_ev_natom~ spacegroup + number_of_total_atoms + percent_atom_al + percent_atom_ga + percent_atom_in +
                 lattice_vector_3_ang+lattice_angle_alpha_degree+lattice_angle_beta_degree + lattice_angle_gamma_degree, data = conductors.train)
lm.bandenergy = lm(bandgap_energy_ev~ spacegroup + number_of_total_atoms + percent_atom_al + percent_atom_ga + percent_atom_in +
                     lattice_vector_3_ang+lattice_angle_alpha_degree+lattice_angle_beta_degree + lattice_angle_gamma_degree, data = conductors.train)

##Linear Regression model RMSEP. 
library(randomForest)

multmodel.cv = function(X,y,p=.667,B=100) 
{
  n = length(y)
  X = scale(X)
  data = data.frame(X,y)
  cv <- rep(0,B)
  for (i in 1:B) 
  {
    ss <- floor(n*p)
    sam <- sample(1:n,ss,replace=F)
    fit2 <- lm(y~., data=data[sam,])
    fit3 <- randomForest(formula = y ~ ., data = data[sam,], importance = T)
    ynew <- predict(fit2, newdata=data[-sam,])
    ynew2 <- predict(fit3, newdata=data[-sam,])
    ##Plotting the graph. 
    plot(y[-sam]~ynew)
    cv[i] <- sqrt(mean((y[-sam]-ynew)^2))
    cv[i] <- sqrt(mean((y[-sam]-ynew2)^2))
  }
  cv
}

randomForest(formation_energy_ev_natom~ spacegroup + number_of_total_atoms + percent_atom_al + percent_atom_ga + percent_atom_in +
lattice_vector_3_ang+lattice_angle_alpha_degree+lattice_angle_beta_degree + lattice_angle_gamma_degree, importance = T, data = conductors.train)








##All Predictors from our models  in this section #####      
lmpred.form = predict(lm.formenergy, conductors.test)
lmpred.band = predict(lm.bandenergy,conductors.test)

lmdf.predform = data.frame(ID=conductors.test$ID,ypred=lmpred.form)
lmdf.predband = data.frame(ID=conductors.test$ID,ypred=lmpred.band)

###### Exploratory Data viz'z.
cor.train = cor(conductors.train[,2:14])
corrplot::corrplot(cor.train)
pairs.plus(spacegroup~., data = conductors.train)

#### Create model. 
create_model <- function(trainData,target) 
{
  set.seed(120)
  myglm <- glm(target ~ ., data=trainData, family = "gaussian")
  return(myglm)
}

##predictions
modelglm.formation <- create_model(conductors.train,conductors.train$formation_energy_ev_natom)
summary(modelglm.formation)
plot(md)

score <- predict(modelglm.formation, newdata = conductors.train, type ="response")
library(AUC)
help(auc)
auc(conductors.train$formation_energy_ev_natom, score)
auc(conductors.train$formation_energy_ev_natom, score)

sum(is.na(score))
sum(is.na(modelglm.formation))
sum(is.na(conductors.train))
sum(is.na(RMSEP.train))

#######
score_train <- predict(myglm, newdata = complete, type = "response")
RMSEP.train = sqrt(conductors.train$formation_energy_ev_natom - score)

#######

RMSEP = rep(0,M)
MAEP = rep(0,M)
MAPEP = rep(0,M)
n = nrow(X)
for (i in 1:M) 
  {
  ss = floor(n*p)
  sam = sample(1:n,ss,replace=F)
  fit = glmnet(X[sam,],y[sam],lambda=lambda,alpha=alpha)
  ypred = predict(fit,newx=X[-sam,])
  RMSEP[i] = sqrt(mean((y[-sam]-ypred)^2))
  MAEP[i] = mean(abs(y[-sam]-ypred))
  yp = ypred[y[-sam]!=0]
  ya = y[-sam][y[-sam]!=0]
  MAPEP[i]=mean(abs(yp-ya)/ya)
}
###################### GLM MODELS ###########################
create_model <- function(trainData,target) 
 {
  set.seed(120)
  myglm <- glm(target ~ ., data=trainData, family = "gaussian")
  }
######### splitting data and Response and Predictors ########################
bandenergy = as.matrix(conductors.train[,14])
formenergy = as.matrix(conductors.train[,13])
predictors = as.matrix(conductors.train[,2:12])

View(formenergy)
View(conductors.train)


sum(is.na(bandenergy))

sum(is.na(predictors))
sum(is.na(formenergy))


############### PRINCIPAL COMPONENT ANALYSIS MODELS ####################

findcomps.form = plsr(bandenergy~predictors, data = conductors.train)

help(plsr)
#PLS
formpls.monte = pls.cv(predictors, ncomp=3, formenergy,  p = .80)
bandpls.monte = pls.cv(predictors, ncomp=3, bandenergy,  p = .80)

#PCR
formpcr.monte = pcr.cv(predictors, ncomp=3, formenergy,  p = .80)
bandpcr.monte = pcr.cv(predictors, ncomp=3, bandenergy,  p = .80)

## change spacegroup to a factor. 
conductors.train$spacegroup = as.factor(conductors.train$spacegroup)

mean(formpcr.monte)
mean(formpls.monte)
max(formpcr.monte)
min(formpcr.monte)


################                  Random Forest             #############
set.seed(4)

library(randomForest)
install.packages("xgboost")

### Split the data 

sam = sample(1:1567,floor(1567*.6667),replace=FALSE)
conductors.trainsubset = conductors.train[sam,]
conductors.validsubset = conductors.train[-sam,]

nrow(Conductors.trainsubset)
nrow(Conductors.validsubset)

View(Conductors.trainsubset)





###Monte carlo multiple model function, 
models.cv = function(X,y,p=.667,B=100) 
{
  n = length(y)
  X = scale(X)
  data = data.frame(X,y)
  cv <- rep(0,B)
  for (i in 1:B) 
  {
    ss <- floor(n*p)
    sam <- sample(1:n,ss,replace=F)
    fit2 <- randomForest(ly ~ ., data = data[sam], ntree = 100, nodesize = 50)
    dtrain <- xgb.DMatrix(data = train$data, label = train$label)
   # fit3 <- xgboost(data = dtrain, max_depth = 2, eta = 1, nthread = 2, nrounds = 2, objective = "regression")
    ynew <- predict(fit2, newdata=data[-sam,])
   # ynew2 <- predict(fit3, newdata=data[-sam,])
    ##Plotting the graph. 
    plot(y[-sam]~ynew)
  #  plot(y[-sam]~ynew2)
    cv[i] <- sqrt(mean((y[-sam]-ynew)^2))
 #  cv[i] <- sqrt(mean((y[-sam]-ynew2)^2))
  }
  cv
}

## Random forest models for form energy using 100 n-trees. 
rf1.formenergy = randomForest(formation_energy_ev_natom~ spacegroup + number_of_total_atoms + percent_atom_al + percent_atom_ga + percent_atom_in +
                     lattice_vector_3_ang+lattice_angle_alpha_degree+lattice_angle_beta_degree + lattice_angle_gamma_degree,ntree=50, data = conductors.trainsubset)


rf2.formenergy = randomForest(formation_energy_ev_natom~ spacegroup + number_of_total_atoms + percent_atom_al + percent_atom_ga + percent_atom_in +
                               lattice_vector_3_ang+lattice_angle_alpha_degree+lattice_angle_beta_degree + lattice_angle_gamma_degree,ntree=75, data = conductors.trainsubset)

rf3.formenergy = randomForest(formation_energy_ev_natom~ spacegroup + number_of_total_atoms + percent_atom_al + percent_atom_ga + percent_atom_in +
                                lattice_vector_3_ang+lattice_angle_alpha_degree+lattice_angle_beta_degree + lattice_angle_gamma_degree,ntree=90, data = conductors.trainsubset)



formenergy1.sscv = rf.sscv(rf1.formenergy,conductors.trainsubset[,2:13] , mtry = 3)
formenergy2.sscv = rf.sscv(rf2.formenergy,conductors.trainsubset[,2:13] , mtry = 3)
formenergy3.sscv = rf.sscv(rf3.formenergy,conductors.trainsubset[,2:13] , mtry = 3)

rf1.predval = predict(rf1.formenergy, newdata = conductors.validsubset)
rf2.predval = predict(rf2.formenergy, newdata = conductors.validsubset)
rf3.predval = predict(rf3.formenergy, newdata = conductors.validsubset)

rf1.predtest = predict(rf1.formenergy, newdata = conductors.test)
rf2.predtest = predict(rf2.formenergy, newdata = conductors.test)
rf3.predtest = predict(rf3.formenergy, newdata = conductors.test)

View(rf1.predtest)
View(data.frame(rf2.predtest))

rf.testpredall = data.frame(rf1.predtest,rf2.predtest,rf3.predtest)

View(rf.testpredall)



as.(rf1.predval)
randomforestmodels = data.frame(rf1.predval,rf2.predval,rf3.predval, conductors.validsubset$formation_energy_ev_natom)

rf.predictionscsv = write.csv(rf.testpredall, file = "predictions.csv")

training = write.csv(conductors.trainsubset, file = "conductors.trainsubset.csv")
validation = write.csv(conductors.validsubset ,  file = "conductors.validsubset.csv")

kapilspredictions = read.csv(file.choose(), header = TRUE)

View(rf.testpredall)

allpredictions1 = as.data.frame(kapilspredictions)



View(rf.testpredall)

formenergy = conductors.trainsubset$formation_energy_ev_natom


### Randomforest models for band energy using 100 n-trees

rf1.bandenergy = randomForest(bandgap_energy_ev~ spacegroup + number_of_total_atoms + percent_atom_al + percent_atom_ga + percent_atom_in +
                               lattice_vector_3_ang+lattice_angle_alpha_degree+lattice_angle_beta_degree + lattice_angle_gamma_degree,ntree=100, data = conductors.trainsubset)
rf2.bandenergy = randomForest(bandgap_energy_ev~ spacegroup + number_of_total_atoms + percent_atom_al + percent_atom_ga + percent_atom_in +
                                lattice_vector_3_ang+lattice_angle_alpha_degree+lattice_angle_beta_degree + lattice_angle_gamma_degree,ntree=200, data = conductors.trainsubset)
rf3.bandenergy = randomForest(bandgap_energy_ev~ spacegroup + number_of_total_atoms + percent_atom_al + percent_atom_ga + percent_atom_in + lattice_vector_3_ang+lattice_angle_alpha_degree+lattice_angle_beta_degree + lattice_angle_gamma_degree, ntree=300 , data = conductors.trainsubset)


plot(rf1.bandenergy)

plot(rf2.bandenergy)

plot(rf3.bandenergy)


rf1.predtest = predict(rf1.bandenergy, newdata = conductors.test)
rf2.predtest = predict(rf2.bandenergy, newdata = conductors.test)
rf3.predtest = predict(rf3.bandenergy, newdata = conductors.test)




allbandpredictions = data.frame(rf1.predtest, rf2.predtest, rf3.predtest)

View(allbandpredictions)
csvbandpredictions = write.csv(allbandpredictions, "allbandpredictions.csv") 


formenergy1.sscv = rf.sscv(rf2.bandenergy,conducband.trainsubset , mtry = 3)
formenergy.sscv = rf.sscv(rf.bandenergy,conducband.trainsubset[,2:13] , mtry = 3)




#change the name
pred = predict(rf.formenergy,newdata = conductors.test)
write.csv(pred)
formenergy.sscv

s(conductors.train)














###############
lmpred.form = predict(lm.formenergy, conductors.test)
lmpred.band = predict(lm.bandenergy,conductors.test)

lmdf.predform = data.frame(ID=conductors.test$ID,ypred=lmpred.form)
lmdf.predband = data.frame(ID=conductors.test$ID,ypred=lmpred.band)

predictors = as.matrix(predictors)
formenergy = as.matrix(formenergy)
View(formenergy)
View(conductors.train)

formenergy.rf = randomForest(formenergy~predictors,data=conductors.train, ntree=100)

View(doc)

bandenergy.sscv = rf.sscv(formenergy.rf,Concrete.train)

