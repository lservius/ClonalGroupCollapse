###### LOAD LIBRARIES ###### 
library(tidyverse)
library(plyr)
library(caret)
library(randomForest)
library(foreach)
library(doParallel)
print("libraries loaded")

###### FUNCTIONS ###### 
source("../functions/func_data-balance.R")
print("functions loaded")

###### DEFINE ENVIRONMENT ###### 
set.seed(346)
numCores <- 29 - 1

print("environment defined")

###### PROCESS DATA ###### 
covid <- read.delim("../data/processed/covid-cg.txt")

# remove pid
covid <- select(covid, -c(PatientID,CloneID))

print("data processed")

###### MODEL SET-UP ###### 
# make characters into factors
cat <- c("IsCrossClass", "Vgene", "Jfamily", "Dfamily", "Dgene", "Sex")
covid[cat] <- lapply(covid[cat], factor)

# create folds
n_k <- 5 
fold <- createFolds(covid$IsCrossClass,k=n_k)

###### RF ###### 
# parallelisation set-up
registerDoParallel(numCores)  

xval <- data.frame(data = "data", model = "RF", iteration = 1, data_balance = 1, sensitivity = 1, specificity = 1, uar = 1)
RF_hyPar <- data.frame(data = character(), model = character(), iteration = numeric(), m = numeric(), min = numeric(), N = numeric(), sens = numeric(), uar = numeric())

for(i in 1:n_k){
  # create list of training sets
  train <- do.call(rbind, lapply(fold[-i], function(y){
    tr <- covid[y,]
    
    return(tr)
  }))
  
  # create test set with 1 fold
  test <- do.call(rbind, lapply(fold[i], function(y){
    t <- covid[y,]
    return(t)
  }))
  
  # break train into train-test for tuning
  tune_fold <- createFolds(train$IsCrossClass, k = 10)
  
  hp_cv <- list()
  
  for(j in 1:10){
    
    tune_train <- do.call(rbind, lapply(tune_fold[-j], function(y){
      t4 <- train[y,]
      return(t4)
    }))
    tune_train <- balance.50(tune_train)
    
    tune_test <- do.call(rbind, lapply(tune_fold[j], function(y){
      t1 <- train[y,]
      return(t1)
    }))
    
    
    # create cluster
    cl <- makeCluster(numCores)
    
    # tuning hyperparamters
    hp <- expand.grid(m = 2^seq(0,6,1), 
                      min = 1,
                      N = seq(250,1000,250)
    )
    clusterExport(cl, c("tune_train", "tune_test", "hp"))
    
    clusterEvalQ(cl, {
      library(caret)
      library(randomForest)
    })
    result <- parApply(cl, hp, 1,function(x){
      model <- randomForest(IsCrossClass ~ .,
                            mtry = x[1], 
                            nodesize = x[2],
                            ntree = x[3],
                            data = tune_train)
      y <- predict(model, tune_test)
      cm <- confusionMatrix(data = y, reference = tune_test$IsCrossClass, positive = "TRUE")
      sens <- cm$byClass[[1]]
      uar <- 0.5*(cm$byClass[[1]] + cm$byClass[[2]])
      g.uar <- sqrt(cm$byClass[[1]] * cm$byClass[[2]]) * 100 # geometric uar
      
      f <- c(x[1],x[2],x[3],sens,g.uar)
      return(f)
    })
    
    stopCluster(cl)
    
    result <- data.frame(t(result))
    colnames(result) <- c("m","min","N","sens","g.uar")
    
    hp_cv[[j]] <- result
    
  }
  
  hp_cv <- do.call(rbind,hp_cv)
  
  # calculate mean sens and uar for CV
  hp_mean <- ddply(hp_cv, c("m","min","N"), summarise, sens.10fold = mean(sens), guar.10fold = mean(g.uar))
  
  hp_f <- hp_mean[which.max(hp_mean$guar.10fold),]
  
  train <- balance.50(train)
  
  db <- sum(train$IsCrossClass == TRUE)/nrow(train) * 100    # save out balance in response for train data
  
  # training the model
  rf <- randomForest(IsCrossClass ~ .,
                     mtry = hp_f$m, 
                     nodesize = hp_f$min,
                     ntree = hp_f$N, 
                     data = train)
  
  # getting expected response, y, using the trained model
  y <- predict(rf,test)
  
  # getting performance
  cm <- confusionMatrix(data = y, reference = test$IsCrossClass, positive = "TRUE")
  
  # committing performance metrics to data.frame
  sen <- cm$byClass[[1]] * 100                                 # sensitivity
  spec <- cm$byClass[[2]] * 100                                # specificity
  uar <- 0.5*(cm$byClass[[1]] + cm$byClass[[2]]) * 100         # unweighted accuracy recall
  
  
  # save out xval
  xval[i,] <- data.frame(data = "COVID", model= "RF", iteration = i, data_balance = db, sensitivity = sen, specificity = spec, uar = uar)
  RF_hyPar[i,] <- c("COVID","RF",i, hp_f)                          # save hyperparameter
}

print("RF done")

###### PRINT RESULTS ######
xval
RF_hyPar

write.csv(xval, file = "../results/CP_RF.csv")

print("results printed")

