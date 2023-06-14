###### LOAD LIBRARIES ###### 
library(tidyverse)
library(plyr)
library(caret)
library(e1071)
library(foreach)
library(doParallel)
print("libraries loaded")

###### ITERATION SET ######
args <- commandArgs(trailingOnly = TRUE)
print(args)

args <- as.numeric(args)
print(args)

which.i<- args[1]

###### FUNCTIONS ###### 

# artificially re-balance data for classifier
balance.50 <- function(df, response = "IsCrossClass", minority = TRUE){
  
  df.min <- df[df[[response]] == minority,]
  df.maj <- df[df[[response]] != minority,]
  
  # get 50% data balance
  n.1 <- nrow(df.min)
  n.0 <- nrow(df.maj)
  
  delta <- n.0-n.1
  
  df.add <- df.min[sample(n.1, size = delta, replace = TRUE),]
  
  df.bal <- rbind(df,df.add)
  
  return(df.bal)
}

print("functions loaded")

###### DEFINE ENVIRONMENT ###### 
set.seed(346)
numCores <- 37 - 1
setwd("/scratch/users/k21111073/")

print("environment defined")

###### PROCESS DATA ###### 
rsv <- read.delim("data/rsv-cg.txt")

# remove pid
rsv <- select(rsv, -c(PatientID,CloneID,Sex))

print("data processed")

###### MODEL SET-UP ###### 
# make characters into factors
cat <- c("IsCrossClass", "Vgene", "Jfamily", "Dfamily", "Dgene", "Age")
rsv[cat] <- lapply(rsv[cat], factor)

# create folds
n_k <- 5 
fold <- createFolds(rsv$IsCrossClass,k=n_k)

###### SVM ######
xval <- data.frame(data = character(), model = character(), iteration = numeric(), data_balance = numeric(), sensitivity = numeric(), specificity = numeric(), uar = numeric())

SVM_hyPar <- data.frame(data = character(), model = character(), iteration = numeric(), g = numeric(), c = numeric(), sens = numeric(), uar = numeric())

# grid for hyperparameter search
hp <- expand.grid(g = 2^seq(-15,1,3), 
                  c = 2^seq(0,10,2)
)

for (i in which.i){
  # create list of training sets
  train <- do.call(rbind, lapply(fold[-i], function(y){
    tr <- rsv[y,]
    
    return(tr)
  }))
  
  # create test set with 1 fold
  test <- do.call(rbind, lapply(fold[i], function(y){
    t <- rsv[y,]
    return(t)
  }))
  
  # break train into train-test for tuning
  tune_fold <- createFolds(train$IsCrossClass, k = 10)
  
  hp_cv <- list()
    
  # tuning hyperparamters  
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
    
    clusterExport(cl, c("tune_train", "tune_test", "hp"))
    
    clusterEvalQ(cl, {
      library(caret)
      library(e1071)
    })
    
    result <- parApply(cl, hp, 1,function(x){
      model <- svm(formula = IsCrossClass ~ .,
                   data = tune_train,
                   type = 'C-classification',
                   kernel = 'radial',
                   gamma = x[1],
                   cost = x[2])
      
      y <- predict(model, tune_test)
      cm <- confusionMatrix(data = y, reference = tune_test$IsCrossClass, positive = "TRUE")
      sens <- cm$byClass[[1]]
      uar <- 0.5*(cm$byClass[[1]] + cm$byClass[[2]])
      g.uar <- sqrt(cm$byClass[[1]] * cm$byClass[[2]]) * 100 # geometric uar
      
      f <- c(x[1],x[2],sens,g.uar)
      
      rm(model)
      return(f)
    })
    
    stopCluster(cl)
    
    result <- data.frame(t(result))
    colnames(result) <- c("g","c","sens","g.uar")
    
    hp_cv[[j]] <- result
    
  }
  
  hp_cv <- do.call(rbind,hp_cv)
  
  # calculate mean sens and uar for CV
  hp_mean <- ddply(hp_cv, c("g","c"), summarise, sens.10fold = mean(sens), guar.10fold = mean(g.uar))
  
  # select highest mean uar
  hp_f <- hp_mean[which.max(hp_mean$guar.10fold),]
  
  # balance train
  train <- balance.50(train)
  db <- sum(train$IsCrossClass == TRUE)/nrow(train) * 100    # save out balance in response for train data
  
  # training the model
  svm <- svm(formula = IsCrossClass ~ .,
             data = train,
             type = 'C-classification',
             kernel = 'radial',
             gamma = hp_f$g,
             cost = hp_f$c)
  
  # getting expected response, y, using the trained model
  y <- predict(svm,test)
  
  # getting performance
  cm <- confusionMatrix(data = y, reference = test$IsCrossClass, positive = "TRUE")
  
  # committing performance metrics to data.frame
  sen <- cm$byClass[[1]] * 100                         # sensitivity
  spec <- cm$byClass[[2]] * 100                        # specificity
  uar <- 0.5*(cm$byClass[1] + cm$byClass[[2]]) * 100   # unweighted accuracy recall
  
  xval[i,] <- c("RSV","SVM",i,db,sen,spec,uar)           # save performance metrics
  SVM_hyPar[i,] <- c("RSV","SVM",i, hp_f)                # save hyperparameters
  
  # checkpoint
  write.csv(xval[i,] , paste0("output/", Sys.Date(),".",i,".RSV_SVM(xNIC).csv"), row.names = FALSE)
  write.csv(SVM_hyPar[i,] , paste0("output/", Sys.Date(),".",i,".RSV_SVMhyPar(xNIC).csv"), row.names = FALSE)
  
  rm(svm)
}

print("SVM done")

###### PRINT RESULTS ######
xval
SVM_hyPar

print("results printed")

