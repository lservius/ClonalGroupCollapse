###### LOAD LIBRARIES ###### 
library(tidyverse)
library(plyr)
library(caret)
library(e1071)
library(foreach)
library(doParallel)
print("libraries loaded")
###### FUNCTIONS ###### 
source("../functions/func_data-balance.R")
print("functions loaded")

###### DEFINE ENVIRONMENT ###### 
set.seed(346)
numCores <- 33 - 1

print("environment defined")

###### PROCESS DATA ###### 
rsv <- read.delim("../data/processed/rsv-cg.txt")

# remove pid
rsv <- select(rsv, -c(CloneID,Sex))

print("data processed")

###### MODEL SET-UP ###### 
# make characters into factors
cat <- c("IsCrossClass", "Vgene", "Jfamily", "Dfamily", "Dgene", "Age")
rsv[cat] <- lapply(rsv[cat], factor)

# patient separation
pid <- unique(rsv$PatientID)
rsv <- lapply(pid, function(x){
  p <- subset(rsv, PatientID == x , select = -c(PatientID))
})

# patient combinations is now (d-1)vs1 
n_s <- length(pid)  

###### SVM ######
xval <- data.frame(data = character(), model = character(), iteration = numeric(), data_balance = numeric(), sensitivity = numeric(), specificity = numeric(), uar = numeric())

SVM_hyPar <- data.frame(data = character(), model = character(), iteration = numeric(), g = numeric(), c = numeric(), sens = numeric(), uar = numeric())

# grid for hyperparameter search
hp <- expand.grid(g = 2^seq(-15,1,3), 
                  c = 2^seq(-4,10,2)
)
for (i in 1:n_s){
  # create list of training sets
  train <- do.call(rbind, lapply(rsv[-i], function(y){
    tr <- y
    return(tr)
  }))
  
  # create test set with remaining patient
  test <- do.call("rbind", rsv[i])
  
  # break train into train-test for tuning
  tune_d <- rsv[-i]
  nfolds <- length(tune_d)
  
  hp_cv <- list()
    
  # tuning hyperparamters
  
  for(j in 1:nfolds){
    
    tune_train <- do.call(rbind, lapply(tune_d[-j], function(y){
      t4 <- y
      return(t4)
    }))
    tune_train <- balance.50(tune_train)
    
    tune_test <- do.call(rbind, tune_d[j])
    
    
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
      return(f)
    })
    
    stopCluster(cl)
    
    result <- data.frame(t(result))
    colnames(result) <- c("g","c","sens","g.uar")
    
    hp_cv[[j]] <- result
    
  }
  
  hp_cv <- do.call(rbind,hp_cv)
  
  # calculate mean sens and uar for CV
  hp_mean <- ddply(hp_cv, c("g","c"), summarise, sens.nfold = mean(sens), guar.nfold = mean(g.uar))
  
  # select highest mean uar
  hp_f <- hp_mean[which.max(hp_mean$guar.nfold),]
  
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
  
  # commiting performance metrics to data.frame
  sen <- cm$byClass[[1]] * 100                         # sensitivity
  spec <- cm$byClass[[2]] * 100                        # specificity
  uar <- 0.5*(cm$byClass[1] + cm$byClass[[2]]) * 100   # unweighted accuracy recall
  
  xval[i,] <- c("RSV","SVM",i,db,sen,spec,uar)           # save performance metrics
  SVM_hyPar[i,] <- c("RSV","SVM",i, hp_f)                # save hyperparameters
}

# checkpoint
  write.csv(xval[i,] , paste0("../results/", Sys.Date(),".",i,".RSV_DvsD_SVM.csv"), row.names = FALSE)
  write.csv(SVM_hyPar[i,] , paste0("../results/", Sys.Date(),".",i,".RSV_DvsD_SVMhyPar.csv"), row.names = FALSE)

print("SVM done")

###### PRINT RESULTS ######
xval
SVM_hyPar

print("results printed")
