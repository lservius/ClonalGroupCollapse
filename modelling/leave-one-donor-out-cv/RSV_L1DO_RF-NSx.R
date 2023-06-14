###### LOAD LIBRARIES ###### 
library(tidyverse)
library(plyr)
library(caret)
library(randomForest)
library(foreach)
library(doParallel)
print("libraries loaded")

###### FUNCTIONS ###### 
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
numCores <- 124 - 1
setwd("/scratch/users/k21111073/")

print("environment defined")

###### PROCESS DATA ###### 
rsv <- read.delim("data/rsv-cg.txt")

# remove pid
rsv <- select(rsv, -c(CloneID))

print("data processed")

###### MODEL SET-UP ###### 
# make characters into factors
cat <- c("IsCrossClass", "Vgene", "Jfamily", "Dfamily", "Dgene", "Age", "Sex")
rsv[cat] <- lapply(rsv[cat], factor)

# patient separation
pid <- unique(rsv$PatientID)
rsv <- lapply(pid, function(x){
  p <- subset(rsv, PatientID == x , select = -c(PatientID))
})

# patient combinations is now (d-1)vs1 
n_s <- length(pid)  

###### RF ###### 
# parallelisation set-up
registerDoParallel(numCores)  

xval <- data.frame(data = "data", model = "RF-NSx", iteration = 1, data_balance = 1, sensitivity = 1, specificity = 1, uar = 1)
RF_hyPar <- data.frame(data = character(), model = character(), iteration = numeric(), m = numeric(), min = numeric(), N = numeric(), sens = numeric(), uar = numeric())

for(i in 1:n_s){
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
  
  for(j in 1:nfolds){
    
    tune_train <- do.call(rbind, lapply(tune_d[-j], function(y){
      t4 <- y
      return(t4)
    }))
    tune_train <- balance.50(tune_train)
    
    tune_test <- do.call(rbind, tune_d[j])
    
    # create cluster
    cl <- makeCluster(numCores)
    
    # tuning hyperparamters
    hp <- expand.grid(m = 2^seq(0,6,1), 
                      min = 10^seq(0,5),
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
  hp_mean <- ddply(hp_cv, c("m","min","N"), summarise, sens.nfold = mean(sens), guar.nfold = mean(g.uar))
  hp_f <- hp_mean[which.max(hp_mean$guar.nfold),]
  
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
  xval[i,] <- data.frame(data = "RSV", model= "RF-NSx", iteration = i, data_balance = db, sensitivity = sen, specificity = spec, uar = uar)
  RF_hyPar[i,] <- c("RSV","RF-NSx",i, hp_f)                          # save hyperparameters
}

print("RF-NSx done")

###### PRINT RESULTS ######
xval
RF_hyPar

write.csv(xval, file = paste0("output/", Sys.Date(), ".RSV_DvsD_RF-NSx(xNIC).csv"))
write.csv(RF_hyPar, file = paste0("output/", Sys.Date(), ".RSV_DvsD_hyPar_RF-NSx(xNIC).csv"))

print("results printed")
