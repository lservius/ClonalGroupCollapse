###### LOAD LIBRARIES ###### 
library(tidyverse)
library(plyr)
library(caret)
library(foreach)
library(doParallel)
print("libraries loaded")

###### FUNCTIONS ###### 
predict.TF <- function(prediction, threshold = 0.5){
  predict_TF <- NULL
  for (i in 1:length(prediction)){
    if (is.na(prediction[i]) == TRUE){
      predict_TF[i] <- NA
    } else if (prediction[i] > threshold){
      predict_TF[i] <- TRUE
    } else {
      predict_TF[i] <-  FALSE
    }
  }
  return(predict_TF)
}

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
numCores <- 17 - 1
setwd("/scratch/users/k21111073/")

print("environment defined")

###### PROCESS DATA ###### 
covid <- read.delim("data/covid-cg.txt")

# remove unnecessary cols
covid <- select(covid, -c(CloneID))

print("data processed")

###### MODEL SET-UP ###### 
# make characters into factors
cat <- c("IsCrossClass", "Vgene", "Jfamily", "Dfamily", "Dgene", "Sex")
covid[cat] <- lapply(covid[cat], factor)

# patient separation
pid <- unique(covid$PatientID)
covid <- lapply(pid, function(x){
  p <- subset(covid, PatientID == x, select = -c(PatientID))
})

# patient combinations is now (d-1)vs1 
n_s <- length(pid)  

###### LASSO ###### 
# parallelisation set-up
registerDoParallel(numCores)  

xval <- data.frame(model = "LR", data = "COVID", iteration = 1, data_balance = 1, sensitivity = 1, specificity = 1, uar = 1)


xval_lr <- foreach(i=1:n_s,
                   .packages = c('caret','tidyverse','glmnet')) %dopar% {
                     # create list of training sets
                     train <- do.call(rbind, lapply(covid[-i], function(y){
                       tr <- y
                       tr <- balance.50(tr)  # balance dataset
                       return(tr)
                     }))
                     
                     # create test set with remaining patient
                     test <- do.call("rbind", covid[i])
                     
                     db <- sum(train$IsCrossClass == TRUE)/nrow(train) * 100    # save out balance in response for train data
                     
                     initial <- glm(IsCrossClass~1,
                                    family=binomial(),
                                    data=train)
                     
                     model <- step(initial,
                                   scope=formula(train),
                                   direction="both")
                     
                     # print model
                     cat("Iteration:", i, "\n", "COVID model summary:","\n")
                     print(summary(model))
                     
                     # getting expected response, y, using the trained model
                     y <- predict(model, test, type = "response")
                     
                     # getting performance
                     y <- factor(predict.TF(y))
                     pred <- merge(y, test$IsCrossClass, by = "row.names")
                     cm <- confusionMatrix(data = y, reference = test$IsCrossClass, positive = "TRUE")
                     
                     # performance metrics
                     sen <- cm$byClass[[1]] * 100
                     
                     spec <- cm$byClass[[2]] * 100
                     
                     uar <- 0.5*(cm$byClass[[1]] + cm$byClass[[2]]) * 100
                     
                     # save out xval
                     xval <- data.frame(model = "LR", data = "COVID", iteration = i, data_balance = db, sensitivity = sen, specificity = spec, uar = uar)
                     
                     return(xval)
                   }

stopImplicitCluster()

xval <- xval_lr

print("LR done")

###### PRINT RESULTS ######
xval

write.csv(xval, file = paste0("output/", Sys.Date(), ".CP_DvsD_LR(xNIC).csv"))

print("results printed")

