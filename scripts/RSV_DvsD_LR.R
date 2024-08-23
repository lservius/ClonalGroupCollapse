###### LOAD LIBRARIES ###### 
library(tidyverse)
library(plyr)
library(pROC)
library(caret)
library(foreach)
library(doParallel)
print("libraries loaded")

###### FUNCTIONS ###### 
source("../functions/func_predictThreshold.R")
source("../functions/func_data-balance.R")
print("functions loaded")

###### DEFINE ENVIRONMENT ###### 
set.seed(346)
numCores <- 13 - 1
setwd("/scratch/users/k21111073/")

print("environment defined")

###### PROCESS DATA ###### 
rsv <- read.delim("../data/processed/rsv-cg.txt")

#separate out PatientID
pid_col <- select(rsv, PatientID,IsCrossClass)
# remove unnecessary cols
rsv <- select(rsv, -c(CloneID,Sex,PatientID,IsCrossClass))

print("data processed")


###### MODEL SET-UP ###### 
# scale numeric
scale_num <- names(rsv)[6:64]
rsv[scale_num] <- lapply(rsv[scale_num], scale)

# hot-encode variable region variables and metadata
rsv <- as.data.frame(model.matrix(~., data = rsv))[-1]

# make indicator vars into integers
cat <- names(rsv)[1:84]
rsv[cat] <- lapply(rsv[cat], as.integer)

# add back patientID and response
rsv <- cbind(pid_col,rsv)

# patient separation
pid <- unique(rsv$PatientID)
rsv <- lapply(pid, function(x){
  p <- subset(rsv, PatientID == x, select = -c(PatientID))
})

# patient combinations is now (d-1)vs1 
n_s <- length(pid)  

###### LASSO ###### 
# parallelisation set-up
registerDoParallel(numCores)  

xval <- data.frame(model = "LR", data = "RSV", iteration = 1, data_balance = 1, sensitivity = 1, specificity = 1, uar = 1, g.uar = 1, AUC_ROC = 1, Precision = 1, Recall = 1, F1 = 1) 

system.time(
xval_lr <- foreach(i=1:n_s,
                   .packages = c('caret','tidyverse','glmnet')) %dopar% {
                     # create list of training sets
                     train <- do.call(rbind, lapply(rsv[-i], function(y){
                       tr <- y
                       tr <- balance.50(tr)  # balance dataset
                       return(tr)
                     }))
                     
                     # create test set with remaining patient
                     test <- do.call("rbind", rsv[i])
                     
                     db <- sum(train$IsCrossClass == TRUE)/nrow(train) * 100    # save out balance in response for train data
                     
                     initial <- glm(IsCrossClass~1,
                                    family=binomial(),
                                    data=train)
                     
                     model <- step(initial,
                                   scope=formula(train),
                                   trace = 0,
                                   direction="both")

                     # getting expected response, y, using the trained model
                     y <- predict(model, test, type = "response")
                     
                     # auc
                     resp <- ifelse(test$IsCrossClass == TRUE, 1, 0)
                     roc <- roc(response = resp, predictor = y)
                     auc_roc <- auc(roc)
                     
                     # getting performance
                     y <- factor(predict.TF(y))
                     ref <- factor(test$IsCrossClass)
                     cm <- confusionMatrix(data = y, reference = ref, positive = "TRUE")
                     
                     # performance metrics
                     sen <- cm$byClass[[1]]
                     
                     spec <- cm$byClass[[2]]
                     
                     uar <- 0.5*(cm$byClass[[1]] + cm$byClass[[2]])
                     
                     g.uar <- sqrt(cm$byClass[[1]] * cm$byClass[[2]]) # geometric uar
                     
                     # save out xval
                     xval <- data.frame(model = "LR", data = "RSV", iteration = i,  data_balance = db, sensitivity = sen, specificity = spec, uar = uar, g.uar = g.uar, AUC_ROC = auc_roc, Precision = cm$byClass[[5]], Recall = cm$byClass[[6]], F1 = cm$byClass[[7]]) 
                     
                     return(xval)
                   }
)
stopImplicitCluster()

xval <- xval_lr

print("LR done")

###### PRINT RESULTS ######
xval
write.csv(xval, file = "../results/RSV_DvsD_LR.csv")

print("results printed")

