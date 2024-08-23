###### LOAD LIBRARIES ###### 
library(tidyverse)
library(plyr)
library(lme4)
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
numCores <- 6 - 1

print("environment defined")

###### PROCESS DATA ###### 
covid <- read.delim("../data/processed/covid-cg.txt")

# remove pid
resp <-select(covid, c(IsCrossClass, PatientID))
covid <- select(covid, -c(PatientID,CloneID,IsCrossClass))

print("data processed")

###### MODEL SET-UP ###### 
# scale numeric
scale_num <- names(covid)[append(5,7:length(covid))]
covid[scale_num] <- lapply(covid[scale_num], scale)

# hot-encode variable region variables and metadata
covid <- as.data.frame(model.matrix(~., data = covid))[-1]
covid <- cbind(resp,covid)

# make hotencode into integers
cat <- append(names(covid)[3:85], names(covid)[87])
covid[cat] <- lapply(covid[cat], as.integer)

# create folds
n_k <- 5 
fold <- createFolds(covid$IsCrossClass,k=n_k)

###### GLMM ###### 
# parallelisation set-up
registerDoParallel(numCores)  

xval <- data.frame(model = "LR", data = "COVID", iteration = 1, data_balance = 1, sensitivity = 1, specificity = 1, uar = 1, g.uar = 1, AUC_ROC = 1, Precision = 1, Recall = 1, F1 = 1) 

system.time(
  xval_glmm <- foreach(i=1:n_k,
                       .packages = c('caret','tidyverse','glmnet')) %dopar% {
                         # create list of training sets
                         train <- do.call(rbind, lapply(fold[-i], function(y){
                           tr <- covid[y,]
                           tr <- balance.50(tr)
                           return(tr)
                         }))
                         
                         # create test set with 1 fold
                         test <- do.call(rbind, lapply(fold[i], function(y){
                           t <- covid[y,]
                           return(t)
                         }))
                         
                         db <- sum(train$IsCrossClass == TRUE)/nrow(train) * 100    # save out balance in response for train data
                         
                         
                         # remove PatientID for step
                         train_lr <- select(train, -c(PatientID))
                         
                         initial <- glm(IsCrossClass~1,
                                        family=binomial(),
                                        data=train_lr)
                         
                         cat("begin step function....")
                         
                         model_lr <- step(initial,
                                          scope=formula(train_lr),
                                          direction="both",
                                          trace = FALSE)
                         
                         # formula
                         fx <- formula(model_lr)
                         fx <- update(fx,    ~ . + (1|PatientID))
                         
                         cat("begin fitting....")
                         
                         # use step LR selection for GLMM
                         model <- glmer(fx,
                                        family = binomial(), 
                                        data = train)
                         
                         # print model
                         cat("Iteration:", i, "\n")
                         beta <- coef(model)
                         saveRDS(model, file = paste0(Sys.Date(),".i",i, ".COVID_GLMM.rda"))
                         
                         # getting expected response, y, using the trained model
                         y <- predict(model, test, type = "response", allow.new.levels = TRUE)
                         
                         
                         # auc
                         resp <- test$IsCrossClass
                         roc <- roc(response = resp, predictor = y)
                         auc_roc <- auc(roc)
                         
                         # getting performance
                         y <- factor(predict.TF(y))
                         ref <- factor(ifelse(test$IsCrossClass == 1, TRUE, FALSE))
                         cm <- confusionMatrix(data = y, reference = ref, positive = "TRUE")
                         
                         # performance metrics
                         sen <- cm$byClass[[1]]
                         spec <- cm$byClass[[2]]
                         uar <- 0.5*(cm$byClass[[1]] + cm$byClass[[2]])
                         g.uar <- sqrt(cm$byClass[[1]] * cm$byClass[[2]]) # geometric uar
                         
                         # save out xval
                         xval <- data.frame(model = "LR", data = "COVID", iteration = i,  data_balance = db, sensitivity = sen, specificity = spec, uar = uar, g.uar = g.uar, AUC_ROC = auc_roc, Precision = cm$byClass[[5]], Recall = cm$byClass[[6]], F1 = cm$byClass[[7]]) 
                         
                         return(xval)
                       }
)
stopImplicitCluster()

xval <- do.call(rbind, xval_glmm)

print("GLMM done")

###### PRINT RESULTS ######
xval
write.csv(xval, file = "../results/COVID_GLMM.csv")

print("results printed")