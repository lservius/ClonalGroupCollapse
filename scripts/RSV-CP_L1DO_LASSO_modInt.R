###### LOAD LIBRARIES ###### 
library(tidyverse)
library(plyr)
library(caret)
library(pROC)
library(foreach)
library(doParallel)
library(glmnet)

print("libraries loaded")

###### FUNCTIONS ###### 
source("../functions/func_data-balance.R")
source("../functions/func_LASSO-CV-L1DO.R")
print("functions loaded")

###### DEFINE ENVIRONMENT ###### 
set.seed(346)
numCores <- 17 - 1
i <- 3
j <- 14

print("environment defined")

###### PROCESS DATA ###### 
rsv <- read.delim("../data/processed/rsv-cg.txt")
covid <- read.delim("../data/processed/covid-cg.txt")

# all together now
ch <- c("covid","rsv")
data <- list(covid,rsv)
names(data) <- ch

# remove constant/unique cols
data <- lapply(data, select, -c(CloneID))

print("data processed")

###### MODEL SET-UP ###### 
# make characters into factors
cat <- c("IsCrossClass", "Vgene", "Jfamily", "Dfamily", "Dgene", "Age","Sex")
data <- lapply(data, function(x){
  x[cat] <- lapply(x[cat], factor)
  return(x)
})

# transform COVID age
data$covid$Age <- as.integer(data$covid$Age)

# remove RSV sex
data$rsv <- select(data$rsv, -c("Sex"))

# separate back into different df
rsv <- data$rsv
covid <- data$covid

###### RSV ######
# patient separation
pid <- unique(rsv$PatientID)
rsv <- lapply(pid, function(x){
  p <- subset(rsv, PatientID == x, select = -c(PatientID))
})


# create train set with (d-1) patients
train <- lapply(rsv[-i], function(y){
  tr <- y
  return(tr)
})

# create test set with remaining patient
test <- do.call("rbind", rsv[i])

# find the amount of penalty
cv_lambda <- cv.lambda(train)

# extract best fit lambda
l_1se <- cv_lambda$lambda.1se

# bind training set
train <- do.call(rbind, train)
train <- balance.50(train)
db <- sum(train$IsCrossClass == TRUE)/nrow(train)                               # save out balance in response for train data

# train for LASSO
tr_x <- model.matrix(formula(train), train)[,-1]
tr_y <- train[,1]
tr_y <- ifelse(tr_y == FALSE,0,1)


# test set for LASSO
ts_x <- model.matrix(formula(test), test)[,-1]
ts_y <- test[,1]
ts_y <- ifelse(ts_y == FALSE,0,1)

# fit model
model <- glmnet(x=tr_x, 
                y=tr_y,
                alpha  = 1, 
                lambda = l_1se,
                family = binomial(),
                standardize = TRUE)


saveRDS(model, "RSV_modInt_LASSO.rda")


# test model
p <- predict(model, newx = ts_x, type = "response")

# auc
resp <- ifelse(test$IsCrossClass == TRUE, 1, 0)
roc <- roc(response = resp, predictor = p[,1])
auc_roc <- auc(roc)

newp <- ifelse(p>0.5,1,0)
newp <- factor(newp, levels = c("0","1"))
Y <- factor(ts_y, levels = c("0","1"))
cm <- confusionMatrix(data = newp, 
                      reference = Y, 
                      positive = "1")
# committing performance metrics to data.frame
sen <- cm$byClass[[1]]                           # sensitivity
spec <- cm$byClass[[2]]                          # specificity
uar <- 0.5*(cm$byClass[[1]] + cm$byClass[[2]])   # unweighted accuracy recall
g.uar <- sqrt(cm$byClass[[1]] * cm$byClass[[2]]) # geometric uar

xval_rsv <- data.frame(data = "RSV", model = "LASSO LR", iteration = i, lambda = l_1se, data_balance = db, sensitivity = sen, specificity = spec, uar = uar, g.uar = g.uar, AUC_ROC = auc_roc, Precision = cm$byClass[[5]], Recall = cm$byClass[[6]], F1 = cm$byClass[[7]]) 

###### PRINT RESULTS ######
write.csv(xval_rsv, "../results/figures_tables/table3_RSV_LASSO.csv", row.names = FALSE)
xval_rsv

print("RSV done")

###### COVID ######
# patient separation
pid <- unique(covid$PatientID)
covid <- lapply(pid, function(x){
  p <- subset(covid, PatientID == x, select = -c(PatientID))
})

# create train set with (d-1) patients
train <- lapply(covid[-j], function(y){
  tr <- y
  return(tr)
})

# create test set with remaining patient
test <- do.call("rbind", covid[j])

# find the amount of penalty
cv_lambda <- cv.lambda(train)

# extract best fit lambda
l_1se <- cv_lambda$lambda.1se

# bind training set
train <- do.call(rbind, train)
train <- balance.50(train)
db <- sum(train$IsCrossClass == TRUE)/nrow(train)                               # save out balance in response for train data

# train for LASSO
tr_x <- model.matrix(formula(train), train)[,-1]
tr_y <- train[,1]
tr_y <- ifelse(tr_y == FALSE,0,1)


# test set for LASSO
ts_x <- model.matrix(formula(test), test)[,-1]
ts_y <- test[,1]
ts_y <- ifelse(ts_y == FALSE,0,1)

# fit model
model <- glmnet(x=tr_x, 
                y=tr_y,
                alpha  = 1, 
                lambda = l_1se,
                family = binomial(),
                standardize = TRUE)

saveRDS(model, "COVID_modInt_LASSO.rda")

# test model
p <- predict(model, newx = ts_x, type = "response")
newp <- ifelse(p>0.5,1,0)
newp <- factor(newp, levels = c("0","1"))
Y <- factor(ts_y, levels = c("0","1"))
cm <- confusionMatrix(data = newp, 
                      reference = Y, 
                      positive = "1")

# auc
resp <- ifelse(test$IsCrossClass == TRUE, 1, 0)
roc <- roc(response = resp, predictor = p[,1])
auc_roc <- auc(roc)

# committing performance metrics to data.frame
sen <- cm$byClass[[1]]                           # sensitivity
spec <- cm$byClass[[2]]                          # specificity
uar <- 0.5*(cm$byClass[[1]] + cm$byClass[[2]])   # unweighted accuracy recall
g.uar <- sqrt(cm$byClass[[1]] * cm$byClass[[2]]) # geometric uar

xval_covid <- data.frame(data = "COVID", model = "LASSO LR", iteration = j, lambda = l_1se, data_balance = db, sensitivity = sen, specificity = spec, uar = uar, g.uar = g.uar, AUC_ROC = auc_roc, Precision = cm$byClass[[5]], Recall = cm$byClass[[6]], F1 = cm$byClass[[7]]) 


print("COVID done")

###### PRINT RESULTS ######
write.csv(xval_covid, "../results/figures_tables/table3_CP_LASSO.csv", row.names = FALSE)

xval_covid

print("results printed")
