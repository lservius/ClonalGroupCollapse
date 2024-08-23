###### LOAD LIBRARIES ###### 
library(tidyverse)
library(plyr)
library(glmnet)
library(caret)
library(foreach)
library(doParallel)
print("libraries loaded")

###### FUNCTIONS ###### 
source("../functions/func_data-balance.R")
source("../functions/func_LASSO-CV-AD.R")


print("functions loaded")

###### DEFINE ENVIRONMENT ###### 
set.seed(346)
numCores <- 6 - 1

print("environment defined")

###### PROCESS DATA ###### 
rsv <- read.delim("../data/processed/rsv-cg.txt")
covid <- read.delim("../data/processed/covid-cg.txt")

# all together now
ch <- c("covid","rsv")
data <- list(covid,rsv)
names(data) <- ch

# remove pid
data <- lapply(data, select, -c(PatientID,CloneID))

# create folds
n_k <- 5 
folds <- lapply(data,function(x){
  fold <- createFolds(x$IsCrossClass,k=n_k)
})

print("data processed")


###### MODEL SET-UP ###### 
# make characters into factors
cat <- c("IsCrossClass", "Vgene", "Jfamily", "Dfamily", "Dgene", "Age", "Sex")
data <- lapply(data, function(x){
  x[cat] <- lapply(x[cat], factor)
  return(x)
})

# transform COVID age
data$covid$Age <- as.integer(data$covid$Age)

# remove sex from RSV
data$rsv <- select(data$rsv, -c(Sex))

###### LASSO ###### 
# parallelisation set-up
registerDoParallel(numCores)  

xval <- list(covid = data.frame(data = "data", iteration = 1, lambda = 1, data_balance = 1, sensitivity = 1, specificity = 1, uar = 1, g.uar = 1),
             rsv = data.frame(data = "data", iteration = 1, lambda = 1, data_balance = 1, sensitivity = 1, specificity = 1, uar = 1, g.uar = 1)
)


xval_lasso <- foreach(i=1:n_k,
                      .packages = c('caret','tidyverse','glmnet')) %dopar% {
                        # create list of training sets
                        train <- mapply(function(x,y){
                          tr <- do.call(rbind,lapply(y[-i], function(z){
                            tr <- x[z,]
                            return(tr)
                          }))
                          return(tr)
                        },data,folds,SIMPLIFY = FALSE)
                        
                        # create list of testing sets
                        test <- mapply(function(x,y){
                          tr <- do.call(rbind,lapply(y[i], function(z){
                            t3 <- x[z,]
                          }))
                        },data,folds,SIMPLIFY = FALSE)
                        
                        # test set for LASSO
                        test_x <- lapply(test,function(x){
                          te_x <- model.matrix(formula(x), x)[,-1]
                        })
                        
                        test_y <- lapply(test,function(x){
                          te_y <- x[,1]
                          te_y <- ifelse(te_y == FALSE,0,1)
                        })
                        
                        # find the amount of penalty
                        cv_lambda <- lapply(train,function(x){
                          cv_lambda <- cv.lambda(x,
                                                 nfolds = 10,
                                                 response = "IsCrossClass")
                          return(cv_lambda)
                        })
                        
                        # extract best fit lambda
                        l_1se <- lapply(cv_lambda, function(x){x$lambda.1se})
                        
                        # balance train set
                        train <- lapply(train,balance.50)
                        
                        db <- lapply(train, function(x){
                          sum(x$IsCrossClass == TRUE)/nrow(x) * 100
                        })
                        
                        # train for LASSO
                        train_x <- lapply(train,function(x){
                          tr_x <- model.matrix(formula(x), x)[,-1]
                          return(tr_x)
                        })
                        
                        train_y <- lapply(train,function(x){
                          tr_y <- x[,1]
                          tr_y <- ifelse(tr_y == FALSE,0,1)
                        })
                        
                        # fit models
                        models <- mapply(function(X,Y,Z){
                          glmnet(x=X, 
                                 y=Y,
                                 alpha  = 1, 
                                 lambda = Z,
                                 family = binomial(),
                                 standardize = TRUE)
                        },train_x,train_y,l_1se,SIMPLIFY = FALSE)
                        
                        # save model
                        write_rds(models$rsv, file = paste0("../results/rLASSO_model",i,".rda"))
                        write_rds(models$covid, file = paste0("../results/cLASSO_model",i,".rda"))
                        
                        # test models
                        cm <- mapply(function(M,X,Y){
                          p <- predict(M, newx = X, type = "response")
                          newp <- ifelse(p>0.5,1,0)
                          newp <- factor(newp, levels = c("0","1"))
                          Y <- factor(Y, levels = c("0","1"))
                          cm <- confusionMatrix(data = newp, 
                                                reference = Y, 
                                                positive = "1")
                          return(cm)
                        },models,test_x,test_y,SIMPLIFY = FALSE)
                        
                        # performance metrics
                        sen <- lapply(cm, function(x){
                          x$byClass[[1]] * 100
                        })
                        
                        spec <- lapply(cm, function(x){
                          x$byClass[[2]] * 100
                        })
                        
                        uar <- lapply(cm, function(x){
                          0.5*(x$byClass[[1]] + x$byClass[[2]]) * 100
                        })
                        
                        g.uar <- lapply(cm, function(x){
                          sqrt(x$byClass[[1]] * x$byClass[[2]]) * 100
                        })
                        
                        # save out xval
                        xval$covid <- data.frame(data = "COVID", iteration = i, lambda = l_1se$covid, data_balance = db$covid, sensitivity = sen$covid, specificity = spec$covid, uar = uar$covid, g.uar = g.uar$covid)
                        
                        xval$rsv <- data.frame(data = "RSV", iteration = i, lambda = l_1se$rsv, data_balance = db$rsv, sensitivity = sen$rsv, specificity = spec$rsv, uar = uar$rsv, g.uar = g.uar$rsv)
                        
                        return(xval)
                      }

stopImplicitCluster()

xval <- do.call('rbind',lapply(xval_lasso, function(x){
  do.call('rbind',x)
}))


print("LASSO done")

###### PRINT RESULTS ######
xval

write.csv(xval, file = "../results/CP-RSV_LASSO.csv")

print("results printed")

