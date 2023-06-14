###### LOAD LIBRARIES ###### 
library(tidyverse)
library(plyr)
library(glmnet)
library(caret)
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

cv.lambda <- function(df,
                      nfolds = 10,
                      response = "IsCrossClass"){
  # folds
  fold <- createFolds(df[[response]],k=nfolds)
  
  # extract lambda range
  df.lamb <- balance.50(df)
  x <- model.matrix(formula(df.lamb), df.lamb)[,-1]
  y <- df.lamb[,1]
  y <- ifelse(y == FALSE,0,1)
  
  fit <- glmnet(x=x, 
                y=y,
                alpha  = 1, 
                family = binomial(),
                standardize = TRUE)
  lambda <- fit$lambda
  nlambda <- length(lambda)
  
  #---------------------
  # iteration over folds
  #---------------------
  
  # performance metrics
  sen <- data.frame()
  spec <- data.frame()
  uar <- data.frame()
  g.uar <- data.frame()
  
  for(i in 1:nfolds){
    
    # train set
    train <- do.call(rbind,lapply(fold[-i], function(x){
      tr <- df[x,]
      return(tr)
    }))
    train <- balance.50(train)
    
    # validation set
    val <- do.call(rbind,lapply(fold[i], function(x){
      t3 <- df[x,]
      return(t3)
    }))
    
    # set-up for LASSO
    tr_x <- model.matrix(formula(train), train)[,-1]
    
    tr_y <- train[,1]
    tr_y <- ifelse(tr_y == FALSE,0,1)
    
    val_x <- model.matrix(formula(val), val)[,-1]
    
    val_y <- val[,1]
    val_y <- ifelse(val_y == FALSE,0,1)
    
    for(l in 1:nlambda){
      fit.lambda <- glmnet(x=tr_x, 
                           y=tr_y,
                           alpha  = 1, 
                           lambda = lambda[l],
                           family = binomial(),
                           standardize = TRUE,
                           maxit = 1e6)
      p <- predict(fit.lambda, newx = val_x, type = "response")
      
      newp <- ifelse(p>0.5,1,0)
      newp <- factor(newp, levels = c("0","1"))
      Y <- factor(val_y, levels = c("0","1"))
      cm <- confusionMatrix(data = newp, 
                            reference = Y, 
                            positive = "1")
      
      # performance metrics
      sen[i,l] <- ifelse(is.na(cm$byClass[[1]]) == FALSE, cm$byClass[[1]] * 100,0)  # remove nas before calculating UAR
      spec[i,l] <- ifelse(is.na(cm$byClass[[2]]) == FALSE, cm$byClass[[2]] * 100,0)
      
      uar[i,l] <- 0.5*(sen[i,l] + spec[i,l])
      
      g.uar[i,l] <- sqrt(sen[i,l] * spec[i,l])
    }
    
    
  }
  g.uar[is.na(g.uar) == TRUE] <- 0
  
  #----------------------
  # save out lambda vals
  #----------------------
  sen.mean <- apply(sen, 2, mean)
  uar.mean <- apply(uar, 2, mean)
  uar.sd <- apply(uar,2,sd)
  uar.se <- uar.sd/sqrt(nfolds)
  
  g.uar.mean <- apply(g.uar, 2, mean)
  g.uar.sd <- apply(g.uar, 2, sd)
  g.uar.se <- g.uar.sd/sqrt(nfolds)
  
  lambda.cv <- data.frame(lambda = lambda, 
                          log_lambda = log(lambda),
                          mean_sen = sen.mean,
                          mean_uar = uar.mean,
                          se_uar = uar.se,
                          mean_g.uar = g.uar.mean,
                          se_g.uar = g.uar.se)
  
  
  #------------------------------------------
  # calculate lambda.max (uar) and lambda.1se
  #------------------------------------------
  where.max <- which.max(lambda.cv$mean_g.uar)
  lambda.max <- lambda.cv$lambda[where.max]
  
  lambda.max.se <- lambda.cv$se_g.uar[where.max]
  
  lambda.1se <- max(lambda.cv$lambda[
    lambda.cv$mean_g.uar > max(lambda.cv$mean_g.uar) - lambda.max.se
  ])
  
  # save values to return
  output <- list(lambda.1se,lambda.max,lambda.cv)
  
  # name list entries
  ch <- c("lambda.1se","lambda.max","lambda.cv")
  names(output) <- ch
  
  return(output)
}

print("functions loaded")

###### DEFINE ENVIRONMENT ###### 
set.seed(346)
numCores <- 6 - 1
setwd("/scratch/users/k21111073/")

print("environment defined")

###### PROCESS DATA ###### 
rsv <- read.delim("data/rsv-cg.txt")
covid <- read.delim("data/covid-cg.txt")

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
                        
                        # print models
                        cat("Iteration:", i)
                        
                        cat("RSV l_1se:", l_1se$rsv,"\n", "RSV coeffs:")
                        print(models$rsv$beta)
                        
                        cat("COVID l_1se:", l_1se$covid,"\n", "COVID coeffs:")
                        print(models$covid$beta)
                        
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

write.csv(xval, file = paste0("output/", Sys.Date(), ".CP-RSV_LASSO(xNIC).csv"))

print("results printed")

