###### LOAD LIBRARIES ###### 
library(tidyverse)
library(plyr)
library(caret)
library(foreach)
library(doParallel)
library(glmnet)

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
                      response = "IsCrossClass"){
  # folds
  nfolds <- length(df)
  
  # extract lambda range
  df.lamb <- do.call(rbind,df)
  df.lamb <- balance.50(df.lamb)
  x <- model.matrix(formula(df.lamb), df.lamb)[,-1]
  y <- df.lamb[,1]
  y <- ifelse(y == FALSE,0,1)
  
  fit <- glmnet(x=x, 
                y=y,
                alpha  = 1, 
                family = binomial(),
                standardize = TRUE,
                maxit = 1e6)
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
    train <- do.call(rbind, lapply(df[-i], function(y){
      tr <- y
      return(tr)
    }))
    train <- balance.50(train)
    
    # validation set
    val <- do.call(rbind, lapply(df[i], function(y){
      tr <- y
      return(tr)
    }))
    
    # set-up for LASSO
    tr_x <- model.matrix(formula(train), train)[,-1]
    
    tr_y <- train[,1]
    tr_y <- ifelse(tr_y == FALSE,0,1)
    
    val_x <- model.matrix(formula(val), val)[,-1]
    
    val_y <- val[,1]
    val_y <- ifelse(val_y == FALSE,0,1)
    
    # fit the lambdas
    fit.lambda <- glmnet(x=tr_x, 
                         y=tr_y,
                         alpha  = 1, 
                         lambda = lambda,
                         family = binomial(),
                         standardize = TRUE,
                         maxit = 1e6)
    p <- predict(fit.lambda, newx = val_x, type = "response")
    
    for(l in 1:nlambda){
      newp <- ifelse(p[,l]>0.5,1,0)
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
numCores <- 17 - 1
setwd("/scratch/users/k21111073/")

print("environment defined")

###### PROCESS DATA ###### 
rsv <- read.delim("data/rsv-cg.txt")
covid <- read.delim("data/covid-cg.txt")

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

# patient combinations is now (d-1)vs1 
n_s <- length(pid)                                                      # number of iterations

# parallelisation set-up
registerDoParallel(numCores)  

# starting loop in parallel
xval_rsv <- foreach(i=1:n_s,
                    .packages = c('caret', 'tidyverse', 'glmnet'),
                    .combine = rbind) %dopar% {
                      
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
                      
                      # print models
                      cat("Iteration:", i)
                      
                      cat("RSV l_1se:", l_1se,"\n", "RSV coeffs:")
                      print(model$beta)
                      
                      # test model
                      p <- predict(model, newx = ts_x, type = "response")
                      newp <- ifelse(p>0.5,1,0)
                      newp <- factor(newp, levels = c("0","1"))
                      Y <- factor(ts_y, levels = c("0","1"))
                      cm <- confusionMatrix(data = newp, 
                                            reference = Y, 
                                            positive = "1")
                      
                      # committing performance metrics to data.frame
                      sen <- cm$byClass[[1]] * 100                           # sensitivity
                      spec <- cm$byClass[[2]] * 100                          # specificity
                      uar <- 0.5*(cm$byClass[[1]] + cm$byClass[[2]]) * 100   # unweighted accuracy recall
                      g.uar <- sqrt(cm$byClass[[1]] * cm$byClass[[2]]) * 100 # geometric uar
                      
                      xval <- data.frame(model = "RSV", iteration = i, lambda = l_1se, data_balance = db, sensitivity = sen, specificity = spec, uar = uar, g.uar = g.uar) 
                      
                      return(xval)
                      
                    }

stopImplicitCluster()

print("RSV done")

###### COVID ######
# patient separation
pid <- unique(covid$PatientID)
covid <- lapply(pid, function(x){
  p <- subset(covid, PatientID == x, select = -c(PatientID))
})

# patient combinations is now (d-1)vs1 
n_s <- length(pid)                                                      # number of iterations

# parallelisation set-up
registerDoParallel(numCores)  

# starting loop in parallel
xval_covid <- foreach(i=1:n_s,
                      .packages = c('caret', 'tidyverse', 'glmnet'),
                      .combine = rbind) %dopar% {
                        
                        # create train set with (d-1) patients
                        train <- lapply(covid[-i], function(y){
                          tr <- y
                          return(tr)
                        })
                        
                        # create test set with remaining patient
                        test <- do.call("rbind", covid[i])
                        
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
                        
                        # print models
                        cat("Iteration:", i)
                        
                        cat("COVID l_1se:", l_1se,"\n", "COVID coeffs:")
                        print(model$beta)
                        
                        # test model
                        p <- predict(model, newx = ts_x, type = "response")
                        newp <- ifelse(p>0.5,1,0)
                        newp <- factor(newp, levels = c("0","1"))
                        Y <- factor(ts_y, levels = c("0","1"))
                        cm <- confusionMatrix(data = newp, 
                                              reference = Y, 
                                              positive = "1")
                        
                        # committing performance metrics to data.frame
                        sen <- cm$byClass[[1]] * 100                           # sensitivity
                        spec <- cm$byClass[[2]] * 100                          # specificity
                        uar <- 0.5*(cm$byClass[[1]] + cm$byClass[[2]]) * 100   # unweighted accuracy recall
                        g.uar <- sqrt(cm$byClass[[1]] * cm$byClass[[2]]) * 100 # geometric uar
                        
                        xval <- data.frame(model = "COVID", iteration = i, lambda = l_1se, data_balance = db, sensitivity = sen, specificity = spec, uar = uar, g.uar = g.uar) 
                        
                        return(xval)
                        
                      }

stopImplicitCluster()

print("COVID done")

###### PRINT RESULTS ######
write.csv(xval_rsv, paste0("output/", Sys.Date(), ".RSV_DvsD_LASSO(xNIC).csv"), row.names = FALSE)
write.csv(xval_covid, paste0("output/", Sys.Date(), ".CP_DvsD_LASSO(xNIC).csv"), row.names = FALSE)

xval_rsv
xval_covid

print("results printed")
