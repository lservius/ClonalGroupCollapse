# hyperparameter search for LASSO lambda using L1DO strategy and the maximising the geometric mean for selection

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