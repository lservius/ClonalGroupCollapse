library(plyr)
library(tidyverse)
library(cowplot)
library(randomForest)
library(caret)
library(dplyr)

setwd("results/")

# collapse iterations into 1 df
allin1 <- function(rootName,
                   n,
                   dirPath = ".",
                   noSim = "NA"){
  files <- list.files(path = dirPath)
  files <- files[grepl(rootName, files)]
  
  if(noSim == "NA"){
    df <- do.call(rbind,lapply(1:n,function(i){
      read.csv(files[i])
    }))
  } else {
    files <- files[!(files %in% files[grepl(noSim, files)])]
    df <- do.call(rbind,lapply(1:n,function(i){
      read.csv(files[i])
    }))
  }
  
  df <- na.omit(df)
  df <- df[order(df$iteration),]
  return(df)
}

## ALL DONOR RESULTS
LASSO.AD <- read.csv("CP-RSV_LASSO.csv")
R_LR.AD <- read.csv("RSV_LR.csv")
C_LR.AD <- read.csv("CP_LR.csv")
R_RF.AD <- read.csv("RSV_RF.csv")
C_RF.AD <- read.csv("CP_RF.csv")

R_RFx.AD <- read.csv("RSV_RF-NSx.csv")
C_RFx.AD <- read.csv("CP_RF-NSx.csv")

R_SVM.AD <- allin1(rootName = "RSV_SVM", n = 5)
C_SVM.AD <- allin1(rootName = "CP_SVM", n = 5)

R_LASSO.AD <- subset(LASSO.AD, data == "RSV")
C_LASSO.AD <- subset(LASSO.AD, data == "COVID")

R_GLMM.AD <- read.csv("RSV_GLMM.csv")
C_GLMM.AD <- read.csv("COVID_GLMM.csv")

### LEAVE ONE DONOR OUT RESULTS
R_LASSO.DD <- read.csv("RSV_DvsD_LASSO.csv")
C_LASSO.DD <- read.csv("CP_DvsD_LASSO.csv")

R_LR.DD <- read.csv("RSV_DvsD_LR.csv")
C_LR.DD <- read.csv("CP_DvsD_LR.csv")

R_RF.DD <- read.csv("RSV_DvsD_RF.csv")
C_RF.DD <- read.csv("CP_DvsD_RF.csv")

R_RFx.DD <- read.csv("RSV_DvsD_RF-NSx.csv")
C_RFx.DD <- read.csv("CP_DvsD_RF-NSx.csv")

R_SVM.DD <- allin1(rootName = "RSV_DvsD_SVM", n = 12, noSim = "SVMhyPar")
C_SVM.DD <- allin1(rootName = "CP_DvsD_SVM", n = 16, noSim = "SVMhyPar")

R_GLMM.DD <- read.csv("RSV_L1DO_GLMM.csv")
C_GLMM.DD <- read.csv("COVID_L1DO_GLMM.csv")


## TRAINFIT RF RESULTS
R_RF_trainfit.AD <- read.csv("RSV_RF_trainfit.csv")
C_RF_trainfit.AD <- read.csv("CP_RF_trainfit.csv")

R_RFx_trainfit.AD <- read.csv("RSV_RF-NSx_trainfit.csv")
C_RFx_trainfit.AD <- read.csv("CP_RF-NSx_trainfit.csv")

R_RF_trainfit.DD <- read.csv("RSV_DvsD_RF_trainfit.csv")
C_RF_trainfit.DD <- read.csv("CP_DvsD_RF_trainfit.csv")

R_RFx_trainfit.DD <- read.csv("RSV_DvsD_RF-NSx_trainfit.csv")
C_RFx_trainfit.DD <- read.csv("CP_DvsD_RF-NSx_trainfit.csv")


# make uniform
R_GLMM.AD$model <- "GLMM"
C_GLMM.AD$model <- "GLMM"

R_LASSO.AD$model <- "LASSO LR"
C_LASSO.AD$model <- "LASSO LR"

R_LASSO.DD$data <- "RSV"
C_LASSO.DD$data <- "COVID"

R_LASSO.DD$model <- "LASSO LR"
C_LASSO.DD$model <- "LASSO LR"

R_LR.AD[,6:13] <- R_LR.AD[,6:13] * 100
C_LR.AD[,6:13] <- C_LR.AD[,6:13] * 100

R_LR.DD[,6:13] <- R_LR.DD[,6:13] * 100
C_LR.DD[,6:13] <- C_LR.DD[,6:13] * 100

R_GLMM.AD[,6:13] <- R_GLMM.AD[,6:13] * 100
C_GLMM.AD[,6:13] <- C_GLMM.AD[,6:13] * 100

R_GLMM.DD[,6:13] <- R_GLMM.DD[,6:13] * 100
C_GLMM.DD[,6:13] <- C_GLMM.DD[,6:13] * 100



# make lists of performance for each test per dataset
RSV.AD <- list(R_LASSO.AD,R_LR.AD,R_RF.AD
               ,R_RFx.AD
               ,R_SVM.AD
               ,R_GLMM.AD[1:8]
)
COVID.AD <- list(C_LASSO.AD,C_LR.AD,C_RF.AD
                 ,C_RFx.AD
                 ,C_SVM.AD
                 ,C_GLMM.AD[1:8]
)

RSV.DD <- list(R_LASSO.DD,R_LR.DD
               ,R_RF.DD
               ,R_RFx.DD
               ,R_SVM.DD
               ,R_GLMM.DD[1:8]
)
COVID.DD <- list(C_LASSO.DD, C_LR.DD
                 ,C_RF.DD
                 , C_RFx.DD
                 ,C_SVM.DD
                 ,C_GLMM.DD[1:8]
)

# calculate gmeans, precision and F1
calc_additional_revResp <- function(x,PN){
  sum <- lapply(x, function(y){
    cols <- names(y)
    
    # set vars
    y$TP <- NA
    y$TN <- NA
    y$FP <- NA
    
    # order based on iteration
    y <- y[order(y$iteration),]
    
    TP <- PN$Positives * y$sensitivity/100
    TN <- PN$Negatives * y$specificity/100
    FP <- TN * ( (y$specificity/100) **(-1) - 1 )
    
    y$TP <- TP
    y$FP <- FP
    y$FN <- PN$Positives - TP
    
    ddply(y, cols, summarise, gmean = sqrt(sensitivity * specificity), precision = TP / (TP + FP) * 100, F1 = 2*TP / (2*TP + FP + FN) * 100)
  })
  
  return(sum)
}

# import positives/negatives
cPN_AD <- read.csv("CP_AD_pn.csv")
cPN_DD <- read.csv("CP_L1DO_pn.csv")

rPN_AD <- read.csv("RSV_AD_pn.csv")
rPN_DD <- read.csv("RSV_L1DO_pn.csv")

RSV.AD <- calc_additional_revResp(RSV.AD,rPN_AD)
COVID.AD <- calc_additional_revResp(COVID.AD,cPN_AD)

RSV.DD <- calc_additional_revResp(RSV.DD,rPN_DD)
COVID.DD <- calc_additional_revResp(COVID.DD,cPN_DD)

## SUMMARISE PERFORMANCE
sumPerf <- function(x){
  sum <- lapply(x, function(y){
    ddply(y, c("data", "model"), summarise, db = mean(data_balance), sens.mu = mean(sensitivity), sens.se = sd(sensitivity)/sqrt(length(sensitivity)), spec.mu = mean(specificity), spec.se = sd(specificity)/sqrt(length(specificity)),prec.mu = mean(precision), prec.se = sd(precision)/sqrt(length(precision)), F1.mu = mean(F1), F1.se = sd(F1)/sqrt(length(F1)), uar.mu = mean(uar), uar.se = sd(uar)/sqrt(length(uar)), gmean.mu = mean(gmean), gmean.se = sd(gmean)/sqrt(length(gmean)))
  })
  
  sumCondense <- do.call(rbind, sum)
  return(sumCondense)
}

RSV.AD <- sumPerf(RSV.AD)
COVID.AD <- sumPerf(COVID.AD)

RSV.DD <- sumPerf(RSV.DD)
COVID.DD <- sumPerf(COVID.DD)


AD <- rbind(
  RSV.AD,
  COVID.AD)

L1DO <- rbind(
  RSV.DD,
  COVID.DD)

# Table S2/3
write.csv(AD, "figures_tables/tableS2.csv")
write.csv(L1DO, "figures_tables/tableS3.csv")


## TRAINFIT
# bind same CV and data, different model
cRFs_AD_trainfit <- rbind(C_RF_trainfit.AD,C_RFx_trainfit.AD)
cRFs_L1DO_trainfit <- rbind(C_RF_trainfit.DD,C_RFx_trainfit.DD)

rRFs_AD_trainfit <- rbind(R_RF_trainfit.AD,R_RFx_trainfit.AD)
rRFs_L1DO_trainfit <- rbind(R_RF_trainfit.DD,R_RFx_trainfit.DD)


# make uniform
cRFs_AD_trainfit$CV <- "AD"
rRFs_AD_trainfit$CV <- "AD"

cRFs_L1DO_trainfit$CV <- "L1DO"
rRFs_L1DO_trainfit$CV <- "L1DO"

# bind same data
cRFs_trainfit <- rbind(cRFs_AD_trainfit, cRFs_L1DO_trainfit)
rRFs_trainfit <- rbind(rRFs_AD_trainfit, rRFs_L1DO_trainfit)

# summarise training fit performance
sumPerfTrain <- function(x){
  sum <- ddply(x, c("data","CV","model"), summarise, sens.mu = mean(sensitivity), sens.se = sd(sensitivity)/sqrt(length(sensitivity)), spec.mu = mean(specificity), spec.se = sd(specificity)/sqrt(length(specificity)), uar.mu = mean(uar), uar.se = sd(uar)/sqrt(length(uar)))
  return(sum)
}

cRFs_trainfit <- sumPerfTrain(cRFs_trainfit)
rRFs_trainfit <- sumPerfTrain(rRFs_trainfit)

# Table S4/5
write.csv(cRFs_trainfit, "figures_tables/tableS4.csv")
write.csv(rRFs_trainfit, "figures_tables/tableS5.csv")

## CG SUMMARY
clean_covid <- read.delim("../data/clean_noCGC/covid-clean-0.txt")
clean_rsv <- read.delim("../data/clean_noCGC/rsv-clean-0.txt")

clean_covid$Data <- "COVID"
clean_rsv$Data <- "RSV"

covid <- read.delim("../data/processed/covid-cg.txt")
rsv <- read.delim("../data/processed/rsv-cg.txt")

# calculate the number of clonal groups (CG) observed at day 0
sumCG_CP <- ddply(clean_covid, c("Data", "PatientID"), summarise, totalCG = length(unique(CloneID)))
distCG_CP <- ddply(covid, c("PatientID"), summarise, meanSize = round(mean(nobs_d0)), sdCG = round(sd(nobs_d0)))
sumCG_CP <- merge(sumCG_CP, distCG_CP, by = "PatientID" )


sumCG_RSV <- ddply(clean_rsv, c("Data","PatientID"), summarise, totalCG =length(unique(CloneID)))
distCG_RSV <- ddply(rsv, c("PatientID"), summarise, meanSize = round(mean(nobs_d0)), sdCG = round(sd(nobs_d0)))
sumCG_RSV <- merge(sumCG_RSV, distCG_RSV, by = "PatientID" )


CG_sum <- rbind(sumCG_CP,sumCG_RSV)

# Table S1
write.csv(CG_sum, "figures_tables/tableS1.csv")

## MODEL INTERPRETATION PERFORMANCE
# Table 3
source("../scripts/RSV-CP_L1DO_LASSO_modInt.R")
# the two models below do take a while to run, so instead the results will be pulled from the CSV, please feel free to uncomment.
#source("../scripts/RSV_L1DO_RFx_modInt.R") 
#source("../scripts/CP_L1DO_RFx_modInt.R")

row1 <- read.csv("CP_modInt_L1DO_LASSO.csv")
row2 <- read.csv("CP_modInt_L1DO_RFx.csv")[2:12]
row3 <- read.csv("RSV_modInt_L1DO_LASSO.csv")
row4 <- read.csv("RSV_modInt_L1DO_RFx.csv")[2:12]


cols <- c('data', 'model', 'sensitivity', 'specificity', 'Precision', 'F1', 'uar', 'g.uar', 'AUC_ROC')
  
tab3 <- rbind(row1[cols],row2[cols],row3[cols],row4[cols])
tab3[cols[3:9]] <- round(tab3[cols[3:9]] * 100)

write.csv(tab3, "figures_tables/table3.csv")

######################################## FIGURE 2 ###################################################

## BOXPLOT
cols <- c("data", "model", "iteration", "sensitivity", "specificity", "uar")

RSV.AD <- list(R_LASSO.AD[cols],R_LR.AD[cols],R_RF.AD[cols]
               ,R_RFx.AD[cols]
               ,R_SVM.AD[cols]
) %>% bind_rows()
COVID.AD <- list(C_LASSO.AD[cols],C_LR.AD[cols],C_RF.AD[cols]
                 ,C_RFx.AD[cols]
                 ,C_SVM.AD[cols]
) %>% bind_rows()

RSV.DD <- list(R_LASSO.DD[cols],R_LR.DD[cols]
               ,R_RF.DD[cols]
               ,R_RFx.DD[cols]
               ,R_SVM.DD[cols]
) %>% bind_rows()
COVID.DD <- list(C_LASSO.DD[cols], C_LR.DD[cols]
                 ,C_RF.DD[cols]
                 , C_RFx.DD[cols]
                 ,C_SVM.DD[cols]
) %>% bind_rows()

# add row for test type and combine
RSV.AD$test <- "AD"
RSV.DD$test <- "L1DO"
COVID.AD$test <- "AD"
COVID.DD$test <- "L1DO"

RSV <- rbind(RSV.AD,RSV.DD)
COVID <- rbind(COVID.AD,COVID.DD)

plt.r <- ggplot(RSV, aes(fill = test)) +
  geom_boxplot(aes(x = model, y = uar)) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey")

plt.c <- ggplot(COVID, aes(fill = test)) +
  geom_boxplot(aes(x = model, y = uar)) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey")


uar.r <- plt.r + ylim(40,100) + labs(x = "" , y = "UAR (%)") + theme_bw() + scale_fill_brewer(palette = "Pastel1", name = "Experiment") + theme(axis.text = element_text(size=13)) + ylab("")

uar.c <- plt.c + ylim(40,100) + labs(x = "" , y = "UAR (%)") + theme_bw() + scale_fill_brewer(palette = "Pastel1", name = "Experiment") + theme(axis.text = element_text(size=13), legend.position = "none") 

plot_grid(uar.c,uar.r, nrow = 1,labels=c("A", "B"), rel_widths = c(1,1.25))

ggsave("figures_tables/figure2.jpeg", width = 10, height = 5)



######################################## FIGURES 3/4 ###################################################

## MODEL VARS
# extract variables
extractVariables_LASSO <- function(n_iter, Data, model_filename_root){
  all_coeffs <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("Variable",   "Estimate", "Data", "Iteration"))))
  if(n_iter == 1){
    modelname <- paste0(model_filename_root,".rda")
    model <- readRDS(modelname)
    
    coefs <- data.frame(Variable=rownames(model$beta), Estimate = matrix(model$beta) )
    coefs$Data <- Data
    coefs$Iteration <- 1
    
    all_coeffs <- rbind(all_coeffs, coefs)
  } else {
    for(i in 1:n_iter){
      modelname <- paste0(model_filename_root,i,".rda")
      model <- readRDS(modelname)
      coefs <- data.frame(Variable=rownames(model$beta), Estimate = matrix(model$beta) )
      
      coefs$Data <- Data
      coefs$Iteration <- i
      
      all_coeffs <- rbind(all_coeffs, coefs)
    }
    
    # remove ` in variable names
    all_coeffs$Variable <- gsub("`", "", all_coeffs$Variable)
    
    return(all_coeffs)
  }
}

cLASSO <- extractVariables_LASSO(1, "COVID", "COVID_modInt_LASSO")
rLASSO <- extractVariables_LASSO(1, "RSV", "RSV_modInt_LASSO")

# how many picked in LASSO
cLASSO <- cLASSO[cLASSO$Estimate != 0, ]
rLASSO <- rLASSO[rLASSO$Estimate != 0, ]

# order so most impact are above
cLASSO <- cLASSO[order(abs(cLASSO$Estimate),decreasing = TRUE),]

# order so most impactful are above
rLASSO <- rLASSO[order(abs(rLASSO$Estimate),decreasing = TRUE),]

# LASSO Analysis
cLASSO[,"Variable"] <- cLASSO[,"Variable"] %>% 
  str_replace("Vgene", "") %>%
  str_replace("Dgene", "") %>%
  str_replace("Jfamily", "")

rLASSO[,"Variable"] <- rLASSO[,"Variable"] %>% 
  str_replace("Vgene", "") %>%
  str_replace("Dgene", "") %>%
  str_replace("Dfamily", "") %>%
  str_replace("Jfamily", "")

# add variable type
varLab <- function(varVec){
  # CG feature average
  varVec$type[grepl(".mu",varVec$Variable)] <- "CG feature average" 
  
  # CG feature variation
  varVec$type[grepl(".sd",varVec$Variable)] <- "CG feature variation" 
  
  # CG descriptor
  varVec$type[grepl("IGH",varVec$Variable)] <- "CG descriptor" 
  varVec$type[grepl("nobs_d0",varVec$Variable)] <- "CG descriptor" 
  
  # metadata
  varVec$type[varVec$Variable == "Age" | varVec$Variable == "Sex" | varVec$Variable == "SexMale" | varVec$Variable == "AgeYoung"] <- "Metadata"
  
  # Baseline
  varVec$type[grepl("Intercept",varVec$Variable)] <- "Baseline" 
  
  varVec$type <- factor(varVec$type, levels = c("CG feature variation", "CG descriptor","CG feature average", "Metadata", "Baseline"))
  
  return(varVec)
}

cLASSO <- varLab(cLASSO)
rLASSO <- varLab(rLASSO)

library(viridis)
library(wesanderson)
library(ggsci)
library(ggbreak)

## LASSO Analysis
plt.c <- ggplot(cLASSO[1:20,], aes(Estimate, reorder(Variable, +Estimate), fill = type)) +
  geom_col(position = "dodge") + 
  theme_bw() + 
  ylab("") +
  xlab("Coefficient") + 
  scale_fill_npg()  +
  theme(legend.title = element_blank(), text = element_text(size=14), legend.key.size = unit(.5,'cm'), legend.spacing.x = unit(.5, 'cm'))
# extract a legend that is laid out horizontally
legend <- get_legend(plt.c)
# remove legend
plt.c <- plt.c + theme(legend.position = "'none")

plt.r <- ggplot(rLASSO[1:20,], aes(Estimate, reorder(Variable, +Estimate), fill = type)) +
  geom_col(position = "dodge") +
  ylab("") +
  xlab("Coefficient") + 
  theme_bw() + 
  theme(legend.title = element_blank(), text = element_text(size=14),legend.position = "none") +
  scale_fill_npg()

plot_grid(plt.c,plt.r, legend, nrow = 1,labels=c("A", "B"), rel_widths = c(1,1,.45))
ggsave("figures_tables/figure3.jpeg", width = 13, height = 5)


## RFx Analysis
# load models
cmod <- readRDS("COVID_modInt_RFx.rda")
rmod <- readRDS("RSV_modInt_RFx.rda")

library(randomForest)

cRFx <- data.frame(name = rownames(importance(cmod)), MeanDecreaseGini=importance(cmod))
rRFx <- data.frame(name = rownames(importance(rmod)), MeanDecreaseGini=importance(rmod))

# order so most impact are above
cRFx <- cRFx[order(abs(cRFx$MeanDecreaseGini),decreasing = TRUE),]

# order so most impactful are above
rRFx <- rRFx[order(abs(rRFx$MeanDecreaseGini),decreasing = TRUE),]

# add variable type
varLab <- function(varVec){
  # CG feature average
  varVec$type[grepl(".mu",varVec$name)] <- "CG feature average" 
  
  # CG feature variation
  varVec$type[grepl(".sd",varVec$name)] <- "CG feature variation" 
  
  # CG descriptor
  varVec$type[grepl("gene",varVec$name)] <- "CG descriptor" 
  varVec$type[grepl("family",varVec$name)] <- "CG descriptor" 
  varVec$type[grepl("nobs_d0",varVec$name)] <- "CG descriptor" 
  
  # metadata
  varVec$type[varVec$name == "Age" | varVec$name == "Sex" | varVec$name == "SexMale" | varVec$name == "AgeYoung"] <- "Metadata"
  
  # Baseline
  varVec$type[grepl("Intercept",varVec$name)] <- "Baseline" 
  
  varVec$type <- factor(varVec$type, levels = c("CG feature variation", "CG descriptor","CG feature average", "Metadata", "Baseline"))
  
  return(varVec)
}

cRFx <- varLab(cRFx)
rRFx <- varLab(rRFx)

plt.c <- ggplot(cRFx[1:20,], aes(MeanDecreaseGini, reorder(name, +MeanDecreaseGini), fill = type)) +
  geom_col(position = "dodge") + 
  theme_bw() + 
  ylab("") +
  xlab("Mean Decrease Gini") + 
  #scale_fill_brewer(palette = "Pastel1") + 
  scale_fill_npg() +
  theme(legend.title = element_blank(), text = element_text(size=14), legend.key.size = unit(.5,'cm'), legend.spacing.x = unit(.5, 'cm'))

# extract a legend that is laid out horizontally
legend <- get_legend(plt.c)
# remove legend
plt.c <- plt.c + theme(legend.position = "'none")

plt.r <- ggplot(rRFx[1:20,], aes(MeanDecreaseGini, reorder(name, +MeanDecreaseGini), fill = type)) +
  geom_col(position = "dodge") +
  ylab("") +
  xlab("Mean Decrease Gini") + 
  theme_bw() + 
  theme(legend.title = element_blank(), text = element_text(size=14),legend.position = "none") +
  scale_fill_npg()

plot_grid(plt.c,plt.r, legend, nrow = 1,labels=c("A", "B"), rel_widths = c(1,1,.57))
ggsave("figures_tables/figure4.jpeg", width = 10, height = 5)



######################################## FIGURES S1/S2 ###################################################

# run LASSO
source("../scripts/RSV-CP_LASSO.R")
source("../scripts/RSV-CP_L1DO_LASSO.R")


cEst <- extractVariables_LASSO(5, "COVID", "cLASSO_model")
rEst <- extractVariables_LASSO(5, "RSV", "rLASSO_model")
cEst_L1DO <- extractVariables_LASSO(16, "COVID", "cLASSO_L1DO_model")
rEst_L1DO <- extractVariables_LASSO(12, "RSV", "rLASSO_L1DO_model")

# most consistently picked?
cEst.freq <- cEst[cEst$Estimate != 0,] %>% count(Variable)
cEst_L1DO.freq <- cEst_L1DO[cEst_L1DO$Estimate != 0,] %>% count(Variable)
rEst.freq <- rEst[rEst$Estimate != 0,] %>% count(Variable)
rEst_L1DO.freq <- rEst_L1DO[rEst_L1DO$Estimate != 0,] %>% count(Variable)

# only select the variables in all loops
cEst.top <- cEst.freq %>% 
  arrange(desc(n)) %>% 
  slice(1:32) %>%
  select(-c(n)) 

cEst_L1DO.top <- cEst_L1DO.freq %>% 
  arrange(desc(n)) %>% 
  slice(1:15) %>% 
  select(-c(n)) 

rEst.top <- rEst.freq %>% 
  arrange(desc(n)) %>% 
  slice(1:21) %>%
  select(-c(n))

rEst_L1DO.top <- rEst_L1DO.freq %>% 
  arrange(desc(n)) %>% 
  slice(1:38) %>% 
  select(-c(n))


# find the intersection of those always picked for both AD and L1DO
cInt <- intersect(cEst.top,cEst_L1DO.top) 

cInt_var <- intersect(cEst.top,cEst_L1DO.top)$Variable %>% 
  str_replace("Vgene.+", "Vgene") %>%
  str_replace("Dgene.+", "Dgene") %>%
  str_replace("Jfamily.+", "Jfamily") %>%
  unique() 

rInt <- intersect(rEst.top,rEst_L1DO.top) 


rInt_var <- intersect(rEst.top,rEst_L1DO.top)$Variable %>% 
  str_replace("Vgene.+", "Vgene") %>%
  str_replace("Dgene.+", "Dgene") %>%
  str_replace("Dfamily", "") %>%
  str_replace("Jfamily.+", "Jfamily") %>%
  unique() 


# combine estimates
cEst$Experiment <- "AD"
cEst_L1DO$Experiment <- "L1DO"
cEst.bind <- rbind(cEst[cEst$Variable %in% cInt$Variable,],cEst_L1DO[cEst_L1DO$Variable %in% cInt$Variable,])

rEst$Experiment <- "AD"
rEst_L1DO$Experiment <- "L1DO"
rEst.bind <- rbind(rEst[rEst$Variable %in% rInt$Variable,],rEst_L1DO[rEst_L1DO$Variable %in% rInt$Variable,])


# clean naming of variables
cEst.bind[,"Variable"] <- cEst.bind[,"Variable"] %>% 
  str_replace("Vgene", "") %>%
  str_replace("Dgene", "") %>%
  str_replace("Dfamily", "") %>%
  str_replace("Jfamily", "")

rEst.bind[,"Variable"] <- rEst.bind[,"Variable"] %>% 
  str_replace("Vgene", "") %>%
  str_replace("Dgene", "") %>%
  str_replace("Dfamily", "") %>%
  str_replace("Jfamily", "")

coeffs <- rbind(cEst.bind,rEst.bind)

coeff.dist <- ddply(coeffs, c("Data", "Experiment","Variable"), summarise, est.median = median(Estimate), est.iqr = IQR(Estimate))

plt.c <- ggplot(coeff.dist[coeff.dist$Data == "COVID",], aes(est.median, Variable, fill = Experiment)) +
  geom_col(position = "dodge") + 
  #geom_errorbar(aes(ymin = est.median - est.iqr, ymax = est.median + est.iqr), width=0, position=position_dodge(.9)) +
  theme_bw() + 
  ylab("Estimate") + 
  theme(axis.title.x = element_blank(), legend.position = "none") +
  scale_x_continuous(breaks = -1:4) +
  #scale_fill_brewer(palette = "Pastel1") + 
  scale_fill_npg()

plt.r <- ggplot(coeff.dist[coeff.dist$Data == "RSV",], aes(est.median, Variable, fill = Experiment)) +
  geom_col(position = "dodge") +
  #geom_errorbar(aes(ymin = est.median - est.iqr, ymax = est.median + est.iqr), width=0, position=position_dodge(.9)) +
  ylab("") +
  theme_bw() + 
  theme(axis.title.x = element_blank(), legend.title = element_blank()) +
  scale_x_continuous(breaks = -1:4) +
  #scale_fill_brewer(palette = "Pastel1") + 
  scale_fill_npg()

plot_grid(plt.c,plt.r, nrow = 1,labels=c("A", "B"), rel_widths = c(1,1.2))
ggsave("figures_tables/figureS1.jpeg", width = 10, height = 5)



cHP <- read.csv("../reduced-scripts/hyperparam/CP_hyPar_RF-NSx.csv")
cHP_L1DO <- read.csv("../reduced-scripts/hyperparam/CP_L1DO_hyPar_RF-NSx.csv")
rHP <- read.csv("../reduced-scripts/hyperparam/RSV_hyPar_RF-NSx.csv")
rHP_L1DO <- read.csv("../reduced-scripts/hyperparam/RSV_L1DO_hyPar_RF-NSx.csv")


# Create formula to calc. importance as scripts

# balancing data
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


extractRFimp <- function(data,fold,nfolds,hyPar,imp_vec){
  for(i in 1:nfolds){
    # create list of training sets
    train <- do.call(rbind, lapply(fold[-i], function(y){
      tr <- data[y,]
      return(tr)
    }))
    
    hp_f <- hyPar[i,]
    train <- balance.50(train)
    
    # training the model
    rf <- randomForest(IsCrossClass ~ .,
                       mtry = hp_f$m, 
                       nodesize = hp_f$min,
                       ntree = hp_f$N, 
                       data = train)
    
    # variable importance
    imp <- importance(rf) %>% as.data.frame() %>% rownames_to_column()
    imp_vec <- rbind(imp_vec, imp)
  }
  
  return(imp_vec)
}


extractRFimp_L1DO <- function(data,nfolds,hyPar_L1DO,imp_L1DO_vec){
  for(i in 1:nfolds){
    # create list of training sets
    train <- do.call(rbind, lapply(data[-i], function(y){
      tr <- y
      return(tr)
    }))
    
    hp_f <- hyPar_L1DO[i,]
    
    train <- balance.50(train)
    
    # training the model
    rf <- randomForest(IsCrossClass ~ .,
                       mtry = hp_f$m, 
                       nodesize = hp_f$min,
                       ntree = hp_f$N, 
                       data = train)
    
    # variable importance
    imp_L1DO <- importance(rf) %>% as.data.frame() %>% rownames_to_column()
    imp_L1DO_vec <- rbind(imp_L1DO_vec, imp_L1DO)
  }
  
  return(imp_L1DO_vec)
}


## COVID Analysis

set.seed(346)
###### PROCESS DATA ###### 
covid <- read.delim("../data/processed/covid-cg.txt")
# remove pid
covid <- select(covid, -c(PatientID,CloneID))
###### MODEL SET-UP ###### 
# make characters into factors
cat <- c("IsCrossClass", "Vgene", "Jfamily", "Dfamily", "Dgene", "Sex")
covid[cat] <- lapply(covid[cat], factor)

# create folds
n_k <- 5 
fold <- createFolds(covid$IsCrossClass,k=n_k)

cRFx.imp <- data.frame()
cRFx.imp <- extractRFimp(covid,fold,n_k,cHP,cRFx.imp)

# Save out result
write.csv(cRFx.imp, file = "cRFx_AD_imp_vec.csv")


set.seed(346)
###### PROCESS DATA ###### 
covid <- read.delim("../data/processed/covid-cg.txt")
# remove pid
covid <- select(covid, -c(CloneID))

###### MODEL SET-UP ###### 
# make characters into factors
cat <- c("IsCrossClass", "Vgene", "Jfamily", "Dfamily", "Dgene", "Sex")
covid[cat] <- lapply(covid[cat], factor)

# patient separation
pid <- unique(covid$PatientID)
covid <- lapply(pid, function(x){
  p <- subset(covid, PatientID == x , select = -c(PatientID))
})

# patient combinations is now (d-1)vs1 
n_s <- length(pid)  

cRFx_L1DO.imp <- data.frame()
cRFx_L1DO.imp <- extractRFimp_L1DO(covid,n_s,cHP_L1DO,cRFx_L1DO.imp)

# Save out result
write.csv(cRFx_L1DO.imp, file = "cRFx_L1DO_imp_vec.csv")


# Combine Results
cRFx.imp$Experiment <- "AD"
cRFx_L1DO.imp$Experiment <- "L1DO"

cRFx.bind <- rbind(cRFx.imp,cRFx_L1DO.imp)
names(cRFx.bind)[1] <- "Variable"

cRFx.bind.med <- ddply(cRFx.bind, c("Experiment", "Variable"), summarise, iqr = IQR(MeanDecreaseGini), MeanDecreaseGini = median(MeanDecreaseGini))



## RSV Analysis
set.seed(346)
###### PROCESS DATA ###### 
rsv <- read.delim("../data/processed/rsv-cg.txt")
# remove pid
rsv <- select(rsv, -c(PatientID,CloneID))

print("data processed")

###### MODEL SET-UP ###### 
# make characters into factors
cat <- c("IsCrossClass", "Vgene", "Jfamily", "Dfamily", "Dgene", "Age", "Sex")
rsv[cat] <- lapply(rsv[cat], factor)

# create folds
n_k <- 5 
fold <- createFolds(rsv$IsCrossClass,k=n_k)

rRFx.imp <- data.frame()
rRFx.imp <- extractRFimp(rsv,fold,n_k,rHP,rRFx.imp)

# Save out result
write.csv(rRFx.imp, file = "rRFx_AD_imp_vec.csv")

set.seed(346)
###### PROCESS DATA ###### 
rsv <- read.delim("../data/processed/rsv-cg.txt")
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

rRFx_L1DO.imp <- data.frame()
rRFx_L1DO.imp <- extractRFimp_L1DO(rsv,n_s,rHP_L1DO,rRFx_L1DO.imp)

# Save out result
write.csv(rRFx_L1DO.imp, file = "rRFx_L1DO_imp_vec.csv")

# Combine results
rRFx.imp$Experiment <- "AD"
rRFx_L1DO.imp$Experiment <- "L1DO"

rRFx.bind <- rbind(rRFx.imp,rRFx_L1DO.imp)
names(rRFx.bind)[1] <- "Variable"

rRFx.bind.med <- ddply(rRFx.bind, c("Experiment", "Variable"), summarise, iqr = IQR(MeanDecreaseGini), MeanDecreaseGini = median(MeanDecreaseGini))


# union of top 10 median mean decrease in Gini between CV strategies
cAD.top10 <- cRFx.bind.med[cRFx.bind.med$Experiment == "AD",] %>% slice_max(MeanDecreaseGini, n = 10)
rAD.top10 <- rRFx.bind.med[rRFx.bind.med$Experiment == "AD",] %>% slice_max(MeanDecreaseGini, n = 10)

cL1DO.top10 <- cRFx.bind.med[cRFx.bind.med$Experiment == "L1DO",] %>% slice_max(MeanDecreaseGini, n = 10)
rL1DO.top10 <- rRFx.bind.med[rRFx.bind.med$Experiment == "L1DO",] %>% slice_max(MeanDecreaseGini, n = 10) 

ctop10.union <- cRFx.bind.med[cRFx.bind.med$Variable %in% union(cAD.top10$Variable,cL1DO.top10$Variable),]
rtop10.union <- rRFx.bind.med[rRFx.bind.med$Variable %in% union(rAD.top10$Variable,rL1DO.top10$Variable),]

# add variable type
varLab <- function(varVec){
  # CG feature average
  varVec$type[grepl(".mu",varVec$Variable)] <- "CG feature average" 
  
  # CG feature variation
  varVec$type[grepl(".sd",varVec$Variable)] <- "CG feature variation" 
  
  # CG descriptor
  varVec$type[grepl("gene",varVec$Variable)] <- "CG descriptor" 
  varVec$type[grepl("nobs_d0",varVec$Variable)] <- "CG descriptor" 
  
  # metadata
  varVec$type[varVec$Variable == "Age" | varVec$Variable == "Sex"] <- "metadata"
  
  varVec$type <- factor(varVec$type, levels = c("CG feature variation", "CG descriptor","CG feature average", "metadata"))
  
  return(varVec)
}

ctop10.union <- varLab(ctop10.union)
rtop10.union <- varLab(rtop10.union)

print(ctop10.union)
print(rtop10.union)

plt.c <- ggplot(ctop10.union, aes(MeanDecreaseGini,Variable, fill = type)) +
  geom_col() + 
  geom_errorbarh(aes(xmin = MeanDecreaseGini - iqr, xmax = MeanDecreaseGini + iqr), height=0) +
  facet_grid(~Experiment, scales = "free_x") +
  theme_bw() + 
  scale_fill_npg() +
  #scale_fill_npg() +
  #scale_fill_manual(values = wes_palette("Rushmore", n = 4)) +
  #scale_fill_viridis(discrete = TRUE) +
  #scale_fill_brewer(palette = "Set2") +
  scale_y_discrete(limits=ctop10.union[order(ctop10.union$MeanDecreaseGini[ctop10.union$Experiment == "AD"]),"Variable"]) +
  ylab("") +
  theme(legend.title = element_blank(),axis.title = element_text(size=10), axis.title.x = element_blank(), legend.position = "bottom", legend.text = element_text(size=10), legend.key.size = unit(.5,'cm'), legend.spacing.x = unit(.5, 'cm'))

# extract a legend that is laid out horizontally
legend_b <- get_legend(plt.c)
# remove legend
plt.c <- plt.c + theme(legend.position = "'none")

plt.r <- ggplot(rtop10.union, aes(MeanDecreaseGini,Variable, fill = type)) +
  geom_col() + 
  geom_errorbarh(aes(xmin = MeanDecreaseGini - iqr, xmax = MeanDecreaseGini + iqr), height=0) +
  facet_grid(~Experiment, scales = "free_x") +
  theme_bw() + 
  scale_fill_npg() +
  #scale_fill_npg() +
  #scale_fill_manual(values = wes_palette("Rushmore", n = 4)) +
  #scale_fill_viridis(discrete = TRUE) +
  #scale_fill_brewer(palette = "Set2") +
  scale_y_discrete(limits=rtop10.union[order(rtop10.union$MeanDecreaseGini[rtop10.union$Experiment == "AD"]),"Variable"]) +
  xlab("Mean Decrease Gini") +
  ylab("") +
  theme(legend.position = "none") 

plot_grid(plt.c,plt.r,legend_b, nrow = 3,labels=c("A", "B"), rel_heights = c(1,1,.1))
ggsave("figures_tables/figureS2.jpeg" , width = 8, height = 8)
