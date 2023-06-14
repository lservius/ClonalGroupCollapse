library(tidyverse)
library(plyr)
library(doParallel)
setwd("/scratch/users/k21111073/")
source("functions_dataProcessing.R")

set.seed(346)

### IMPORT ###

mix <- read.delim("data/AntibodyRepertoireCV+RSV+EBOV+YFV.txt")


### SPLIT DATA ###

names(mix)[names(mix) == "Gender"] <- "Sex"

# list of challenges
ch <- unique(mix$SampleType)

# split by challenge
split <- lapply(ch, function(x){
  data <- subset(mix, SampleType == x)
})

covid <- split[[1]]
rsv <- rbind(split[[3]],split[[4]])

# all together now
ch <- c("covid","rsv")
data <- list(covid,rsv)
names(data) <- ch

# limit to day 0 for model comparison sake
data <- lapply(data,function(x){
  day0 <- subset(x, TimePoint == 0)
})

### PROCESS DATA ###

# filter out LC observations and undetermined subclass
data <- lapply(data, heavy.chain)

# select features to use for subset                                               
f <- c("IsCrossClass","PatientID","TimePoint","Class","Subclass","Vgene","Jfamily","Dfamily","Dgene","CloneID","Seq_ID","NumInClone",
       "V.REGION.identity..","J.REGION.identity..","Sequence","TotalN","TotalP","TotalD","Num_AAs","Tiny_ACGST",
       "Small_ABCDGNPSTV","Aliphatic_AILV",
       "Aromatic_FHWY","NonPolar_ACFGILMPVWY","Polar_DEHKNQRSTZ","Charged_BDEHKRZ","Basic_HKR","Acidic_BDEZ","Aliphatic_Index",
       "Boman","pI_EMBOSS","Hydrophobicity","Instability","kidera1.helix.bend.pref","kidera2.side.chain.size",
       "kidera3.extended.str.pref","kidera4.hydrophobicity","kidera5.double.bend.pref","kidera6.partial.spec.vol",
       "kidera7.flat.ext.pref","kidera8.occurrence.alpha.reg","kidera9.pK.C","kidera10.surrounding.hydrop","Age","Sex")

# subset and clean data
numCores <- 2
cl <- makeCluster(numCores)

clusterExport(cl, c("data", "f"))

clusterEvalQ(cl, {
      library(tidyverse)
      library(plyr)
})

data <- parLapply(cl, data, sub.feature,
                  features = f, 
                  unique = "Sequence", 
                  clonal_group = TRUE)

stopCluster(cl)

# make characters into factors
cat <- c("Class","Subclass","Vgene","Jfamily","Dfamily", "Sex","Age")
data <- lapply(data, function(x){
  x[cat] <- lapply(x[cat], factor)
  return(x)
})

# remove NumInClone == 1
data <- lapply(data, function(x){subset(x, NumInClone>1)})

save.image("clonalGroupCollapse")

### COLLAPSE CLONE ###

# collapse based on CloneID
col <- c("IsCrossClass","PatientID","Vgene","Jfamily","Dfamily","Dgene","CloneID","NumInClone","Age","Sex")

data.cg <- lapply(data, function(x){
  ddply(x, col, summarise, 
        # average the continuous variables over the clonal groups
        
        # mean of values
        V.REGION.id.mu = mean(V.REGION.identity..),
        J.REGION.id.mu = mean(J.REGION.identity..),
        TotalN.mu = mean(TotalN),
        TotalP.mu = mean(TotalP),
        TotalD.mu = mean(TotalD),
        Num_AAs.mu = mean(Num_AAs),
        Tiny.mu = mean(Tiny_ACGST),
        Small.mu = mean(Small_ABCDGNPSTV),
        Aromatic.mu = mean(Aromatic_FHWY),
        NonPolar.mu = mean(NonPolar_ACFGILMPVWY),
        Polar.mu = mean(Polar_DEHKNQRSTZ),
        Charged.mu = mean(Charged_BDEHKRZ),
        Basic.mu = mean(Basic_HKR),
        Acidic.mu = mean(Acidic_BDEZ),
        Aliphatic.mu = mean(Aliphatic_Index),
        Boman.mu = mean(Boman),
        pI_EMBOSS.mu = mean(pI_EMBOSS),
        Hydrophobicity.mu = mean(Hydrophobicity),
        Instability.mu = mean(Instability),
        kidera1.mu = mean(kidera1.helix.bend.pref),
        kidera2.mu = mean(kidera2.side.chain.size),
        kidera3.mu = mean(kidera3.extended.str.pref),
        kidera4.mu = mean(kidera4.hydrophobicity),
        kidera5.mu = mean(kidera5.double.bend.pref),
        kidera6.mu = mean(kidera6.partial.spec.vol),
        kidera7.mu = mean(kidera7.flat.ext.pref),
        kidera8.mu = mean(kidera8.occurrence.alpha.reg),
        kidera9.mu = mean(kidera9.pK.C),
        kidera10.mu = mean(kidera10.surrounding.hydrop),
        
        # sd of values
        V.REGION.id.sd = sd(V.REGION.identity..),
        J.REGION.id.sd = sd(J.REGION.identity..),
        TotalN.sd = sd(TotalN),
        TotalP.sd = sd(TotalP),
        TotalD.sd = sd(TotalD),
        Num_AAs.sd = sd(Num_AAs),
        Tiny.sd = sd(Tiny_ACGST),
        Small.sd = sd(Small_ABCDGNPSTV),
        Aromatic.sd = sd(Aromatic_FHWY),
        NonPolar.sd = sd(NonPolar_ACFGILMPVWY),
        Polar.sd = sd(Polar_DEHKNQRSTZ),
        Charged.sd = sd(Charged_BDEHKRZ),
        Basic.sd = sd(Basic_HKR),
        Acidic.sd = sd(Acidic_BDEZ),
        Aliphatic.sd = sd(Aliphatic_Index),
        Boman.sd = sd(Boman),
        pI_EMBOSS.sd = sd(pI_EMBOSS),
        Hydrophobicity.sd = sd(Hydrophobicity),
        Instability.sd = sd(Instability),
        kidera1.sd = sd(kidera1.helix.bend.pref),
        kidera2.sd = sd(kidera2.side.chain.size),
        kidera3.sd = sd(kidera3.extended.str.pref),
        kidera4.sd = sd(kidera4.hydrophobicity),
        kidera5.sd = sd(kidera5.double.bend.pref),
        kidera6.sd = sd(kidera6.partial.spec.vol),
        kidera7.sd = sd(kidera7.flat.ext.pref),
        kidera8.sd = sd(kidera8.occurrence.alpha.reg),
        kidera9.sd = sd(kidera9.pK.C),
        kidera10.sd = sd(kidera10.surrounding.hydrop),
        
        # see how many are averaged for each entry for future duplicate processing
        n = length(PatientID),

        # new numinclone, count unique cell IDs for each clone summary
        nobs_d0 = length(unique(Seq_ID))
  )})

### CHECK DUPLICATES ###
dup1 <- lapply(data.cg, function(x){
  n_occur <- ddply(x, c("PatientID", "CloneID"), summarise, freq = length(CloneID))
  n_occur[n_occur$freq > 1,]
})
dup1


# remove count n from data frame
data.clean <- lapply(data.cg, function(x){
  select(x, -c(n))
})


# check for duplicate CloneIDs
dup2 <- lapply(data.clean, function(x){
  n_occur <- ddply(x, c("PatientID", "CloneID"), summarise, freq = length(CloneID))
  n_occur[n_occur$freq > 1,]
})
dup2

# compare NumInClone with nobs_d0
proportion <- lapply(data.clean, function(x){
ddply(x, c("PatientID","CloneID"), summarise, prop = nobs_d0/NumInClone)
})
proportion

prop_check <- lapply(proportion, function(x){
  x[x$prop > 1,]
})
prop_check


# remove NAs for the sd where the number of observations was only 1
data.clean <- lapply(data.clean, function(x){
  replace(x,is.na(x),0)
})


### SAVE ###
write.table(data.clean$covid, file = "covid-cg.txt", sep = "\t", row.names = FALSE)
write.table(data.clean$rsv, file = "rsv-cg.txt", sep = "\t", row.names = FALSE)
