# filter annotated antibody repertoire to only have heavy chain observations 
heavy.chain <- function(data){
  s = subset(data, (Class != "K" & Class != "L" & Subclass != "IgA" & Subclass != "IgG" & grepl("IGH", Vfamily)))
  return(s)
}


# artificially re-balance data for classifier: oversampling rare class
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

# subset of features in order to clean the dataset by removing duplicate sequences
sub.feature <- function(
    df,
    classifier = TRUE,
    clonal_group = FALSE,
    clone_no = c("CloneID","Vgene","Jfamily","Dfamily","Dgene"),
    response = "IsCrossClass",
    unique,
    patient = "PatientID",
    features
){
  print("Subsetting for functionality ...") # outputting progress
  
  df <- subset(df, V.DOMAIN.Functionality == "productive") # filter productive cells
  
  print("Grouping ...")
  
  df <- ddply(df, features, summarise, n = length(Seq_ID)) # grouping the sequences based on the features
  
  print("Splitting ...")
  
  pid <- unique(df[[patient]])
  
  clean_split <- do.call("rbind", lapply(pid, function(p){
    
    cat("patient", p, "...")
    
    df_p <- subset(df, PatientID = p) # subset for singular patient/donor
    
    df_seq <- unique(df_p[[unique]]) # make list of unique sequences in the dataset
    
    clean_split_p <- lapply(df_seq, function(s){ # split table by the cell column to select the most representative 
      subset(df_p, Sequence == s)
    })
    
    clean_split_p <- do.call("rbind", lapply(clean_split_p, function(x){
      x <- x[order(x$n, decreasing = TRUE),]
      x[1,]
    }))
    
    return(clean_split_p)
  }))
  
  print("Removing NAs...")
  
  for(i in 1:length(features)){
    clean_split <- clean_split[!is.na(clean_split[features[i]]),] # remove NA rows for features
  }
  
  clean <- clean_split[!(clean_split$Dfamily == "None"),] # remove "None" from Dfamily
  
  
  if (clonal_group == TRUE){
    print("Selecting most abundant genes for clonal group representation...")
    
    
    clean <- do.call("rbind", lapply(pid, function(y){  # do the following on a per patient basis
      
      cat("patient", y, "...")
      
      data <- subset(clean_split, PatientID == y)
      
      cg_abundance <- ddply(data, clone_no, summarise, n = length(CloneID))
      
      all_cg <- unique(data$CloneID)   # get all cloneIDs
      
      cg_list <- lapply(all_cg, function(x){  # order observations by clonal group and variable region categories. select most abundant
        ls <- subset(cg_abundance, CloneID == x)
        ls <- ls[order(ls$n, decreasing = TRUE),]
        
        return(ls[1,])
      })
      
      cg_unique <- do.call("rbind", lapply(cg_list, function(x){
        i1 <- x[[clone_no[1]]]    # values of most abundant
        i2 <- x[[clone_no[2]]] 
        i3 <- x[[clone_no[3]]] 
        i4 <- x[[clone_no[4]]] 
        i5 <- x[[clone_no[5]]] 
        
        result <- data[(data[[clone_no[1]]]  == i1) & 
                         (data[[clone_no[2]]]  == i2) & 
                         (data[[clone_no[3]]]  == i3) & 
                         (data[[clone_no[4]]]  == i4) & 
                         (data[[clone_no[5]]]  == i5),]
        
        return(result)
      }))
      
    }))
  }
  
  if (classifier == TRUE){
    clean[[response]] <- as.factor(clean[[response]]) # make response variable a factor
  }
  
  clean <- select(clean, -c("n", all_of(unique))) # remove count and sequence
  return(clean)
  
}

