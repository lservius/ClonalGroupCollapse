# convert probility output into a binary with a threshold of 0.5
predict.TF <- function(prediction, threshold = 0.5){
  predict_TF <- NULL
  for (i in 1:length(prediction)){
    if (is.na(prediction[i]) == TRUE){
      predict_TF[i] <- NA
    } else if (prediction[i] > threshold){
      predict_TF[i] <- TRUE
    } else {
      predict_TF[i] <-  FALSE
    }
  }
  return(predict_TF)
}