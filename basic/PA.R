PreAccuracy <- function(data, theta){
  X <- data[,-1]
  Y <- data[, 1]
    pred <- (2*Y-1)*(X%*%theta)
    PA <- sum(pred>0)/nrow(data)
  return(PA)
}