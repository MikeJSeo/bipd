
#' Find performance metrics
#'
#' Find performance metrics (MSE, MAE, R-squared) taking into account the study level.
#'
#' @param testoutcome Outcome of the testing data to calculate performance metrics
#' @param predictions Predictions from model developed
#' @param aggregation Aggregation method to be used to combine performance metrics across studies. 
#' There are three options: "simple" gives you a simple average across studies; 
#' "weighted" gives you a weighted average taking into account the sample size of each study. 
#' "ignore" gives you a performance measure ignoring study level. Default is set to "weighted".
#' @return 
#' \item{performancemetrics}{A vector of MSE, MAE, and R-squared values}
#' @export

findPerformance <- function(testoutcome = NULL, predictions = NULL, aggregation = "weighted"){
  
  if(is.null(testoutcome) | is.null(predictions)){
    stop("testoutcome or predictions is not specified")
  }
  
  if(!aggregation %in% c("weighted", "simple", "ignore")){
    stop("aggregation should be either weighted, simple or ignore")
  }
  
  if(aggregation == "ignore"){
    
    testoutcome_unlisted <- unlist(testoutcome)
    predictions_unlisted <- unlist(predictions)
    
    performancemetrics <- rep(NA, 3)
    performancemetrics[1] <- findMSE(testoutcome_unlisted, predictions_unlisted)
    performancemetrics[2] <- findMAE(testoutcome_unlisted, predictions_unlisted)
    performancemetrics[3] <- findRsquared(testoutcome_unlisted, predictions_unlisted)
    
    names(performancemetrics) <- c("MSE", "MAE","Rsquared")
    
  } else if(aggregation %in% c("weighted", "simple")){
    
    performances <- matrix(NA, nrow = 3, ncol = length(predictions))
    samplesize <- rep(NA, length(predictions))
    
    for(ii in 1:length(predictions)){
      
      performances[1, ii] <- findMSE(testoutcome[[ii]], predictions[[ii]])
      performances[2, ii] <- findMAE(testoutcome[[ii]], predictions[[ii]])
      performances[3, ii] <- findRsquared(testoutcome[[ii]], predictions[[ii]])
      
      samplesize[ii] <- length(testoutcome[[ii]])
    }
    
    if(aggregation == "weighted"){
      product_store <- performances * samplesize
      samplesize_sum <- sum(samplesize)
      
      final_store <- product_store/samplesize_sum
      performancemetrics <- apply(final_store, 1, sum)  
    } else if(aggregation == "simple"){
      performancemetrics <- apply(performances, 1, mean)
    }
    names(performancemetrics) <- c("MSE", "MAE", "Rsquared")
  }
  
  return(performancemetrics)
}


findRsquared <- function(y, pred){
  total <- (y - mean(y, na.rm = TRUE))^2
  tss <- sum(total[!is.na(total)])
  residual <- (y - pred)^2
  rss <- sum(residual[!is.na(residual)])
  rsquared <- 1 - rss/tss
  rsquared
}


findMAE <- function(y, pred){
  err_MAE <- abs(pred - y)
  err_MAE <- err_MAE[!is.na(err_MAE)]
  mean(err_MAE)
}


findMSE <- function(y, pred){
  err_mse <- (pred-y)^2
  err_mse <- err_mse[!is.na(err_mse)]
  mean(err_mse)
}
