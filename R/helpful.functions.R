
#' Calculate patient-specific treatment effect
#'
#' Convenient function for calculating the patient-specific treatment effect.
#' Patient-specific treatment effect includes the main effect of treatment and 
#' treatment-covariate interaction effect (i.e. effect modification). 
#' Reports odds ratio for the binary outcome.
#' 
#' @param ipd IPD object created from running model functions i.e. `ipdma.onestage.model`
#' @param samples MCMC samples found from running `ipd.run`
#' @param newpatient Covariate values of patients that you want to predict treatment effect on. Must have length equal to total number of covariates.
#' @param reference Reference group used for finding patient-specific treatment effect. This is only used for "deft" approach
#' @param quantiles Quantiles finding credible interval of the patient-specific treatment effect
#' @references Seo M, White IR, Furukawa TA, et al. Comparing methods for estimating patient-specific treatment effects in individual patient data meta-analysis. \emph{Stat Med}. 2021;40(6):1553-1573. \doi{10.1002/sim.8859}
#' @references Riley RD, Debray TP, Fisher D, et al. Individual participant data meta-analysis to examine interactions between treatment effect and participant-level covariates: Statistical recommendations for conduct and planning. \emph{Stat Med}. 2020:39(15):2115-2137. [\url{https://doi.org/10.1002/sim.8516}] 
#' @export

treatment.effect <- function(ipd = NULL, samples = NULL, newpatient = NULL, 
                             reference = NULL, quantiles = c(0.025, 0.5, 0.975)){

  
  if(class(ipd) == "ipdma.onestage"){
    
    newpatient <- (newpatient - ipd$scale_mean)/ipd$scale_sd
      
    index0 <- which(colnames(samples[[1]]) == "delta[2]") 
    index1 <- grep("gamma", colnames(samples[[1]]))
    index <- c(index0, index1)
    samples2 <- samples[,index]
    
    merged <- samples2[[1]]
    for(i in 2:length(samples2)){
      merged <- rbind(merged, samples2[[i]])
    }
    
    pred <- merged %*% c(1, newpatient)

    if(ipd$response == "normal"){
      CI <- quantile(pred, probs = quantiles)
    } else if(ipd$response == "binomial"){
      CI <- exp(quantile(pred, probs = quantiles))
    }
    names(CI) <- quantiles
  } else if (class(ipd) == "ipdma.onestage.deft"){
    
    if(is.null(reference)){
      stop("Need to specify reference group for deft approach")
    }
    newpatient <- newpatient - reference
    
    index0 <- which(colnames(samples[[1]]) == "delta[2]") 
    index1 <- grep("gamma.within", colnames(samples[[1]]))
    index <- c(index0, index1)
    samples2 <- samples[,index]
    
    merged <- samples2[[1]]
    for(i in 2:length(samples2)){
      merged <- rbind(merged, samples2[[i]])
    }
    
    pred <- merged %*% c(1, newpatient)

    if(ipd$response == "normal"){
      CI <- quantile(pred, probs = quantiles)
    } else if(ipd$response == "binomial"){
      CI <- exp(quantile(pred, probs = quantiles)) 
    }
    names(CI) <- quantiles
    
  } else{
    stop("Calculating patient specific treatment effect is not yet implemented for this method")
  }
  
  return(CI)
}



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


