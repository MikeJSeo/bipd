
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
#' @param quantile Quantiles finding confidence interavl of the patient-specific treatment effect
#' @references Seo M, White IR, Furukawa TA, et al. Comparing methods for estimating patient-specific treatment effects in individual patient data meta-analysis. \emph{Stat Med}. 2021;40(6):1553-1573. \doi{10.1002/sim.8859}
#' @references Riley RD, Debray TP, Fisher D, et al. Individual participant data meta-analysis to examine interactions between treatment effect and participant-level covariates: Statistical recommendations for conduct and planning. \emph{Stat Med}. 2020:39(15):2115-2137. [\url{https://doi.org/10.1002/sim.8516}] 
#' @export

treatment.effect <- function(ipd = NULL, samples = NULL, newpatient = NULL, 
                             reference = NULL, quantile = c(0.025, 0.5, 0.975)){

  
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
    mean1 <- mean(pred)
    se1 <- sd(pred)
    
    if(ipd$response == "normal"){
      CI <- mean1 + qnorm(quantile)* se1
    } else if(ipd$response == "binomial"){
      CI <- exp(mean1 + qnorm(quantile)* se1)    
    }
    names(CI) <- quantile
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
    mean1 <- mean(pred)
    se1 <- sd(pred)
    
    if(ipd$response == "normal"){
      CI <- mean1 + qnorm(quantile)* se1
    } else if(ipd$response == "binomial"){
      CI <- exp(mean1 + qnorm(quantile)* se1)    
    }
    names(CI) <- quantile
    
  } else{
    stop("Calculating patient specific treatment effect is not yet implemented for this method")
  }
  
  return(list(CI = CI, pred = pred))
}

