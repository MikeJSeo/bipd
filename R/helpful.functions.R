
#' Calculate patient-specific treatment effect
#'
#' Function for calculating the patient-specific treatment effect.
#' Patient-specific treatment effect includes the main effect of treatment and 
#' treatment-covariate interaction effect (i.e. effect modification). 
#' Reports odds ratio for the binary outcome.
#' 
#' @param ipd IPD object created from running ipdma.model type function
#' @param samples MCMC samples found from running ipd.run function
#' @param newpatient covariate values of patients that you want to predict treatment effect on. Must have length equal to total number of covariates.
#' @param scale_mean option to specify different overall mean compared to what was calculated in IPD object. can be useful when using multiple imputation.
#' @param scale_sd option to specify different overall standard deviation compared to what was calculated in IPD object.
#' @param reference reference group used for finding patient-specific treatment effect. This is only used for deft approach
#' @param quantiles quantiles for credible interval of the patient-specific treatment effect
#' @return patient-specific treatment effect with credible interval at specified quantiles
#' @examples
#' ds <- generate_ipdma_example(type = "continuous")
#' ipd <- with(ds, ipdma.model.onestage(y = y, study = studyid, treat = treat, X = cbind(z1, z2), 
#' response = "normal", shrinkage = "none"))
#' \donttest{
#' samples <- ipd.run(ipd, pars.save = c("beta", "gamma", "delta"), n.chains = 3, n.burnin = 500, 
#' n.iter = 5000)
#' treatment.effect(ipd, samples, newpatient = c(1,0.5))
#' }
#' @references Seo M, White IR, Furukawa TA, et al. Comparing methods for estimating patient-specific treatment effects in individual patient data meta-analysis. \emph{Stat Med}. 2021;40(6):1553-1573. \doi{10.1002/sim.8859}
#' @references Riley RD, Debray TP, Fisher D, et al. Individual participant data meta-analysis to examine interactions between treatment effect and participant-level covariates: Statistical recommendations for conduct and planning. \emph{Stat Med}. 2020:39(15):2115-2137. \doi{10.1002/sim.8516}
#' 
#' @export

treatment.effect <- function(ipd = NULL, samples = NULL, newpatient = NULL, scale_mean = NULL, scale_sd = NULL,
                             reference = NULL, quantiles = c(0.025, 0.5, 0.975)){

  if(is.null(scale_mean)){
    scale_mean <- ipd$scale_mean
  }
  if(is.null(scale_sd)){
    scale_sd <- ipd$scale_sd
  }
    
  if(inherits(ipd, "ipdma.onestage")){
    
    newpatient <- (newpatient - scale_mean)/scale_sd
      
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
  } else if(inherits(ipd, "ipdma.onestage.deft")){
    
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
    
  } else if(inherits(ipd, "ipdnma.onestage")){

    newpatient <- (newpatient - scale_mean)/scale_sd
    
    store_result <- list()
    for(ii in 2:ipd$data.JAGS$Ntreat){
      
      index0 <- which(colnames(samples[[1]]) == paste0("delta[", ii, "]"))  
      index1 <- grep(paste0("gamma\\[", ii, ","), colnames(samples[[1]]))
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
      store_result[[paste0("treatment ", ii)]] <- CI
    }
    return(store_result)
  }
  
  return(CI)
}


#' Convenient function to add results (i.e. combine mcmc.list)
#'
#' This is a convenient function to add results (i.e. combine mcmc.list). 
#' This can be useful when combining results obtained from multiple imputation
#' 
#' @param x first result in a format of mcmc.list
#' @param y second result in a format of mcmc.list
#' @examples
#' ds <- generate_ipdma_example(type = "continuous")
#' ds2 <- generate_ipdma_example(type = "continuous")
#' ipd <- with(ds, ipdma.model.onestage(y = y, study = studyid, treat = treat, X = cbind(z1, z2), 
#' response = "normal", shrinkage = "none"))
#' ipd2 <- with(ds2, ipdma.model.onestage(y = y, study = studyid, treat = treat, X = cbind(z1, z2), 
#' response = "normal", shrinkage = "none"))
#' \donttest{
#' samples <- ipd.run(ipd, pars.save = c("beta", "gamma", "delta"), n.chains = 3, n.burnin = 500, 
#' n.iter = 5000)
#' samples2 <- ipd.run(ipd2, pars.save = c("beta", "gamma", "delta"), n.chains = 3, n.burnin = 500, 
#' n.iter = 5000)
#' combined <- add.mcmc(samples, samples2)
#' }
#' 
#' @export

add.mcmc <- function(x,y){
  
  n.chains <- length(x)
  n.var <- nvar(x)
  newobjects <- vector("list", length = n.chains)
  
  for(i in 1:n.chains){
    newobjects[[i]] <- matrix(NA, nrow = 0, ncol = n.var, dimnames = list(NULL, dimnames(x[[1]])[[2]]))
    newobjects[[i]] <- rbind(x[[i]], y[[i]])
    newobjects[[i]] <- mcmc(newobjects[[i]])
  }
  mcmc.list(newobjects)
}

