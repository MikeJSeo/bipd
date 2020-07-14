##########################################################################
#helpful functions that were used to run/evaluate IPD-MA models (one-stage)
##########################################################################

#' Generate a fake IPD-MA data for demonstration
#'
#' Generate a fake IPD-MA data for demonstration
#' @param type "continuous" for continuous outcome and "binary" for binary outcome
#'
#' @export

generate_ipdma_example <- function(type = "continuous"){

  if(type == "continuous"){
  
  # continuous outcome
  N <- 100  #number of patients per trial
  N.trials <- 6 #number of trials
  alpha <- c(11, 8, 10.5, 9.6, 12.9, 15.8) #study effects
  delta <- c(-2.95, -2.97, -2.89, -2.91, -2.93, -2.90) #treatment effects
  beta1 <- c(0.24, 0.21, 0.20, 0.18, 0.25, 0.22) #prognostic effects of z1
  gamma1 <- c(-0.9, -0.5, -0.6, -0.7, -0.1, -0.3) #interaction effects of z1
  beta2 <- c(0.15, 0.21, 0.30, 0.38, 0.45, 0.42) #prognostic effects of z2
  gamma2 <- c(0.9, 0.5, 0.5, 0.7, 0.1, 0.5) #interaction effects of z2
  studyid <- c(1:6)
  ds <- as.data.frame(array(NA, dim=c(N*6, 5)))
  colnames(ds) <- c("studyid", "treat", "z1","z2", "y")
  for (i in 1:N.trials) {
    treat <- rbinom(N,1,0.5)
    z1 <- rnorm(N, mean=0, sd=1)
    z2 <- rnorm(N, mean=0, sd=1)
    y <- round(alpha[i] + delta[i]*treat + beta1[i]*z1 + beta2[i]*z2 + gamma1[i]*z1*treat + gamma2[i]*z2*treat)
    ds[(((i-1)*N)+1):(i*N), ] <- cbind(studyid[i], treat, z1,z2, y)
  }
  ds$studyid <- as.factor(ds$studyid)

  } else if(type == "binary"){
  
  # dichotomous outcome
  N <- 100  #number of patients per trial
  N.trials <- 6 #number of trials
  alpha <- c(0.11, 0.8, 1.05, 0.96, 0.129, 0.158) #study effects
  delta <- c(-0.95, -0.97, -0.89, -0.91, -0.93, -0.90) #treatment effects
  beta1 <- c(0.24, 0.21, 0.20, 0.18, 0.25, 0.22) #prognostic effects of w1
  gamma1 <- c(-0.9, -0.5, -0.6, -0.7, -0.1, -0.3) #interaction effects of w1
  beta2 <- c(0.15, 0.21, 0.30, 0.38, 0.45, 0.42) #prognostic effects of w2
  gamma2 <- c(0.9, 0.5, 0.5, 0.7, 0.1, 0.5) #interaction effects of w2
  studyid <- c(1:6)
  ds <- as.data.frame(array(NA, dim=c(N*6, 5)))
  colnames(ds) <- c("studyid", "treat", "w1","w2", "y")
  
  expit <- function(x){
    exp(x)/(1+exp(x))
  }
  
  for (i in 1:N.trials) {
    treat <- rbinom(N,1,0.5)
    w1 <- rnorm(N, mean=0, sd=1) 
    w2 <- rnorm(N, mean=0, sd=1)
    linearpredictor <- alpha[i] + delta[i]*treat + beta1[i]*w1 + beta2[i]*w2 + gamma1[i]*w1*treat + gamma2[i]*w2*treat
    y <- rbinom(N, 1, expit(linearpredictor))
    ds[(((i-1)*N)+1):(i*N), ] <- cbind(studyid[i], treat, w1, w2, y)
  }
  ds$studyid <- as.factor(ds$studyid)
    
  }
  
  return(ds)
}



#' Calculate patient-specific treatment effect
#'
#' Convenient function for calculating the patient-specific treatment effect.
#' Patient-specific treatment effect includes the main effect of treatment and 
#' treatment-covariate interaction effect (i.e. effect modification)
#' 
#' @param ipd ipd object created from ipd.data function
#' @param samples mcmc samples found from running `ipd.run`
#' @param newpatient Covariate values of patients that you want to predict treatment effect on. Must have length equal to total number of covariates.
#' @param response "normal" for continuous outcome and "binomial" for binary outcome. Reports odds ratio for the binary outcome.
#' @param quantile quantile for the confidence interval
#'
#' @export

treatment.effect <- function(ipd = NULL, samples = NULL, newpatient = NULL, 
                             response = "normal", quantile = c(0.025, 0.5, 0.975)){

  if(!is.null(ipd$scale_mean)) {
    newpatient <- (newpatient - ipd$scale_mean)/ipd$scale_sd
  }
  
  index0 <- which(colnames(samples[[1]]) == "delta[2]") #only works for IPD-MA
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
  
  if(response == "normal"){
    CI <- mean1 + qnorm(quantile)* se1
  } else if(response == "binomial"){
    CI <- exp(mean1 + qnorm(quantile)* se1)    
  }
  names(CI) <- quantile
  return(CI)
}


