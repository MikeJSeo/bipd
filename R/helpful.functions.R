##########################################################################
#helpful functions that were used to run/evaluate models
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
#' treatment-covariate interaction effect (i.e. effect modification). 
#' Reports odds ratio for the binary outcome.
#' 
#' @param ipd ipd object created from ipd.data function
#' @param samples mcmc samples found from running `ipd.run`
#' @param newpatient Covariate values of patients that you want to predict treatment effect on. Must have length equal to total number of covariates.
#' @param reference reference group used for finding patient-specific treatment effect; only used for "deft" approach
#' @param quantile quantile for the confidence interval
#' @references Riley RD, Debray TP, Fisher D, et al. Individual participant data meta-analysis to examine interactions between treatment effect and participant-level covariates: Statistical recommendations for conduct and planning. \emph{Stat Med}. 2020:39(15):2115-2137. [\url{https://doi.org/10.1002/sim.8516}] 
#' @export

treatment.effect <- function(ipd = NULL, samples = NULL, newpatient = NULL, 
                             reference = NULL, quantile = c(0.025, 0.5, 0.975)){

  if(!is.null(ipd$scale_mean)) {
    newpatient <- (newpatient - ipd$scale_mean)/ipd$scale_sd
  }
  
  if(ipd$approach == "deft"){
    if(is.null(reference)){
      stop("Need to specify reference group for deft approach")
    }
    newpatient <- newpatient - reference
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
  
  if(ipd$response == "normal"){
    CI <- mean1 + qnorm(quantile)* se1
  } else if(ipd$response == "binomial"){
    CI <- exp(mean1 + qnorm(quantile)* se1)    
  }
  names(CI) <- quantile
  return(CI)
}



#' Generate a fake IPD-MA data with systematically missing variables for demonstration
#'
#' Generate a fake IPD-MA data with systematically missing variables for demonstration
#'
#' @export

generate_ipdma_example_with_systematic_missing <- function(){
  
  set.seed(1)
  Nstudies <- 50
  Npatients <- 1500
  Npatients.tot <- Nstudies*Npatients
  study <- rep(1:Nstudies, each = Npatients)
  treat <- rbinom(Npatients.tot, 1, 0.5)
  
  a <- 1
  b <- c(0.7, 1, 0.7, 0.5)
  c <- c(0.1, 0.5, 0.2, 0.2)
  d <- 0.5
  Sigma <- matrix(c(0.2^2, -0.1*0.2*0.2, -0.1*0.2*0.2, 0.2^2), nrow = 2)
  
  u <- rmvnorm(Nstudies, rep(0, 2), matrix(c(0.2^2, 0.2 * 0.2 * -0.1, 0.2 * 0.2 * -0.1, 0.2^2), nrow = 2) )
  
  #generate x1
  gamma1 <- rnorm(Nstudies, 0, 0.5)
  e1 <- rnorm(Npatients.tot, 0, 0.2)
  x1 <- rep(gamma1, each = Npatients)+ e1
  
  #generate x2
  gamma2 <- rmvnorm(Nstudies, c(0, 0), matrix(c(0.5^2, 0.2 * 0.5 * 0.5, 0.2 * 0.5 * 0.5, 0.5^2), nrow = 2))
  e2 <- rnorm(Npatients.tot, 0, 0.1)
  x2 <- 0.3 * x1 + gamma2[,1] + gamma2[,2]*x1 + e2
  #  gamma2 <- rnorm(Nstudies, 0, 0.5)
  #  x2 <- rep(gamma2, each = Npatients) + e2
  
  #generate x3
  p3 <- runif(Nstudies, 0.05, 0.15)
  x3 <- rbinom(Npatients.tot, 1, rep(p3, each = Npatients))
  
  #generate x4
  p4 <- runif(Nstudies, 0.15, 0.25)
  x4 <- rbinom(Npatients.tot, 1, rep(p4, each = Npatients))
  
  #generate y
  y <- rep(a, Npatients.tot) + b[1] * x1 + b[2] * x2 + b[3] * x3 + b[4] * x4 + 
    c[1] * x1 * treat + c[2] * x2 * treat + c[3] * x3 * treat + c[4] * x4 * treat +
    d * treat + u[,1] + u[,2] * treat
  
  #generate systematically missing framework in x2
  pi_sys <- 0.2
  study_missing <- as.logical(rep(rbinom(Nstudies, 1, 0.2), each = Npatients))
  x2[study_missing] <- NA
  
  #generate sporadically missing framework for all variables
  X <- cbind(x1, x2, x3, x4)
  X[as.logical(rbinom(Npatients.tot * 4, 1, 0.05))] <- NA
  
  return(as_tibble(cbind(y, X, study, treat)))
}