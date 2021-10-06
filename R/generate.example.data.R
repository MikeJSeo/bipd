
#' Generate a fake IPD-MA data with systematically missing covariates 
#'
#' Generate a fake IPD-MA data with systematically missing covariates
#' @param Nstudies Number of studies. Default is 10.
#' @param Ncov Number of covariates in total. Options are 5 or 10 studies. Default is set to 5.
#' @param sys_missing_prob Proportion of systematically missing studies for each covariate. Default is set to 0.3.
#' @param signal Signal to noise ratio between predictors and the outcome by changing the variance of the random error. Options are "small" or "large".
#' @param interaction Whether to include treatment indicator and treatment 
#'
#' @export

generate_sysmiss_ipdma_example <- function(Nstudies = 10, Ncov = 5, sys_missing_prob = 0.3,
                                          signal = "small", sign = "different", interaction = FALSE) {

  Npatients <- sample(150:500, Nstudies, replace = TRUE)
  Npatients.tot <- sum(Npatients)
  study <- rep(1:Nstudies, times = Npatients)
  
  a <- runif(Nstudies, 0.5, 1.5)
  a <- rep(a, times = Npatients)
  
  ### generate X
  rho <- 0.2
  Omega <- diag(1, Ncov)
  for(i in 1:Ncov){
    for(j in 1:Ncov){
      Omega[i,j] <- rho^abs(i - j) 
    }
  }
  sigma2 <- 1
  
  X <- NULL
  for(i in 1:Nstudies){
    mu <- runif(Ncov, -0.5, 0.5)
    X <- rbind(X, rmvnorm(Npatients[i], mu, Omega * sigma2))
  }
  
  #categorize predictors
  if(Ncov == 5){
    X[,2] <- ifelse(X[,2] > 0, 1, 0)
    X[,3] <- ifelse(X[,3] > 1, 1, 0)
  } else if(Ncov == 10){
    X[,2] <- ifelse(X[,2] > 0, 1, 0)
    X[,3] <- ifelse(X[,3] > 0, 1, 0)
    X[,8] <- ifelse(X[,8] > 0, 1, 0)
    X[,9] <- ifelse(X[,9] > 0.5, 1, 0)
    X[,10] <- ifelse(X[,10] > 1, 1, 0)
  }
  
  if(signal == "small"){
    e_vec <- rnorm(Npatients.tot, 0, 1.0) # R squared around 0.1
  } else if(signal == "large"){
    e_vec <- rnorm(Npatients.tot, 0, 0.2) # R squared around 0.6
  }
  
  b <- matrix(NA, Npatients.tot, Ncov)
  
  for(i in 1:Ncov){
    if(sign == "different"){
      b_dummy <- abs(rnorm(Nstudies, 0.5, 0.1))
      if(Nstudies == 10){
        b_dummy[c(1,3,5,7,9)] <- -b_dummy[c(1,3,5,7,9)]
      } else{
        stop("Need to change code to allow Nstudies other than 10")
      }
    } else if (sign == "same"){
      b_dummy <- abs(rnorm(Nstudies, 0.5, 0.1))
    }
    b_dummy <- rep(b_dummy, times = Npatients)
    b[,i] <- b_dummy
  }
  
  y <- a + apply(X * b, 1, sum) + e_vec  
  
  # introduce systematically missing; first two predictors are always observed
  for(i in 3:Ncov){
    systematic_missing_study <- sample(Nstudies, size = sys_missing_prob*Nstudies)
    systematic_missing_study_dummy <- rep(0, Nstudies)
    systematic_missing_study_dummy[systematic_missing_study] <- 1
    
    study_missing <- as.logical(rep(systematic_missing_study_dummy, times = Npatients))  
    X[,i][study_missing] <- NA
  }    
  
  
  # Define dataset to return
  dataset <- data.frame(y = y, X = X[,1:Ncov], study = study)
  colnames(dataset) <- c("y", paste0("x", 1:Ncov), "study")
  dataset <- as_tibble(dataset)
  
  if(Ncov == 5){
    dataset <- dataset %>% mutate(x2 = as.factor(x2),
                                  x3 = as.factor(x3))
  } else if(Ncov == 10){
    dataset <- dataset %>% mutate(x2 = as.factor(x2),
                                  x3 = as.factor(x3),
                                  x8 = as.factor(x8),
                                  x9 = as.factor(x9),
                                  x10 = as.factor(x10)
    )
  }  
  
  dataset <- as_tibble(dataset)
  return(list(y = y , X = X, study = study, dataset = dataset))
  
}
  
  
  



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