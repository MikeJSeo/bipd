
#' Generate a simulated IPD-MA data with systematically missing covariates 
#'
#' Generate a simulated IPD-MA data with systematically missing covariates
#' @param Nstudies number of studies. Default is 10.
#' @param Ncov number of covariates in total. Options are 5 or 10 studies. Default is set to 5.
#' @param sys_missing_prob probability of systematically missing studies for each covariates. Default is set to 0.3.
#' @param magnitude magnitude of the regression estimates (the mean). Default is set to 0.2.
#' @param heterogeneity heterogeneity of regression estimates across studies. Default is set to 0.1.
#' @param interaction whether to include treatment indicator and treatment
#' @return returns simulated IPD-MA data with systematically missing covariates
#' @examples
#' simulated_dataset <- generate_sysmiss_ipdma_example(Nstudies = 10, Ncov = 5, sys_missing_prob = 0.3, 
#' magnitude = 0.2, heterogeneity = 0.1)
#' head(simulated_dataset)
#' @export

generate_sysmiss_ipdma_example <- function(Nstudies = 10, Ncov = 5, sys_missing_prob = 0.1, magnitude = 0.3,
                                           heterogeneity = 0.1, interaction = TRUE) {

  Npatients <- sample(150:300, Nstudies, replace = TRUE)
  Npatients.tot <- sum(Npatients)
  study <- rep(1:Nstudies, times = Npatients)
  
  a <- stats::runif(Nstudies, 0.5, 1.5)
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
    mu <- stats::runif(Ncov, -1, 1)
    X <- rbind(X, mvtnorm::rmvnorm(Npatients[i], mu, Omega * sigma2))
  }
  
  #categorize predictors
  if(Ncov == 5){
    X[,2] <- ifelse(X[,2] > 0, 1, 0)
    X[,3] <- ifelse(X[,3] > 0.5, 1, 0)
  } else if(Ncov == 10){
    X[,2] <- ifelse(X[,2] > 0, 1, 0)
    X[,3] <- ifelse(X[,3] > 0, 1, 0)
    X[,8] <- ifelse(X[,8] > 0, 1, 0)
    X[,9] <- ifelse(X[,9] > 0.5, 1, 0)
    X[,10] <- ifelse(X[,10] > 0.5, 1, 0)
  }
  
  e_vec <- stats::rnorm(Npatients.tot, 0, 1) 
  b <- matrix(NA, Npatients.tot, Ncov)
  #b[,1] <- rep(0.2, Npatients.tot)

  for(i in 1:Ncov){
    b_dummy <- stats::rnorm(Nstudies, magnitude, heterogeneity)
    b_dummy <- rep(b_dummy, times = Npatients)
    b[,i] <- b_dummy
  }
  
  if(interaction == TRUE){
    treat <- stats::rbinom(Npatients.tot, 1, 0.5)
    Xinteraction <- X[,1:Ncov] * treat
    cvec <- matrix(NA, Npatients.tot, Ncov)
    
    for(i in 1:Ncov){
      cvec_dummy <- stats::rnorm(Nstudies, magnitude/2, heterogeneity)
      cvec_dummy <- rep(cvec_dummy, times = Npatients)
      cvec[,i] <- cvec_dummy
    }
    
    d <- stats::rnorm(Nstudies, 1, 0.5)
    d <- rep(d, times = Npatients)
  }
  
  if(interaction == FALSE){
    y <- a + apply(X * b, 1, sum) + e_vec    
  } else if(interaction == TRUE){
    y <- a + apply(X * b, 1, sum) + e_vec + d *treat + apply(Xinteraction * cvec, 1, sum)
  }
  
  # introduce systematically missing; first two predictors are always observed
  for(j in 3:Ncov){
    for(i in 1:Nstudies){
      if(stats::rbinom(1, 1, sys_missing_prob) == 1){
        X[study == i,j] <- NA  
      }
    }
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
  
  if(interaction == TRUE){
    dataset <- cbind(dataset, treat = treat)
    dataset <- as_tibble(dataset)
    #return(list(y = y , X = X, study = study, treat = treat, dataset = dataset))
  } else if(interaction == FALSE){
    dataset <- as_tibble(dataset)
    #return(list(y = y , X = X, study = study, dataset = dataset))
  }
  return(dataset)
}
  

#' Generate a simulated IPD-MA data for demonstration
#'
#' Generate a simulated IPD-MA data for demonstration
#' @param type "continuous" for continuous outcome and "binary" for binary outcome
#' @return returns simulated IPD-MA data
#' @examples
#' ds <- generate_ipdma_example(type = "continuous")
#' head(ds)
#' @export

generate_ipdma_example <- function(type = "continuous"){
  
  N <- 100  #number of patients per trial
  N.trials <- 6 #number of trials
  studyid <- rep(1:6, each = N)
  
  if(type == "continuous"){
    
    alpha <- rep(c(11, 8, 10.5, 9.6, 12.9, 15.8), each = N) #study effects
#    delta <- rep(c(-2.95, -2.97, -2.89, -2.91, -2.93, -2.90), each = N) #treatment effects
    delta <- rep(rnorm(6, -2.5, 2), each = N) #treatment effects
    beta1 <- rep(c(0.24, 0.21, 0.20, 0.18, 0.25, 0.22), each = N) #prognostic effects of z1
    gamma1 <- rep(c(-0.9, -0.5, -0.6, -0.7, -0.1, -0.3), each = N) #interaction effects of z1
    beta2 <- rep(c(0.15, 0.21, 0.30, 0.38, 0.45, 0.42), each = N) #prognostic effects of z2
    gamma2 <- rep(c(0.9, 0.5, 0.5, 0.7, 0.1, 0.5), each = N) #interaction effects of z2
    
    z1 <- stats::rnorm(N*N.trials)
    z2 <- stats::rnorm(N*N.trials)
    treat <- rbinom(N*N.trials, 1, 0.5)
    
    y <- round(alpha + delta*treat + beta1*z1 + beta2*z2 + gamma1*z1*treat + gamma2*z2*treat) 
    ds <- as.data.frame(cbind(studyid, treat, z1, z2, y))
    ds$studyid <- as.factor(studyid)
  } else if (type == "binary"){
    
    alpha <- rep(c(0.11, 0.8, 1.05, 0.96, 0.129, 0.158), each = N) #study effects
#    delta <- rep(c(-0.95, -0.97, -0.89, -0.91, -0.93, -0.90), each = N) #treatment effects
    delta <- rep(rnorm(6, -1.0, 2.0), each = N) #treatment effects
    beta1 <- rep(c(0.24, 0.21, 0.20, 0.18, 0.25, 0.22), each = N) #prognostic effects of w1
    gamma1 <- rep(c(-0.9, -0.5, -0.6, -0.7, -0.1, -0.3), each = N) #interaction effects of w1
    beta2 <- rep(c(0.15, 0.21, 0.30, 0.38, 0.45, 0.42), each = N) #prognostic effects of w2
    gamma2 <- rep(c(0.9, 0.5, 0.5, 0.7, 0.1, 0.5), each = N) #interaction effects of w2
    
    w1 <- stats::rnorm(N*N.trials)
    w2 <- stats::rnorm(N*N.trials)
    treat <- rbinom(N*N.trials, 1, 0.5)
    
    expit <- function(x){
      exp(x)/(1+exp(x))
    }
    linearpredictor <- round(alpha + delta*treat + beta1*w1 + beta2*w2 + gamma1*w1*treat + gamma2*w2*treat)
    y <- stats::rbinom(N, 1, expit(linearpredictor))
    ds <- as.data.frame(cbind(studyid, treat, w1, w2, y))
    ds$studyid <- as.factor(studyid)
  }  

  return(ds)
}


#' Generate a simualted IPD-NMA data for demonstration
#' 
#' Generate a simulated IPD-NMA data for demonstration
#' @param type "continuous" for continuous outcome and "binary" for binary outcome
#' @return return simulated IPD-NMA data
#' ds <- generate_ipdnma_example(type = "continuous")
#' head(ds)
#' @export

generate_ipdnma_example <- function(type = "continuous"){
  
  N <- 100  #number of patients per trial
  N.trials <- 10 #number of trials
  studyid <- rep(1:10, each = N)
  
  treat <- c(sample(c(1,2), size = 300, replace = TRUE), sample(c(1,3), size = 200, replace = TRUE),
             sample(c(2,3), size = 200, replace = TRUE), sample(c(1,2,3), size = 300, replace = TRUE))
  
  if(type == "continuous"){
    
    alpha <- rep(c(11, 8, 10.5, 9.6, 12.9, 15.8, 10.3, 11.2, 10.2, 9.2), each = N) #study effects
    delta12 <- rep(c(-2.95, -2.97, -2.89, -2.88, -2.73, -2.59, -2.89, -3.21, -3.01, -2.94), each = N) #treatment effects
    delta13 <- rep(c(-1.21, -1.1, -1.05, -1.11, -1.30, -0.95, -1.0, -1.2, -1.17, -1.01), each = N)
    
    beta1 <- rep(c(0.24, 0.21, 0.20, 0.18, 0.25, 0.22, 0.19, 0.22, 0.23, 0.25), each = N) #prognostic effects of z1
    gamma1_12 <- rep(c(-0.9, -0.5, -0.6, -0.7, -0.1, -0.3, -0.5, -0.6, -0.7, -0.8), each = N) #interaction effects of z1
    gamma1_13 <- rep(c(-0.1, -0.2, -0.4, -0.6, -0.1, -0.3, -0.2, -0.5, -0.4, -0.3), each = N)
    
    beta2 <- rep(c(0.15, 0.21, 0.30, 0.38, 0.45, 0.42, 0.33, 0.34, 0.36, 0.37), each = N) #prognostic effects of z2
    gamma2_12 <- rep(c(0.9, 0.5, 0.5, 0.7, 0.1, 0.5, 0.4, 0.8, 0.5, 0.6), each = N) #interaction effects of z2
    gamma2_13 <-rep(c(0.8, 0.5, 0.4, 0.5, 0.2, 0.5, 0.1, 0.8, 0.2, 0.1), each = N)
    
    z1 <- stats::rnorm(N*N.trials)
    z2 <- stats::rnorm(N*N.trials)
    
    y <- round(alpha + delta12*(treat==2) + delta13*(treat==3) + beta1*z1 + beta2*z2 + 
                 gamma1_12*z1*(treat==2) + gamma1_13*z1*(treat==3) + gamma2_12*z2*(treat==2) + gamma2_13*z2*(treat==3)) 
    ds <- as.data.frame(cbind(studyid, treat, z1, z2, y))
    ds$studyid <- as.factor(studyid)
    
  } else if(type == "binary"){
    
    alpha <- rep(c(0.11, 0.8, 1.05, 0.96, 0.129, 0.158, 0.13, 0.15, 0.33, 0.50), each = N) #study effects
    delta12 <- rep(c(-0.95, -0.97, -0.89, -0.91, -0.93, -0.90, -1.01, -0.99, -1.05, -0.95), each = N) #treatment effects
    delta13 <- rep(c(-0.55, -0.53, -0.49, -0.76, -0.53, -0.52, -0.52, -0.44, -0.50, -0.39), each = N)
    
    beta1 <- rep(c(0.24, 0.21, 0.20, 0.18, 0.25, 0.22, 0.22, 0.12, 0.15, 0.24), each = N) #prognostic effects of w1
    gamma1_12 <- rep(c(-0.9, -0.5, -0.6, -0.7, -0.1, -0.3, -0.5, -0.6, -0.7, -0.8), each = N) #interaction effects of w1
    gamma1_13 <- rep(c(-0.1, -0.2, -0.4, -0.6, -0.1, -0.3, -0.2, -0.5, -0.4, -0.3), each = N)
    
    beta2 <- rep(c(0.15, 0.21, 0.30, 0.38, 0.45, 0.42, 0.33, 0.34, 0.22, 0.13), each = N) #prognostic effects of w2
    gamma2_12 <- rep(c(0.9, 0.5, 0.5, 0.7, 0.1, 0.5, 0.4, 0.8, 0.5, 0.6), each = N) #interaction effects of w2
    gamma2_13 <-rep(c(0.8, 0.5, 0.4, 0.5, 0.2, 0.5, 0.1, 0.8, 0.2, 0.1), each = N)
    
    w1 <- stats::rnorm(N*N.trials)
    w2 <- stats::rnorm(N*N.trials)
    treat <- rbinom(N*N.trials, 1, 0.5)
    
    expit <- function(x){
      exp(x)/(1+exp(x))
    }
    linearpredictor <- round(alpha + delta12*(treat==2) + delta13*(treat==3) + beta1*w1 + beta2*w2 + 
                      gamma1_12*w1*(treat==2) + gamma1_13*w1*(treat==3) + gamma2_12*w2*(treat==2) + gamma2_13*w2*(treat==3)) 
    y <- stats::rbinom(N, 1, expit(linearpredictor))
    ds <- as.data.frame(cbind(studyid, treat, w1, w2, y))
    ds$studyid <- as.factor(studyid)
  }
  return(ds)
}
