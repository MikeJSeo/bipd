## Code to generate simulation
  
generate.simulation <- function(Nstudies = NULL, Ncovariate = NULL, continuous.cov = NULL, pf = NULL, em = NULL,
                                beta = NULL, gamma = NULL, sampleSize = c(50, 100), model.type = "gaussian", 
                                tau = 0.2, tau_beta = 0.2, tau_gamma = 0.3){
  
  #treatment effect
  delta <- 1
  d <- rnorm(Nstudies, delta, tau)
  
  # prognostic effect (covariate effect)
  b <- matrix(NA, nrow = Nstudies, ncol = length(beta))
  for(i in 1:length(beta)){
    b[,i] <- rnorm(Nstudies, beta[i], tau_beta)
  }
  
  # effect modification (treatment-covariate interaction)
  if(!is.null(gamma)){
    c <- matrix(NA, nrow = Nstudies, ncol = length(gamma))
    for(i in 1:length(gamma)){
      c[,i] <- rnorm(Nstudies, gamma[i], tau_gamma)
    }   
  }
  
  studyid <- NULL
  for(i in 1:Nstudies){
    studyid <- c(studyid, rep(i, sample(sampleSize[1]:sampleSize[2], size = 1)))
  }
  
  #study baseline effect
  if(model.type == "gaussian"){
    alpha <- runif(Nstudies, -1, 1)
  } else{
    alpha <- runif(Nstudies, -2, -1) 
  }
  treat <- rbinom(length(studyid), 1, 0.5)
  
  #generating data
  rho <- 0.3
  len <- length(continuous.cov)
  
  cov_matrix <- matrix(NA, nrow = len, ncol = len)
  for(ii in 1:len){
    for(jj in 1:len){
      cov_matrix[ii,jj] <- rho^abs(ii - jj) 
    }
  }
  X <- matrix(NA, nrow = length(studyid), ncol = Ncovariate)
  for(i in 1:length(studyid)){
    X[i,continuous.cov] <-  mvrnorm(n = 1, mu = rep(0, len), cov_matrix)
  }
  X[,-continuous.cov] <- rbinom(length(studyid)* (Ncovariate - length(continuous.cov)), 1, 0.5)
  
  # standardize X: binary and continuous variables in same scale
  X <- apply(X, 2, scale)
  data <- model.matrix(~ -1 + X*treat)
  
  meany <- alpha[studyid] + d[studyid] * treat + apply(X[,pf, drop = FALSE] * b[studyid,], 1, sum)
  
  if(!is.null(em)){
    meany <- meany + apply(X[,em, drop = FALSE] * c[studyid,] * treat, 1, sum)
  }
  
  sigmay <- 0.5
  py <- expit(meany)

  if(model.type == "gaussian"){
    y <- rnorm(length(studyid), meany, sigmay)
  } else if (model.type == "binary"){
    y <- rbinom(length(studyid), 1, py)
  }

  data <- cbind(y = y, data = data, studyid = studyid)
  data <- as.data.frame(data)
  data$studyid <- as.factor(data$studyid)
  return(data)
}

expit <- function(x){
  exp(x)/(1+exp(x))
}

calc_mse <- function(a, b){
  mean((a - b)^2)
}



find_performance <- function(val, correct_values, correct_em, data_subset){
  
  val_without_treat <- val[-length(val)]
  val_treat <- val[length(val)]
  c(calc_mse(val_without_treat[correct_em != 1], correct_values[correct_em != 1]),
    calc_mse(val_without_treat[correct_em == 1], correct_values[correct_em == 1]),
    calc_mse(val_treat, 1),
    mean((data_subset %*% val - data_subset %*% c(correct_values, 1))^2)
    )
}

find_performance1 <- function(val, correct_values, data_subset){
  
  bias <- mean((data_subset %*% val - data_subset %*% c(correct_values, 1)))
  
  diff <- data_subset %*% val - data_subset %*% c(correct_values, 1)
  variance <- mean((diff - mean(diff))^2)
#  variance <- var((data_subset %*% val - data_subset %*% c(correct_values, 1)))
  
  c(bias,bias^2, variance)
}


find_performance2 <- function(val, correct_em, continuous.cov){
  
  val_treat <- val[length(val)]
  
  continuous.indicator <- rep(0, length(correct_em))
  continuous.indicator[continuous.cov] <- 1
  
  val_without_treat <- val[-length(val)]
  true_em_value_continuous <- val_without_treat[correct_em == 1 & continuous.indicator == 1]
  true_em_value_continuous <- true_em_value_continuous[true_em_value_continuous != 0]
  
  true_em_value_binary <- val_without_treat[correct_em == 1 & continuous.indicator == 0]
  true_em_value_binary <- true_em_value_binary[true_em_value_binary != 0]
  
  c(ifelse(length(true_em_value_continuous) == 0, NA, mean(true_em_value_continuous)),
    ifelse(length(true_em_value_binary) == 0, NA, mean(true_em_value_binary)),
    val_treat)
}


bootstrap_function  <- function(model_data, ndraws, p.fac, family, alpha = 1) {
  
  coeff_mtx <- matrix(0, nrow = ndraws, ncol = length(col_labels))
  
  for (ii in 1:ndraws) {
    
    bootstrap_ids <- sample(seq(nrow(model_data)), nrow(model_data), replace = TRUE)
    bootstrap_data <- model_data[bootstrap_ids,]
    
    bootstrap_model <- cv.glmnet(as.matrix(bootstrap_data[,-1]), as.matrix(bootstrap_data[,1]), penalty.factor = p.fac, family = family, standardize = FALSE, type.measure = "deviance", alpha = alpha)  
    aa <- coef(bootstrap_model, s = "lambda.min")
    coeff_mtx[ii,]   <- sapply(col_labels, function(x) ifelse(x %in% rownames(aa)[aa[,1] != 0], aa[x,1], 0))  
  }
  se <- apply(coeff_mtx, 2, sd, na.rm = TRUE)
  return(se)
}

