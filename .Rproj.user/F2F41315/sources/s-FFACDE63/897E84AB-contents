#########
# run.simulation runs glmm (full and oracle) and stepwise regression
# run.simulation2 runs lasso and ridge
# run.simulation3 runs Bayesian methods (Bayesian LASSO and SSVS)

run.simulation <- function(){
  
  step_store_mse <- glmm_full_store_mse <- glmm_oracle_store_mse <- matrix(NA, nrow = niter, ncol = 4)
  step_store_sd <- glmm_full_store_sd <- glmm_oracle_store_sd <- matrix(NA, nrow = niter, ncol = 3)
  step_store_bias <- glmm_full_store_bias <- glmm_oracle_store_bias <- matrix(NA, nrow = niter, ncol = 3)
  
  for(i in seq(niter)){
    
    set.seed(i)
    data <-generate.simulation(Nstudies = Nstudies, Ncovariate = Ncovariate, continuous.cov = continuous.cov, pf = pf, em = em, beta = beta, gamma = gamma, model.type = model.type, tau = tau)
    if(model.type == "gaussian"){
      m1 <- lmer(glmm_oracle_formula, data = data)
    } else if(model.type == "binary"){
      m1 <- glmer(glmm_oracle_formula, data = data, family = binomial(link = "logit"))
    }
    summary(m1)
    
    mean_values <- sapply(col_labels, function(x) ifelse(x %in% names(fixef(m1)), summary(m1)$coef[x,"Estimate"], 0))
    sd_values <- sapply(col_labels, function(x) ifelse(x %in% names(fixef(m1)), summary(m1)$coef[x,"Std. Error"], 0))
    
    glmm_oracle_store_mse[i,] <- find_performance(mean_values, correct_em_values, correct_em, as.matrix(data[,col_labels]))
    glmm_oracle_store_sd[i,] <- find_performance2(sd_values, correct_em, continuous.cov)
    glmm_oracle_store_bias[i,] <- find_performance1(mean_values, correct_em_values, as.matrix(data[,col_labels]))
  }

  for(i in seq(niter)){
    
    set.seed(i)
    data <-generate.simulation(Nstudies = Nstudies, Ncovariate = Ncovariate, continuous.cov = continuous.cov, pf = pf, em = em, beta = beta, gamma = gamma, model.type = model.type, tau = tau)
    if(model.type == "gaussian"){
      m1 <- lmer(glmm_full_formula, data = data)
    } else if(model.type == "binary"){
      m1 <- glmer(glmm_full_formula, data = data, family = binomial(link = "logit"))
    }

    mean_values <- sapply(col_labels, function(x) ifelse(x %in% names(fixef(m1)), summary(m1)$coef[x,"Estimate"], 0))
    sd_values <- sapply(col_labels, function(x) ifelse(x %in% names(fixef(m1)), summary(m1)$coef[x,"Std. Error"], 0))
    
    glmm_full_store_mse[i,] <- find_performance(mean_values, correct_em_values, correct_em, as.matrix(data[,col_labels]))
    glmm_full_store_sd[i,] <- find_performance2(sd_values, correct_em, continuous.cov)
    glmm_full_store_bias[i,] <- find_performance1(mean_values, correct_em_values, as.matrix(data[,col_labels]))
  }
  
  for(i in seq(niter)){
    
    set.seed(i)
    data <-generate.simulation(Nstudies = Nstudies, Ncovariate = Ncovariate, continuous.cov = continuous.cov, pf = pf, em = em, beta = beta, gamma = gamma, model.type = model.type, tau = tau)
    
    if(model.type == "gaussian"){
      m1 <- lm(step_full_formula, data = data)  
    } else if(model.type == "binary"){
      m1 <- glm(step_full_formula, family = binomial(link = "logit"), data = data)
    }
    s1 <- step(m1, scope=list(lower= step_lower_formula))
    
    mean_values <- sapply(col_labels, function(x) ifelse(x %in% variable.names(s1), summary(s1)$coef[x,"Estimate"], 0))
    sd_values <- sapply(col_labels, function(x) ifelse(x %in% variable.names(s1), summary(s1)$coef[x,"Std. Error"], 0))
    
    step_store_mse[i,] <- find_performance(mean_values, correct_em_values, correct_em, as.matrix(data[,col_labels]))
    step_store_sd[i,] <- find_performance2(sd_values, correct_em, continuous.cov)
    step_store_bias[i,] <- find_performance1(mean_values, correct_em_values, as.matrix(data[,col_labels]))
  }
  
  glmm_oracle_store_mse_mean <- apply(glmm_oracle_store_mse, 2, mean)
  glmm_oracle_store_sd_mean <- apply(glmm_oracle_store_sd, 2, mean, na.rm = TRUE)
  glmm_full_store_mse_mean <- apply(glmm_full_store_mse, 2, mean)
  glmm_full_store_sd_mean <- apply(glmm_full_store_sd, 2, mean, na.rm = TRUE)
  step_store_mse_mean <- apply(step_store_mse, 2, mean)
  step_store_sd_mean <- apply(step_store_sd, 2, mean, na.rm = TRUE)

  glmm_oracle_store_bias_mean <- apply(glmm_oracle_store_bias, 2, mean)
  glmm_full_store_bias_mean <- apply(glmm_full_store_bias, 2, mean)
  step_store_bias_mean <- apply(step_store_bias, 2, mean)
  
  result_matrix_mse <- matrix(NA, nrow = 3, ncol = 4)
  colnames(result_matrix_mse) <- c("false em mse", "true em mse","treatment mse", "patient specific trt mse")
  rownames(result_matrix_mse) <-  c("glmm oracle", "glmm full","naive step")
  result_matrix_mse[1,] <- glmm_oracle_store_mse_mean
  result_matrix_mse[2,] <- glmm_full_store_mse_mean
  result_matrix_mse[3,] <- step_store_mse_mean

  result_matrix_sd <- matrix(NA, nrow = 3, ncol = 3)
  colnames(result_matrix_sd) <- c("continuous EM se", "binary EM se","treatment se")
  rownames(result_matrix_sd) <-  c("glmm oracle", "glmm full","naive step")
  result_matrix_sd[1,] <- glmm_oracle_store_sd_mean
  result_matrix_sd[2,] <- glmm_full_store_sd_mean
  result_matrix_sd[3,] <- step_store_sd_mean

  result_matrix_bias <- matrix(NA, nrow = 3, ncol = 3)
  colnames(result_matrix_bias) <- c("bias", "bias^2","variance")
  rownames(result_matrix_bias) <-  c("glmm oracle", "glmm full","naive step")
  result_matrix_bias[1,] <- glmm_oracle_store_bias_mean
  result_matrix_bias[2,] <- glmm_full_store_bias_mean
  result_matrix_bias[3,] <- step_store_bias_mean
  
  cbind(result_matrix_mse, result_matrix_sd, result_matrix_bias)
}


run.simulation2 <- function(){
  
  glmnet_store_mse <- ridge_store_mse <- matrix(NA, nrow = niter, ncol = 4)
  glmnet_store_sd <- ridge_store_sd <- matrix(NA, nrow = niter, ncol = 3)
  glmnet_store_bias <- ridge_store_bias <- matrix(NA, nrow = niter, ncol = 3)
  
  lambdas <- 10^seq(3, -3, by = -.1) #define predefined set of lambdas to cross validate
  
  for(i in seq(niter)){
    
    set.seed(i)
    data <-generate.simulation(Nstudies = Nstudies, Ncovariate = Ncovariate, continuous.cov = continuous.cov, pf = pf, em = em, beta = beta, gamma = gamma, model.type = model.type, tau = tau)
    
    data_glmnet <- model.matrix(step_full_formula, data = data)
    data_glmnet <- data_glmnet[,-1] 
    data_glmnet <- cbind(y = data$y, data_glmnet = data_glmnet)  
    
    p.fac <- c(rep(0, Nstudies - 1), rep(0, Ncovariate), 0, rep(1, Ncovariate))
    
    family <- ifelse(model.type == "gaussian", "gaussian", "binomial")
    
    cvfit <- cv.glmnet(as.matrix(data_glmnet[,-1]), as.matrix(data_glmnet[,1]), penalty.factor = p.fac, family = family, standardize = FALSE, type.measure = "deviance", lambda = lambdas)  
    aa <- coef(cvfit, s = "lambda.min")
    
    mean_values <-  sapply(col_labels, function(x) ifelse(x %in% rownames(aa)[aa[,1] != 0], aa[x,1], 0))
    glmnet_store_mse[i,] <- find_performance(mean_values, correct_em_values, correct_em, data_glmnet[,col_labels])
    glmnet_store_bias[i,] <- find_performance1(mean_values, correct_em_values, data_glmnet[,col_labels])
  }
  
  for(i in seq(niter)){
    
    set.seed(i)
    data <-generate.simulation(Nstudies = Nstudies, Ncovariate = Ncovariate, continuous.cov = continuous.cov, pf = pf, em = em, beta = beta, gamma = gamma, model.type = model.type, tau = tau)
    
    data_glmnet <- model.matrix(step_full_formula, data = data)
    data_glmnet <- data_glmnet[,-1] 
    data_glmnet <- cbind(y = data$y, data_glmnet = data_glmnet)  
    
    p.fac <- c(rep(0, Nstudies - 1), rep(0, Ncovariate), 0, rep(1, Ncovariate))
    
    family <- ifelse(model.type == "gaussian", "gaussian", "binomial")
        
    cvfit <- cv.glmnet(as.matrix(data_glmnet[,-1]), as.matrix(data_glmnet[,1]), penalty.factor = p.fac, family = family, standardize = FALSE, type.measure = "deviance", alpha = 0, lambda = lambdas)  
    aa <- coef(cvfit, s = "lambda.min")
    
    mean_values <-  sapply(col_labels, function(x) ifelse(x %in% rownames(aa)[aa[,1] != 0], aa[x,1], 0))
    ridge_store_mse[i,] <- find_performance(mean_values, correct_em_values, correct_em, data_glmnet[,col_labels])
    ridge_store_bias[i,] <- find_performance1(mean_values, correct_em_values, data_glmnet[,col_labels])
  }
  
  glmnet_store_mse_mean <- apply(glmnet_store_mse, 2, mean)
  glmnet_store_sd_mean <- apply(glmnet_store_sd, 2, mean, na.rm = TRUE)
  ridge_store_mse_mean <- apply(ridge_store_mse, 2, mean)
  ridge_store_sd_mean <- apply(ridge_store_sd, 2, mean, na.rm = TRUE)
  
  glmnet_store_bias_mean <- apply(glmnet_store_bias, 2, mean)
  ridge_store_bias_mean <- apply(ridge_store_bias, 2, mean)
  
  result_matrix_mse <- matrix(NA, nrow = 2, ncol = 4)
  colnames(result_matrix_mse) <- c("false em mse", "true em mse","treatment mse", "patient specific trt mse")
  rownames(result_matrix_mse) <-  c("naive lasso", "ridge")
  result_matrix_mse[1,] <- glmnet_store_mse_mean
  result_matrix_mse[2,] <- ridge_store_mse_mean
  
  result_matrix_sd <- matrix(NA, nrow = 2, ncol = 3)
  colnames(result_matrix_sd) <- c("continuous EM se", "binary EM se","treatment se")
  rownames(result_matrix_sd) <-  c("naive lasso", "ridge")
  result_matrix_sd[1,] <- glmnet_store_sd_mean
  result_matrix_sd[2,] <- ridge_store_sd_mean
  
  result_matrix_bias <- matrix(NA, nrow = 2, ncol = 3)
  colnames(result_matrix_bias) <- c("bias", "bias^2","variance")
  rownames(result_matrix_bias) <-  c("naive lasso", "ridge")
  result_matrix_bias[1,] <- glmnet_store_bias_mean
  result_matrix_bias[2,] <- ridge_store_bias_mean
  
  cbind(result_matrix_mse, result_matrix_sd, result_matrix_bias)
}

run.simulation3 <- function(){

  SSVS_store_mse <- bayesLASSO_store_mse <- matrix(NA, nrow = niter, ncol = 4)
  SSVS_store_sd <- bayesLASSO_store_sd <- matrix(NA, nrow = niter, ncol = 3)
  SSVS_store_bias <- bayesLASSO_store_bias <- matrix(NA, nrow = niter, ncol = 3)
  
  for(i in seq(niter)){
    
    set.seed(i)
    data <-generate.simulation(Nstudies = Nstudies, Ncovariate = Ncovariate, continuous.cov = continuous.cov, pf = pf, em = em, beta = beta, gamma = gamma, model.type = model.type, tau = tau)
    
    data_jags <- with(data,{
      list(Nstudies = length(unique(studyid)),
           X = data[,paste0("X",1:Ncovariate)],
           Np = length(X1),
           Ncovariate = (dim(data)[2] - 3)/2,
           studyid = studyid,
           treat = treat + 1,
           y = y)
    })
    
    if(model.type == "gaussian") {
      samples <- jags.parfit(cl = cl, data = data_jags, params = c("gamma", "delta", "beta"), model = "IPD-MA-bayesLASSO.txt", n.chains = 2, n.adapt = 100, n.update = 200, n.iter = 2000)
    } else if(model.type == "binary"){
      samples <- jags.parfit(cl = cl, data = data_jags, params = c("gamma", "delta", "beta"), model = "IPD-MA-bayesLASSO-binomial.txt", n.chains = 2, n.adapt = 100, n.update = 200, n.iter = 2000)
    }
    
    a <- summary(samples)
    
    g_mean <-  a$statistics[grep("gamma\\[", rownames(a$statistics)), "Mean"]
    treat_mean <- a$statistics["delta[2]", "Mean"]
    names(treat_mean) <- "treat"
    mean_values <- c(g_mean, treat_mean)
    bayesLASSO_store_mse[i,] <- find_performance(mean_values, correct_em_values, correct_em, as.matrix(data[,col_labels]))
    bayesLASSO_store_bias[i,] <- find_performance1(mean_values, correct_em_values, as.matrix(data[,col_labels]))
    
    g_sd <-  a$statistics[grep("gamma\\[", rownames(a$statistics)), "SD"]
    treat_sd <- a$statistics["delta[2]", "SD"]
    names(treat_sd) <- "treat"
    sd_values <- c(g_sd, treat_sd)
    bayesLASSO_store_sd[i,] <- find_performance2(sd_values, correct_em, continuous.cov)
  }

  for(i in seq(niter)){
    
    set.seed(i)
    data <-generate.simulation(Nstudies = Nstudies, Ncovariate = Ncovariate, continuous.cov = continuous.cov, pf = pf, em = em, beta = beta, gamma = gamma, model.type = model.type, tau = tau)
    data_jags <- with(data,{
      list(Nstudies = length(unique(studyid)),
           X = data[,paste0("X",1:Ncovariate)],
           Np = length(X1),
           Ncovariate = (dim(data)[2] - 3)/2,
           studyid = studyid,
           treat = treat + 1,
           y = y)
    })
    
    if(model.type == "gaussian") {
      samples <- jags.parfit(cl = cl, data = data_jags, params = c("gamma", "delta", "beta"), model = "IPD-MA-SSVS.txt", n.chains = 2, n.adapt = 100, n.update = 200, n.iter = 2000)
    } else if(model.type == "binary"){
      samples <- jags.parfit(cl = cl, data = data_jags, params = c("gamma", "delta", "beta"), model = "IPD-MA-SSVS-binomial.txt", n.chains = 2, n.adapt = 100, n.update = 200, n.iter = 2000)
    }
    
    a <- summary(samples)
    
    g_mean <-  a$statistics[grep("gamma\\[", rownames(a$statistics)), "Mean"]
    treat_mean <- a$statistics["delta[2]", "Mean"]
    names(treat_mean) <- "treat"
    mean_values <- c(g_mean, treat_mean)
    SSVS_store_mse[i,] <- find_performance(mean_values, correct_em_values, correct_em, as.matrix(data[,col_labels]))
    SSVS_store_bias[i,] <- find_performance1(mean_values, correct_em_values, as.matrix(data[,col_labels]))
    
    g_sd <-  a$statistics[grep("gamma\\[", rownames(a$statistics)), "SD"]
    treat_sd <- a$statistics["delta[2]", "SD"]
    names(treat_sd) <- "treat"
    sd_values <- c(g_sd, treat_sd)
    SSVS_store_sd[i,] <- find_performance2(sd_values, correct_em, continuous.cov)
  }
  
  bayesLASSO_store_mse_mean <- apply(bayesLASSO_store_mse, 2, mean)
  bayesLASSO_store_sd_mean <- apply(bayesLASSO_store_sd, 2, mean, na.rm = TRUE)
  SSVS_store_mse_mean <- apply(SSVS_store_mse, 2, mean)
  SSVS_store_sd_mean <- apply(SSVS_store_sd, 2, mean, na.rm = TRUE)
  
  bayesLASSO_store_bias_mean <- apply(bayesLASSO_store_bias, 2, mean)
  SSVS_store_bias_mean <- apply(SSVS_store_bias, 2, mean)
  
  result_matrix_mse <- matrix(NA, nrow = 2, ncol = 4)
  colnames(result_matrix_mse) <- c("false em mse", "true em mse","treatment mse", "patient specific trt mse")
  rownames(result_matrix_mse) <-  c("bayesLASSO", "SSVS")
  result_matrix_mse[1,] <- bayesLASSO_store_mse_mean
  result_matrix_mse[2,] <- SSVS_store_mse_mean
  
  result_matrix_sd <- matrix(NA, nrow = 2, ncol = 3)
  colnames(result_matrix_sd) <- c("continuous EM se", "binary EM se","treatment se")
  rownames(result_matrix_sd) <-  c("bayesLASSO", "SSVS")
  result_matrix_sd[1,] <- bayesLASSO_store_sd_mean
  result_matrix_sd[2,] <- SSVS_store_sd_mean
  
  result_matrix_bias <- matrix(NA, nrow = 2, ncol = 3)
  colnames(result_matrix_bias) <- c("bias", "bias^2","variance")
  rownames(result_matrix_bias) <-  c("bayesLASSO", "SSVS")
  result_matrix_bias[1,] <- bayesLASSO_store_bias_mean
  result_matrix_bias[2,] <- SSVS_store_bias_mean
  
  cbind(result_matrix_mse, result_matrix_sd, result_matrix_bias)
}