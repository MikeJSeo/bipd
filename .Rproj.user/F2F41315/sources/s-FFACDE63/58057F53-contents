set.seed(42) #the answer to the question of life, the universe, and everything

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
head(ds)


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
ds2 <- as.data.frame(array(NA, dim=c(N*6, 5)))
colnames(ds2) <- c("studyid", "treat", "w1","w2", "x")

expit <- function(x){
  exp(x)/(1+exp(x))
}

for (i in 1:N.trials) {
  treat <- rbinom(N,1,0.5)
  w1 <- rnorm(N, mean=0, sd=1) 
  w2 <- rnorm(N, mean=0, sd=1)
  linearpredictor <- alpha[i] + delta[i]*treat + beta1[i]*w1 + beta2[i]*w2 + gamma1[i]*w1*treat + gamma2[i]*w2*treat
  x <- rbinom(N, 1, expit(linearpredictor))
  ds2[(((i-1)*N)+1):(i*N), ] <- cbind(studyid[i], treat, w1, w2, x)
}
ds2$studyid <- as.factor(ds2$studyid)
head(ds2)


## GLMM
# continuous outcome
library(lme4) #for fitting glmm
m1 <- lmer(y ~ studyid + (z1+z2)*treat + (-1 + treat|studyid), data = ds)
summary(m1)

# binary outcome
m2 <- glmer(x ~ studyid + (w1+w2)*treat + (-1 + treat|studyid), data = ds2, family = binomial)
summary(m2)
              

## STEP
# continuous outcome
m3 <- glm(y ~ studyid + (z1+z2)*treat, data = ds) #glm model without mixed effects 
s1 <- step(m3, scope=list(lower = ~ z1+z2+treat), direction = "both")
summary(s1)

# binary outcome
m4 <- glm(x ~ studyid + (w1+w2)*treat, family = binomial(link = "logit"), data = ds2)
s2 <- step(m4, scope=list(lower = ~ w1+w2+treat), direction = "both") 
summary(s2)


## LASSO
# continuous outcome
library(glmnet)
p.fac <- c(rep(0, 5), rep(0, 2), 0, rep(1,2)) #No shrinkage for treatment effect and baseline effect; only on effect modifier denoted with 1's at the end
lambdas <- 10^seq(3, -3, by = -.1) #manually specify lambda value to cross validate

data_glmnet <- model.matrix(y~ studyid + (z1+z2)*treat, data = ds)
data_glmnet <- data_glmnet[,-1]
data_glmnet <- cbind(y = ds$y, data_glmnet = data_glmnet)
cvfit = cv.glmnet(as.matrix(data_glmnet[,-1]), as.matrix(data_glmnet[,1]), penalty.factor = p.fac, family = "gaussian", type.measure = "deviance", lambda = lambdas)
coef(cvfit, s = "lambda.min")

# binary outcome
data_glmnet <- model.matrix(x~ studyid + (w1+w2)*treat, data = ds2)
data_glmnet <- data_glmnet[,-1]
data_glmnet <- cbind(x = ds2$x, data_glmnet = data_glmnet)
cvfit = cv.glmnet(as.matrix(data_glmnet[,-1]), as.matrix(data_glmnet[,1]), penalty.factor = p.fac, family = "binomial", type.measure = "deviance", lambda = lambdas)
coef(cvfit, s = "lambda.min")

## RIDGE
# continuous outcome
data_glmnet <- model.matrix(y~ studyid + (z1+z2)*treat, data = ds)
data_glmnet <- data_glmnet[,-1]
data_glmnet <- cbind(y = ds$y, data_glmnet = data_glmnet)
cvfit = cv.glmnet(as.matrix(data_glmnet[,-1]), as.matrix(data_glmnet[,1]), penalty.factor = p.fac, family = "gaussian", alpha = 0, type.measure = "deviance", lambda = lambdas)
coef(cvfit, s = "lambda.min")

# binary outcome
data_glmnet <- model.matrix(x~ studyid + (w1+w2)*treat, data = ds2)
data_glmnet <- data_glmnet[,-1]
data_glmnet <- cbind(x = ds2$x, data_glmnet = data_glmnet)
cvfit = cv.glmnet(as.matrix(data_glmnet[,-1]), as.matrix(data_glmnet[,1]), penalty.factor = p.fac, family = "binomial", alpha = 0, type.measure = "deviance", lambda = lambdas)
coef(cvfit, s = "lambda.min")


## BAYES-LASSO
library(R2jags)# for Bayesian models
library(coda) # for diagnostic plots

# continuous outcome
X <- cbind(ds$z1, ds$z2)
data_jags <- with(ds, {
  list(Nstudies = length(unique(studyid)),
       Ncovariate = 2,
       X = X,
       Np = dim(X)[1],
       studyid = studyid,
       treat = treat + 1,
       y = y
)})

modelBayesLasso <- function(){
  
  ########## IPD-MA model
  for (i in 1:Np){
    y[i] ~ dnorm(mu[i], sigma)  
    mu[i]<-a[studyid[i]] + inprod(beta[], X[i,]) +
      (1 - equals(treat[i],1)) * inprod(gamma[], X[i,]) + d[studyid[i],treat[i]] 
  }
  sigma ~ dgamma(0.001, 0.001)
  
  #####treatment effect
  for(j in 1:Nstudies){
    d[j,1] <- 0
    d[j,2] ~ dnorm(delta[2], tau)
  }
  
  ## prior distribution for heterogeneity of the treatment effect
  tau <- pow(sd, -2)
  sd ~ dnorm(0,1);T(0,)
  
  ## prior distribution for average treatment effect
  delta[1] <- 0
  delta[2] ~ dnorm(0, 0.001)
  
  ## prior distribution for baseline effect
  for (j in 1:Nstudies){
    a[j] ~ dnorm(0, 0.001)	
  }
  
  tt <- lambda * sigma  # uses conditional prior
  lambda ~ dgamma(2, 0.1) # this diffuse prior has mean 20, variance 200 
  
  for(k in 1:Ncovariate){
    beta[k] ~ dnorm(0, 0.001)
  }
  
  for(k in 1:Ncovariate){
    gamma[k] ~ ddexp(0, tt) # penalization is only on gamma
  }
}

model_bayesLASSO<- jags(data = data_jags, inits = NULL,
                        parameters.to.save =  c("gamma", "delta", "beta", "sd"),
                        n.chains = 3, n.iter = 10000,
                        n.burnin = 1000,DIC=F,
                        model.file = modelBayesLasso)

model_bayesLASSO
jagsfit.mcmc <- as.mcmc(model_bayesLASSO)
jagsfit.mcmc <- jagsfit.mcmc[,-3] #remove delta[1] which is 0
plot(jagsfit.mcmc) #traceplot and posterior of parameters
gelman.plot(jagsfit.mcmc) #gelman diagnostic plot


# dichotomous outcome
X <- cbind(ds2$w1, ds2$w2)
data_jags_binary <- with(ds2, {
  list(Nstudies = length(unique(studyid)),
       Ncovariate = 2,
       X = X,
       Np = dim(X)[1],
       studyid = studyid,
       treat = treat + 1,
       y = x
)})


modelBayesLasso_binary <- function(){
  
  ########## IPD-MA model
  for (i in 1:Np){
    y[i] ~ dbern(p[i])  
    logit(p[i]) <- a[studyid[i]] + inprod(beta[], X[i,]) +
      (1 - equals(treat[i],1)) * inprod(gamma[], X[i,]) + d[studyid[i],treat[i]] 
  }
  
  #####treatment effect
  for(j in 1:Nstudies){
    d[j,1] <- 0
    d[j,2] ~ dnorm(delta[2], tau)
  }
  
  ## prior distribution for heterogeneity of the treatment effect
  tau <- pow(sd, -2)
  sd ~ dnorm(0,1);T(0,)
  
  ## prior distribution for average treatment effect
  delta[1] <- 0
  delta[2] ~ dnorm(0, 0.001)
  
  ## prior distribution for baseline effect
  for (j in 1:Nstudies){
    a[j] ~ dnorm(0, 0.001)	
  }
  
  ## prior distribution for penalization parameter
  #lambda ~ dgamma(0.001, 0.001)
  lambda ~ dgamma(2, 0.1)
  
  for(k in 1:Ncovariate){
    beta[k] ~ dnorm(0, 0.001)
  }
  
  for(k in 1:Ncovariate){
    gamma[k] ~ ddexp(0, lambda) 
  }
}

model_bayesLASSO_binary <- jags(data = data_jags_binary, inits = NULL,
                        parameters.to.save =  c("gamma", "delta", "beta", "sd"),
                        n.chains = 3, n.iter = 10000,
                        n.burnin = 1000,DIC=F,
                        model.file = modelBayesLasso_binary)

model_bayesLASSO_binary
jagsfit.mcmc2 <- as.mcmc(model_bayesLASSO_binary)
jagsfit.mcmc2 <- jagsfit.mcmc2[,-3] #remove delta[1] which is 0
plot(jagsfit.mcmc2) #traceplot and posterior of parameters
gelman.plot(jagsfit.mcmc2) #gelman diagnostic plot



## SSVS

# continuous outcome
modelSSVS <- function(){
  
  ########## IPD-MA model
  for (i in 1:Np){
    y[i] ~ dnorm(mu[i], sigma)  
    mu[i]<- a[studyid[i]] + inprod(beta[], X[i,]) +
      (1 - equals(treat[i],1)) * inprod(gamma[], X[i,]) + d[studyid[i],treat[i]] 
  }
  sigma ~ dgamma(0.001, 0.001)
  
  #####treatment effect
  for(j in 1:Nstudies){
    d[j,1] <- 0
    d[j,2] ~ dnorm(delta[2], tau)
  }
  
  ## prior distribution for heterogeneity of the treatment effect
  tau <- pow(sd, -2)
  sd ~ dnorm(0,1);T(0,)
  
  ## prior distribution for average treatment effect
  delta[1] <- 0
  delta[2] ~ dnorm(0, 0.001)
  
  ## prior distribution for baseline effect
  for (j in 1:Nstudies){
    a[j] ~ dnorm(0, 0.001)	
  }
  
  for(k in 1:Ncovariate){
    beta[k] ~ dnorm(0, 0.001)
  }
  
  for(k in 1:Ncovariate){
    IndA[k] ~ dcat(Pind[])
    Ind[k] <- IndA[k] - 1
    gamma[k] ~ dnorm(0, tauCov[IndA[k]])
  }
  
  zeta <- pow(eta, -2)
  eta ~ dunif(0, 5)
  tauCov[1] <- zeta
  tauCov[2] <- zeta * 0.01  # g = 100
  
  Pind[1] <- 0.5 #P(I_j=1)= 0.5
  Pind[2] <- 0.5 
}

model_SSVS <- jags(data = data_jags, inits = NULL,
                        parameters.to.save =  c("gamma", "delta", "beta", "sd", "Ind"),
                        n.chains = 3, n.iter = 10000,
                        n.burnin = 1000,DIC=F,
                        model.file = modelSSVS)

model_SSVS
jagsfit.mcmc3 <- as.mcmc(model_SSVS)
jagsfit.mcmc3 <- jagsfit.mcmc3[,-3] #remove delta[1] which is 0
plot(jagsfit.mcmc3) #traceplot and posterior of parameters
gelman.plot(jagsfit.mcmc3) #gelman diagnostic plot


# dichotomous outcome
modelSSVS_binary <- function(){
  
  ########## IPD-MA model
  for (i in 1:Np){
    y[i] ~ dbern(p[i])
    logit(p[i]) <- a[studyid[i]] + inprod(beta[], X[i,]) +
      (1 - equals(treat[i],1)) * inprod(gamma[], X[i,]) + d[studyid[i],treat[i]] 
  }
  
  #####treatment effect
  for(j in 1:Nstudies){
    d[j,1] <- 0
    d[j,2] ~ dnorm(delta[2], tau)
  }
  
  ## prior distribution for heterogeneity of the treatment effect
  tau <- pow(sd, -2)
  sd ~ dnorm(0,1);T(0,)
  
  ## prior distribution for average treatment effect
  delta[1] <- 0
  delta[2] ~ dnorm(0, 0.001)
  
  ## prior distribution for baseline effect
  for (j in 1:Nstudies){
    a[j] ~ dnorm(0, 0.001)	
  }
  
  for(k in 1:Ncovariate){
    beta[k] ~ dnorm(0, 0.001)
  }
  
  for(k in 1:Ncovariate){
    IndA[k] ~ dcat(Pind[])
    Ind[k] <- IndA[k] - 1
    gamma[k] ~ dnorm(0, tauCov[IndA[k]])
  }
  
  zeta <- pow(eta, -2)
  eta ~ dunif(0, 5)
  tauCov[1] <- zeta
  tauCov[2] <- zeta * 0.01  # g = 100
  
  Pind[1] <- 0.5 #P(I_j=1)= 0.5
  Pind[2] <- 0.5 
  
}

model_SSVS_binary <- jags(data = data_jags_binary, inits = NULL,
                   parameters.to.save =  c("gamma", "delta", "beta", "sd", "Ind"),
                   n.chains = 3, n.iter = 10000,
                   n.burnin = 1000,DIC=F,
                   model.file = modelSSVS_binary)

model_SSVS_binary
jagsfit.mcmc4 <- as.mcmc(model_SSVS_binary)
jagsfit.mcmc4 <- jagsfit.mcmc4[,-3] #remove delta[1] which is 0
plot(jagsfit.mcmc4) #traceplot and posterior of parameters
gelman.plot(jagsfit.mcmc4) #gelman diagnostic plot

