#' Make a (deft-approach) one-stage individual patient data meta-analysis object containing data, priors, and a JAGS model code
#'
#' This function sets up data and JAGS code that is needed to run (deft-approach) one-stage IPD-MA models in JAGS.
#' 
#' @param y outcome of the study. Can be continuous or binary. 
#' @param study vector indicating which study the patient belongs to. Please change the study names into numbers (i.e. 1, 2, 3, etc)
#' @param treat vector indicating which treatment the patient was assigned to (i.e. 1 for treatment, 0 for placebo)
#' @param X matrix of covariate values for each patient. Dimension would be number of patients x number of covariates.
#' @param response specification of the outcome type. Must specify either "normal" or "binomial".
#' @param type assumption on the treatment effect: either "random" for random effects model or "fixed" for fixed effects model. Default is "random".
#' @param mean.alpha prior mean for the study intercept
#' @param prec.alpha prior precision for the study intercept
#' @param mean.beta prior mean for the regression coefficients of the main effects of the covariates; main effects are assumed to have common effect.
#' @param prec.beta prior precision for the regression coefficients of the main effects of the covariates
#' @param mean.gamma.within prior mean for effect modifiers of within study information.
#' @param prec.gamma.within prior precision for the effect modifiers of within study information.
#' @param mean.gamma.across prior mean for the effect modifiers of across study information; effect modification is assumed to have common effect.
#' @param prec.gamma.across prior precision for the effect modifiers of across study information
#' @param mean.delta prior mean for the average treatment effect
#' @param prec.delta prior precision for the average treatment effect
#' @param hy.prior prior for the heterogeneity parameter. Supports uniform, gamma, and half normal for normal and binomial response
#' It should be a list of length 3, where first element should be the distribution (one of dunif, dgamma, dhnorm) and the next two are the parameters associated with the distribution. For example, list("dunif", 0, 5) gives uniform prior with lower bound 0 and upper bound 5 for the heterogeneity parameter.
#' @return 
#' \item{data.JAGS}{data organized in a list so that it can be used when running code in JAGS}
#' \item{code}{JAGS code that is used to run the model. Use cat(code) to see the code in a readable format}
#' \item{model.JAGS}{JAGS code in a function. This is used when running model in parallel}
#' \item{Xbar}{study specific averages of covariates}
#' @references Fisher DJ, Carpenter JR, Morris TP, et al. Meta-analytical methods to identify who benefits most from treatments: daft, deluded, or deft approach?. \emph{BMJ}. 2017;356:j573 \doi{10.1136/bmj.j573}
#' @examples
#' ds <- generate_ipdma_example(type = "continuous")
#' ipd <- with(ds, ipdma.model.deft.onestage(y = y, study = studyid, treat = treat, X = cbind(z1, z2), 
#' response = "normal", shrinkage = "none"))
#' \donttest{
#' samples <- ipd.run(ipd, pars.save = c("beta", "gamma", "delta"), n.chains = 3, n.burnin = 500, 
#' n.iter = 5000)
#' treatment.effect(ipd, samples, newpatient= c(1,0.5), reference = c(0, 0))
#' }
#' @export

ipdma.model.deft.onestage <- function(y = NULL, study = NULL, treat = NULL, X = NULL, 
                                 response = "normal", type = "random",
                                 mean.alpha = 0, prec.alpha = 0.001, mean.beta = 0, prec.beta = 0.001, 
                                 mean.gamma.within = 0, prec.gamma.within = 0.001, mean.gamma.across = 0, prec.gamma.across = 0.001, mean.delta = 0, prec.delta = 0.001,
                                 hy.prior = list("dhnorm", 0, 1)
){
  
  if(!all(grepl("^-?[0-9.]+$", study))){
    stop("Please change the study names into numbers (i.e. 1,2,3,etc)")
  }
  
  if(length(unique(treat)) > 2){
    stop("There are more than 2 different treatments specified; need to use ipdnma.model.onestage (under development)")
  }
  
  
  Xbar <- NULL
  Xbar <- matrix(NA, length(unique(study)), dim(X)[2])
  # scale with trial specific mean and sd
  for(i in 1:length(unique(study))){
    this.study <- unique(study)[i]
    Xbar[i,] <- apply(X[study == this.study,], 2, mean)
  }
  
  #JAGS data input
  data.JAGS <- 
    list(Nstudies = length(unique(study)),
         Ncovariate = dim(X)[2],
         X = X,
         Np = dim(X)[1],
         studyid = study,
         treat = treat + 1,
         y = y)
  
  
  data.JAGS$Xbar <- Xbar

  ipd <- list(y = y, study = study, treat = treat, X = X, response = response, type = type, 
              mean.alpha = mean.alpha, prec.alpha = prec.alpha, 
              mean.beta = mean.beta, prec.beta = prec.beta, mean.gamma.within = mean.gamma.within, 
              prec.gamma.within = prec.gamma.within, mean.gamma.across = mean.gamma.across, prec.gamma.across = prec.gamma.across, mean.delta = mean.delta, prec.delta = prec.delta,
              hy.prior = hy.prior)
  
  code <- ipdma.onestage.deft.rjags(ipd)
  
  code2 <- substring(code, 10)
  code2 <- sub("T(0,)", ";T(0,)", code2, fixed = T)
  eval(parse(text = paste('model.JAGS <- function() {', code2, sep='')))
  
  ipd <- list(data.JAGS = data.JAGS, code = code, model.JAGS = model.JAGS, response = response, Xbar = Xbar)
  class(ipd) <- "ipdma.onestage.deft"
  return(ipd)
}


