#' Make an one-stage individual patient data meta-analysis (IPD-MA) object containing data, priors, and a JAGS model code
#'
#' This function sets up data and JAGS code that is needed to run one-stage IPD-MA models in JAGS.
#' 
#' @param y Outcome of the study. Can be continuous or binary. 
#' @param study A vector indicating which study the patient belongs to. Please change the study names into numbers (i.e. 1, 2, 3, etc)
#' @param treat A vector indicating which treatment the patient was assigned to (i.e. 1 for treatment, 0 for placebo)
#' @param X A matrix of covariate values for each patient. Dimension would be number of patients x number of covariates.
#' @param response Specification of the outcome type. Must specify either "normal" or "binomial".
#' @param type Assumption on the treatment effect: either "random" for random effects model or "fixed" for fixed effects model. Default is "random".
#' @param shrinkage Shrinkage method applied to the effect modifiers. "none" correspond to no shrinkage.
#' "laplace" corresponds to a adaptive shrinkage with a Laplacian prior (ie often known as Bayesian LASSO).
#' "SSVS" corresponds to the Stochastic Search Variable Selection method. SSVS is not strictly a shrinkage method, but pulls the estimated coefficient toward zero through variable selection in each iteration of the MCMC. 
#' See O'hara et al (2009) for more details.
#' @param scale Indicator for scaling the covariates by the overall average; default is TRUE.
#' @param mean.alpha Prior mean for the study intercept
#' @param prec.alpha Prior precision for the study intercept
#' @param mean.beta Prior mean for the regression coefficients of the main effects of the covariates; main effects are assumed to have common effect.
#' @param prec.beta Prior precision for the regression coefficients of the main effects of the covariates
#' @param mean.gamma Prior mean for the effect modifiers. This parameter is not used if penalization is placed on effect modifiers.
#' @param prec.gamma Prior precision for the effect modifiers. This parameter is not used if penalization is placed on effect modifiers.
#' @param mean.delta Prior mean for the average treatment effect
#' @param prec.delta Prior precision for the average treatment effect
#' @param hy.prior Prior for the heterogeneity parameter. Supports uniform, gamma, and half normal for normal and binomial response
#' It should be a list of length 3, where first element should be the distribution (one of dunif, dgamma, dhnorm) and the next two are the parameters associated with the distribution. For example, list("dunif", 0, 5) give uniform prior with lower bound 0 and upper bound 5 for the heterogeneity parameter.
#' @param lambda.prior (Only for shrinkage = "laplace") Two options for laplace shrinkage. We can put a gamma prior on the lambda (i.e. list("dgamma",2,0.1)) or put a uniform prior on the inverse of lambda (i.e. list("dunif",0,5))
#' @param p.ind (Only for shrinkage = "SSVS") Prior probability of including each of the effect modifiers. Length should be same as the total length of the covariates.
#' @param g (Only for shrinkage = "SSVS") Multiplier for the precision of spike. Default is g = 1000.
#' @param hy.prior.eta (Only for shrinkage = "SSVS") Standard deviation of the slab prior. Currently only support uniform distribution. Default is list("dunif", 0, 5)
#' @return 
#' \item{data.JAGS}{Data organized in a list so that it can be used when running code in JAGS}
#' \item{code}{JAGS code that is used to run the model. Use cat(code) to see the code in a readable format}
#' \item{model.JAGS}{JAGS code in a function. This is used when running model in parallel}
#' \item{scale.mean}{Mean used in scaling covariates}
#' \item{scale.sd}{Standard deviation used in scaling covariates}
#' @references Seo M, White IR, Furukawa TA, et al. Comparing methods for estimating patient-specific treatment effects in individual patient data meta-analysis. \emph{Stat Med}. 2021;40(6):1553-1573. \doi{10.1002/sim.8859}
#' @examples
#' ds <- generate_ipdma_example(type = "continuous")
#' ipd <- with(ds, ipdma.model.onestage(y = y, study = studyid, treat = treat, X = cbind(z1, z2), 
#' response = "normal", shrinkage = "none"))
#' @export

ipdma.model.onestage <- function(y = NULL, study = NULL, treat = NULL, X = NULL, 
                      response = "normal", type = "random", shrinkage = "none", scale = TRUE,
                      mean.alpha = 0, prec.alpha = 0.001, mean.beta = 0, prec.beta = 0.001, 
                      mean.gamma = 0, prec.gamma = 0.001, mean.delta = 0, prec.delta = 0.001,
                      hy.prior = list("dhnorm", 0, 1), lambda.prior = NULL, p.ind = NULL, g = NULL, hy.prior.eta = NULL
                      ){

  if(!all(grepl("^-?[0-9.]+$", study))){
    stop("Please change the study names into numbers (i.e. 1,2,3,etc)")
  }
  
  if(length(unique(treat)) > 2){
    stop("There are more than 2 different treatments specified; need to use ipdnma.model.onestage (under development)")
  }

  if(shrinkage== "none" & (!is.null(lambda.prior) || !is.null(p.ind) || !is.null(g) || !is.null(hy.prior.eta))){
    stop("Shrinkage is set to none but have specified prior for shrinkage parameters")
  }
  
  
  #center the covariates
  scale_mean <- scale_sd <- NULL
  if(scale == TRUE){
    scale_mean <- apply(X, 2, mean)
    scale_sd <- apply(X, 2, sd)
    X <- apply(X, 2, scale) 
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
  
  # default prior assignment
  if(is.null(lambda.prior)) lambda.prior <- list("dunif", 0, 5)
  if(is.null(p.ind)) p.ind <- rep(0.5, dim(X)[2])
  if(is.null(g)) g <- 1000
  if(is.null(hy.prior.eta)) hy.prior.eta <- list("dunif", 0, 5)
  
  if(shrinkage == "SSVS"){
    data.JAGS$p.ind <- p.ind
  }

  ipd <- list(y = y, study = study, treat = treat, X = X, response = response, type = type, 
              shrinkage = shrinkage, mean.alpha = mean.alpha, prec.alpha = prec.alpha, 
              mean.beta = mean.beta, prec.beta = prec.beta, mean.gamma = mean.gamma, 
              prec.gamma = prec.gamma, mean.delta = mean.delta, prec.delta = prec.delta,
              hy.prior = hy.prior, lambda.prior = lambda.prior, p.ind = p.ind, g = g, hy.prior.eta = hy.prior.eta)
  
  code <- ipdma.onestage.rjags(ipd)
  
  code2 <- substring(code, 10)
  code2 <- sub("T(0,)", ";T(0,)", code2, fixed = T)
  eval(parse(text = paste('model.JAGS <- function() {', code2, sep='')))

  ipd <- list(data.JAGS = data.JAGS, code = code, model.JAGS = model.JAGS, response = response, scale_mean = scale_mean, scale_sd = scale_sd)
  class(ipd) <- "ipdma.onestage"
  return(ipd)
}


