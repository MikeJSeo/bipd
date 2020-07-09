#' Make an ipd object containing data, priors, and a JAGS model file
#'
#' This function sets up data and JAGS code that is needed to run ipd models in JAGS
#' 
#' @param y outcome of the study. Can be continuous or binary. 
#' @param study a vector indicating which study the patient belongs to. Please change the study names into numbers (i.e. 1,2,3,etc)
#' @param treat a vector indicating which treatment the patient was assigned to. Please change the study names into numbers (i.e. 1,2,3,etc)
#' @param X a matrix of covariate values for each patient. Dimension would be number of patients x number of covariates.
#' @param response "continuous" for continuous outcome and "binary" for binary outcome
#' @param type Assumption on the treatment effect: either "random" for random effects model or "fixed" for fixed effects model. Default is "random".
#' @param model "onestage" for one-stage model; two-stage are planned to be implemented
#' @param shrinkage shrinkage method applied to the effect modifiers. "none" correspond to no shrinkage.
#' "laplace" corresponds to a adaptive shrinkage with a Laplacian prior (ie often known as Bayesian LASSO).
#' "SSVS" corresponds to the Search Variable Selection method. SSVS is not strictly a shrinkage method, 
#' but pulls the estimated coefficient toward zero through variable selection in each iteration of the MCMC. 
#' See O'hara et al (2009) for more reference.
#' @param mean.alpha Prior mean for the study intercept
#' @param prec.alpha Prior precision for the study intercept
#' @param mean.beta Prior mean for the regression coefficients of the main effects of the covariates 
#' @param prec.beta Prior mean for the regression coefficients of the main effects of the covariates
#' @param mean.delta Prior mean for the treatment effect
#' @param prec.delta Prior precision for the treatment effect
#' @param hy.prior Prior for the heterogeneity parameter. Supports uniform, gamma, and half normal for normal and binomial response
#' It should be a list of length 3, where first element should be the distribution (one of dunif, dgamma, dhnorm) and the next two are the parameters associated with the distribution. For example, list("dunif", 0, 5) give uniform prior with lower bound 0 and upper bound 5 for the heterogeneity parameter.
#' @return 
#' \item{data}{Data organized in a list so that it can be used when running code in JAGS}
#' \item{code}{JAGS code that is used to run the model. Use cat(code) to see the code in a nicer format.}
#'
#' @export

ipd.model <- function(y = NULL, study = NULL, treat = NULL, X = NULL, 
                      response = "continuous", type = "random", model = "onestage", shrinkage = "none"){

  
  data <- 
    list(Nstudies = length(unique(study)),
         Ncovariates = dim(X)[2],
         X = X,
         Np = dim(X)[1],
         studyid = study,
         treat = treat,
         y = y)
         
  ipd <- list(y = y, study = study, treat = treat, X = X, response = response, type = type, model = model, shrinkage = shrinkage)
  
  code <- ipd.rjags(ipd)

  list(data = data, code = code)
}


