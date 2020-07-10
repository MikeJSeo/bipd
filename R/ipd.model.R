#' Make an ipd object containing data, priors, and a JAGS model file
#'
#' This function sets up data and JAGS code that is needed to run ipd models in JAGS
#' 
#' @param y outcome of the study. Can be continuous or binary. 
#' @param study a vector indicating which study the patient belongs to. Please change the study names into numbers (i.e. 1,2,3,etc)
#' @param treat a vector indicating which treatment the patient was assigned to (ie 1 for treatment, 0 for placebo)
#' @param X a matrix of covariate values for each patient. Dimension would be number of patients x number of covariates.
#' @param response Specification of the outcomes type. Must specify either "normal" or "binomial"
#' @param type Assumption on the treatment effect: either "random" for random effects model or "fixed" for fixed effects model. Default is "random".
#' @param model "onestage" for one-stage model; two-stage are planned to be implemented
#' @param shrinkage shrinkage method applied to the effect modifiers. "none" correspond to no shrinkage.
#' "laplace" corresponds to a adaptive shrinkage with a Laplacian prior (ie often known as Bayesian LASSO).
#' "SSVS" corresponds to the Search Variable Selection method. SSVS is not strictly a shrinkage method, 
#' but pulls the estimated coefficient toward zero through variable selection in each iteration of the MCMC. 
#' See O'hara et al (2009) for more reference.
#' @param mean.a Prior mean for the study intercept
#' @param prec.a Prior precision for the study intercept
#' @param mean.beta Prior mean for the regression coefficients of the main effects of the covariates 
#' @param prec.beta Prior precision for the regression coefficients of the main effects of the covariates
#' @param mean.gamma Prior mean for the effect modifiers. This parameter is not used if penalization is placed on effect modifiers. 
#' @param prec.gamma Prior precision for the effect modifiers.
#' @param mean.delta Prior mean for the treatment effect
#' @param prec.delta Prior precision for the treatment effect
#' @param hy.prior Prior for the heterogeneity parameter. Supports uniform, gamma, and half normal for normal and binomial response
#' It should be a list of length 3, where first element should be the distribution (one of dunif, dgamma, dhnorm) and the next two are the parameters associated with the distribution. For example, list("dunif", 0, 5) give uniform prior with lower bound 0 and upper bound 5 for the heterogeneity parameter.
#' @param lambda.prior (Only for shrinkage = "laplace") Two options for laplace shrinkage. We can put a gamma prior on the lambda (i.e. list("dgamma",2,0.1)) or put a uniform prior on the inverse of lambda (i.e. list("dunif",0,5))
#' @param p.ind (Only for shrinkage = "SSVS") Prior probability of including each of the effect modifiers. Length should be same as the total length of the covariates.
#' @return 
#' \item{data}{Data organized in a list so that it can be used when running code in JAGS}
#' \item{code}{JAGS code that is used to run the model. Use cat(code) to see the code in a nicer format.}
#'
#' @export

ipd.model <- function(y = NULL, study = NULL, treat = NULL, X = NULL, 
                      response = "normal", type = "random", model = "onestage", shrinkage = "none",
                      mean.a = 0, prec.a = 0.001, mean.beta = 0, prec.beta = 0.001, 
                      mean.gamma = 0, prec.gamma = 0.001, mean.delta = 0, prec.delta = 0.001,
                      hy.prior = list("dhnorm", 0, 1), lambda.prior = NULL, p.ind = NULL
                      ){

  data <- 
    list(Nstudies = length(unique(study)),
         Ncovariate = dim(X)[2],
         X = X,
         Np = dim(X)[1],
         studyid = study,
         treat = treat + 1,
         y = y)
  
  if(is.null(lambda.prior)) lambda.prior <- list("dunif", 0, 5)
  if(is.null(p.ind)) p.ind <- rep(0.5, dim(X)[2])
  
  if(shrinkage == "SSVS"){
    data$p.ind <- p.ind
  }
         
  ipd <- list(y = y, study = study, treat = treat, X = X, response = response, type = type, 
              model = model, shrinkage = shrinkage, mean.a = mean.a, prec.a = prec.a, 
              mean.beta = mean.beta, prec.beta = prec.beta, mean.gamma = mean.gamma, 
              prec.gamma = prec.gamma, mean.delta = mean.delta, prec.delta = prec.delta,
              hy.prior = hy.prior, lambda.prior = lambda.prior, p.ind = p.ind)

  code <- ipd.rjags(ipd)

  list(data = data, code = code)
}


