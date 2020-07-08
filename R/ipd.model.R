#' Make an ipd object containing data, priors, and a JAGS model file
#'
#' This function sets up data and JAGS code that is needed to run ipd models in JAGS
#' 
#' @param y outcome of the study. Can be continuous or binary. 
#' @param study a vector indicating which study the patient belongs to. Please change the study names into numbers (i.e. 1,2,3,etc)
#' @param treat a vector indicating which treatment the patient was assigned to. Please change the study names into numbers (i.e. 1,2,3,etc)
#' @param X a matrix of covariate values for each patient. Dimension would be number of patients x number of covariates.
#' @param response "continuous" for continuous outcome and "binary" for binary outcome
#' @param model "onestage" for one-stage model; two-stage are planned to be implemented
#' @param shrinkage shrinkage method applied to the effect modifiers. "none" correspond to no shrinkage.
#' "laplace" corresponds to a adaptive shrinkage with a Laplacian prior (ie often known as Bayesian LASSO).
#' "SSVS" corresponds to the Search Variable Selection method. SSVS is not strictly a shrinkage method, 
#' but shrinks the estimated coefficient through variable selection in each iteration of the MCMC. 
#' See O'hara et al (2009) for more reference.
#' @return 
#' \item{data}{Data organized in a list so that it can be used when running code in JAGS}
#' \item{code}{JAGS code that is used to run the model. Use cat(code) to see the code in a nicer format.}
#'
#' @export

ipd.model <- function(y = NULL, study = NULL, treat = NULL, X = NULL, 
                               response = "continuous", model = "onestage", shrinkage = "none"){

  
  data <- 
    list(Nstudies = length(unique(study)),
         Ncovariates = dim(X)[2],
         X = X,
         Np = dim(X)[1],
         studyid = study,
         treat = treat,
         y = y)
         
  ipd <- list(y = y, study = study, treat = treat, X = X, response = response, model = model)
  
  code <- ipd.rjags(ipd)


  list(data = data, code = code)
}


