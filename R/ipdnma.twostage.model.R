#' Make a first part of the two-stage individual patient data network meta-analysis object containing data, priors, and a JAGS model code
#'
#' This is the first part of the two-stage IPD-NMA model that sets up data and JAGS code needed. 
#' Each study is fitted using this function.
#' 
#' @param y outcome of the study. Can be continuous or binary. 
#' @param treat vector indicating which treatment the patient was assigned to. Since this is a network meta-analysis and there would be more than 2 treatments.
#' However this doesn't mean all studies have more than 2 treatments and they can only have two. One thing to be careful is treatment names should be same for all studies.
#' Treatment that is assigned 1 would be the baseline treatment.
#' @param X matrix of covariate values for each patient. Dimension would be number of patients x number of covariates.
#' @param response specification of the outcome type. Must specify either "normal" or "binomial".
#' @param shrinkage shrinkage method applied to the effect modifiers. "none" correspond to no shrinkage.
#' "laplace" corresponds to a adaptive shrinkage with a Laplacian prior (ie often known as Bayesian LASSO).
#' "SSVS" corresponds to the Stochastic Search Variable Selection method. SSVS is not strictly a shrinkage method, but pulls the estimated coefficient toward zero through variable selection in each iteration of the MCMC. 
#' See O'hara et al (2009) for more details.
#' @param scale indicator for scaling the covariates by the overall average; default is TRUE.
#' @param mean.a prior mean for the study intercept
#' @param prec.a prior precision for the study intercept
#' @param mean.b prior mean for the regression coefficients of the main effects of the covariates; main effects are assumed to have common effect.
#' @param prec.b prior precision for the regression coefficients of the main effects of the covariates
#' @param mean.c prior mean for the effect modifiers. This parameter is not used if penalization is placed on effect modifiers.
#' @param prec.c prior precision for the effect modifiers. This parameter is not used if penalization is placed on effect modifiers.
#' @param mean.d prior mean for the average treatment effect
#' @param prec.d prior precision for the average treatment effect
#' @param lambda.prior (only for shrinkage = "laplace") two options for laplace shrinkage. We can put a gamma prior on the lambda (i.e. list("dgamma",2,0.1)) or put a uniform prior on the inverse of lambda (i.e. list("dunif",0,5))
#' @param p.ind (only for shrinkage = "SSVS") prior probability of including each of the effect modifiers. Length should be same as the total length of the covariates.
#' @param g (only for shrinkage = "SSVS") multiplier for the precision of spike. Default is g = 1000.
#' @param hy.prior.eta (only for shrinkage = "SSVS") standard deviation of the slab prior. Currently only support uniform distribution. Default is list("dunif", 0, 5)
#' @return 
#' \item{data.JAGS}{data organized in a list so that it can be used when running code in JAGS}
#' \item{code}{JAGS code that is used to run the model. Use cat(code) to see the code in a readable format}
#' \item{model.JAGS}{JAGS code in a function. This is used when running model in parallel}
#' \item{scale.mean}{mean used in scaling covariates}
#' \item{scale.sd}{standard deviation used in scaling covariates}
#' \item{treat.original}{original treatment names from user}
#' \item{treat.relabeled}{treatment relabeled for the purpose of modelling}
#' @references Dias S, Sutton AJ, Ades AE, et al. A Generalized Linear Modeling Framework for Pairwise and Network Meta-analysis of Randomized Controlled Trials. \emph{Medical Decision Making}. 2013;33(5):607-617. \doi{10.1177/0272989X12458724}
#' @references O'Hara RB, Sillanpaa MJ. A review of Bayesian variable selection methods: what, how and which. \emph{Bayesian Anal}. 2009;4(1):85-117. \doi{10.1214/09-BA403}
#' @references Seo M, White IR, Furukawa TA, et al. Comparing methods for estimating patient-specific treatment effects in individual patient data meta-analysis. \emph{Stat Med}. 2021;40(6):1553-1573. \doi{10.1002/sim.8859}
#' @examples
#' ds <- generate_ipdnma_example(type = "continuous")
#' ipd <- with(ds, ipdnma.model.twostage.first(y = y, treat = treat, X = cbind(z1, z2), 
#' response = "normal", shrinkage = "none"))
#' \donttest{
#' samples <- ipd.run(ipd)
#' }
#' @export


ipdnma.model.twostage.first <- function(y = NULL, treat = NULL, X = NULL, 
                                  response = "normal", shrinkage = "none", scale = TRUE,
                                  mean.a = 0, prec.a = 0.001, mean.b = 0, prec.b = 0.001, 
                                  mean.c = 0, prec.c = 0.001, mean.d = 0, prec.d = 0.001,
                                  lambda.prior = NULL, p.ind = NULL, g = NULL, hy.prior.eta = NULL){
  
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
  
  #relabel treat for fitting model purposes
  #however original treatment names are returned 
  treat.relabeled <- as.numeric(as.factor(treat))
  
  #JAGS data input
  data.JAGS <- 
    list(Ncovariate = dim(X)[2],
         Ntreat = length(unique(treat)),
         X = X,
         Np = dim(X)[1],
         treat = treat.relabeled,
         y = y)
  
  # default prior assignment
  if(is.null(lambda.prior)) lambda.prior <- list("dunif", 0, 5)
  if(is.null(p.ind)) p.ind <- rep(0.5, dim(X)[2])
  if(is.null(g)) g <- 1000
  if(is.null(hy.prior.eta)) hy.prior.eta <- list("dunif", 0, 5)
  
  if(shrinkage == "SSVS"){
    data.JAGS$p.ind <- p.ind
  }
  
  ipd <- list(y = y, treat = treat.relabeled, X = X, response = response, 
              shrinkage = shrinkage, mean.a = mean.a, prec.a = prec.a, 
              mean.b = mean.b, prec.b = prec.b, mean.c = mean.c, 
              prec.c = prec.c, mean.d = mean.d, prec.d = prec.d,
              lambda.prior = lambda.prior, p.ind = p.ind, g = g, hy.prior.eta = hy.prior.eta)
  
  code <- ipdnma.twostage.first.rjags(ipd)
  
  code2 <- substring(code, 10)
  code2 <- sub("T(0,)", ";T(0,)", code2, fixed = T)
  model.JAGS <- NULL
  eval(parse(text = paste('model.JAGS <- function() {', code2, sep='')))
  
  ipd <- list(data.JAGS = data.JAGS, code = code, model.JAGS = model.JAGS, response = response, scale_mean = scale_mean, scale_sd = scale_sd, treat.original = treat, treat.relabeled = treat.relabeled)
  class(ipd) <- "ipdnma.twostage.first"
  return(ipd)
}

#' Make the second part of the two-stage individual patient data network meta-analysis object containing data, priors, and a JAGS model code
#'
#' This is the second part of the two-stage IPD-NMA model that sets up data and JAGS code needed. 
#' This model aggregates results from first stage model.
#' 
#' @param ipd_all list of first part of the two-stage IPD-NMA model object
#' @param samples_all list of samples from running first part of the two-stage IPD-NMA model
#' @param mean.a prior mean for the study intercept
#' @param prec.a prior precision for the study intercept
#' @param mean.b prior mean for the regression coefficients of the main effects of the covariates; main effects are assumed to have common effect.
#' @param prec.b prior precision for the regression coefficients of the main effects of the covariates
#' @param mean.c prior mean for the effect modifiers. This parameter is not used if penalization is placed on effect modifiers.
#' @param prec.c prior precision for the effect modifiers. This parameter is not used if penalization is placed on effect modifiers.
#' @param mean.d prior mean for the average treatment effect
#' @param prec.d prior precision for the average treatment effect
#' @return 
#' \item{data.JAGS}{data organized in a list so that it can be used when running code in JAGS}
#' \item{code}{JAGS code that is used to run the model. Use cat(code) to see the code in a readable format}
#' \item{model.JAGS}{JAGS code in a function. This is used when running model in parallel}
#' @references Dias S, Sutton AJ, Ades AE, et al. A Generalized Linear Modeling Framework for Pairwise and Network Meta-analysis of Randomized Controlled Trials. \emph{Medical Decision Making}. 2013;33(5):607-617. \doi{10.1177/0272989X12458724}
#' @references Seo M, White IR, Furukawa TA, et al. Comparing methods for estimating patient-specific treatment effects in individual patient data meta-analysis. \emph{Stat Med}. 2021;40(6):1553-1573. \doi{10.1002/sim.8859}
#' @examples
#' ds <- generate_ipdnma_example(type = "continuous")
#' #TODO
#' @export

ipdnma.model.twostage.second <- function(ipd_all = NULL, samples_all = NULL,
                                        mean.a = 0, prec.a = 0.001, mean.b = 0, prec.b = 0.001, 
                                        mean.c = 0, prec.c = 0.001, mean.d = 0, prec.d = 0.001,
                                        lambda.prior = NULL, p.ind = NULL, g = NULL, hy.prior.eta = NULL){
  
  
  
}



