% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ipdma.onestage.model.R
\name{ipdma.model.onestage}
\alias{ipdma.model.onestage}
\title{Make an one-stage individual patient data meta-analysis object containing data, priors, and a JAGS model code}
\usage{
ipdma.model.onestage(
  y = NULL,
  study = NULL,
  treat = NULL,
  X = NULL,
  response = "normal",
  type = "random",
  shrinkage = "none",
  scale = TRUE,
  mean.alpha = 0,
  prec.alpha = 0.001,
  mean.beta = 0,
  prec.beta = 0.001,
  mean.gamma = 0,
  prec.gamma = 0.001,
  mean.delta = 0,
  prec.delta = 0.001,
  hy.prior = list("dhnorm", 0, 1),
  lambda.prior = NULL,
  p.ind = NULL,
  g = NULL,
  hy.prior.eta = NULL
)
}
\arguments{
\item{y}{outcome of the study. Can be continuous or binary.}

\item{study}{vector indicating which study the patient belongs to. Please change the study names into numbers (i.e. 1, 2, 3, etc)}

\item{treat}{vector indicating which treatment the patient was assigned to (i.e. 1 for treatment, 0 for placebo)}

\item{X}{matrix of covariate values for each patient. Dimension would be number of patients x number of covariates.}

\item{response}{specification of the outcome type. Must specify either "normal" or "binomial".}

\item{type}{assumption on the treatment effect: either "random" for random effects model or "fixed" for fixed effects model. Default is "random".}

\item{shrinkage}{shrinkage method applied to the effect modifiers. "none" correspond to no shrinkage.
"laplace" corresponds to a adaptive shrinkage with a Laplacian prior (ie often known as Bayesian LASSO).
"SSVS" corresponds to the Stochastic Search Variable Selection method. SSVS is not strictly a shrinkage method, but pulls the estimated coefficient toward zero through variable selection in each iteration of the MCMC. 
See O'hara et al (2009) for more details.}

\item{scale}{indicator for scaling the covariates by the overall average; default is TRUE.}

\item{mean.alpha}{prior mean for the study intercept}

\item{prec.alpha}{prior precision for the study intercept}

\item{mean.beta}{prior mean for the regression coefficients of the main effects of the covariates; main effects are assumed to have common effect.}

\item{prec.beta}{prior precision for the regression coefficients of the main effects of the covariates}

\item{mean.gamma}{prior mean for the effect modifiers. This parameter is not used if penalization is placed on effect modifiers.}

\item{prec.gamma}{prior precision for the effect modifiers. This parameter is not used if penalization is placed on effect modifiers.}

\item{mean.delta}{prior mean for the average treatment effect}

\item{prec.delta}{prior precision for the average treatment effect}

\item{hy.prior}{prior for the heterogeneity parameter. Supports uniform, gamma, and half normal for normal and binomial response
It should be a list of length 3, where first element should be the distribution (one of dunif, dgamma, dhnorm) and the next two are the parameters associated with the distribution. For example, list("dunif", 0, 5) gives uniform prior with lower bound 0 and upper bound 5 for the heterogeneity parameter.}

\item{lambda.prior}{(only for shrinkage = "laplace") two options for laplace shrinkage. We can put a gamma prior on the lambda (i.e. list("dgamma",2,0.1)) or put a uniform prior on the inverse of lambda (i.e. list("dunif",0,5))}

\item{p.ind}{(only for shrinkage = "SSVS") prior probability of including each of the effect modifiers. Length should be same as the total length of the covariates.}

\item{g}{(only for shrinkage = "SSVS") multiplier for the precision of spike. Default is g = 1000.}

\item{hy.prior.eta}{(only for shrinkage = "SSVS") standard deviation of the slab prior. Currently only support uniform distribution. Default is list("dunif", 0, 5)}
}
\value{
\item{data.JAGS}{data organized in a list so that it can be used when running code in JAGS}
\item{code}{JAGS code that is used to run the model. Use cat(code) to see the code in a readable format}
\item{model.JAGS}{JAGS code in a function. This is used when running model in parallel}
\item{scale.mean}{mean used in scaling covariates}
\item{scale.sd}{standard deviation used in scaling covariates}
}
\description{
This function sets up data and JAGS code that is needed to run one-stage IPD-MA models in JAGS.
}
\examples{
ds <- generate_ipdma_example(type = "continuous")
ipd <- with(ds, ipdma.model.onestage(y = y, study = studyid, treat = treat, X = cbind(z1, z2), 
response = "normal", shrinkage = "none"))
\donttest{
samples <- ipd.run(ipd)
treatment.effect(ipd, samples, newpatient= c(1,0.5))
} 
}
\references{
O'Hara RB, Sillanpaa MJ. A review of Bayesian variable selection methods: what, how and which. \emph{Bayesian Anal}. 2009;4(1):85-117. \doi{10.1214/09-BA403}

Seo M, White IR, Furukawa TA, et al. Comparing methods for estimating patient-specific treatment effects in individual patient data meta-analysis. \emph{Stat Med}. 2021;40(6):1553-1573. \doi{10.1002/sim.8859}
}
