% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ipdma.onestage.deft.model.R
\name{ipdma.model.deft.onestage}
\alias{ipdma.model.deft.onestage}
\title{Make a (Deft-approach) one-stage individual patient data meta-analysis (IPD-MA) object containing data, priors, and a JAGS model code}
\usage{
ipdma.model.deft.onestage(
  y = NULL,
  study = NULL,
  treat = NULL,
  X = NULL,
  response = "normal",
  type = "random",
  mean.alpha = 0,
  prec.alpha = 0.001,
  mean.beta = 0,
  prec.beta = 0.001,
  mean.gamma.within = 0,
  prec.gamma.within = 0.001,
  mean.gamma.across = 0,
  prec.gamma.across = 0.001,
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
\item{y}{Outcome of the study. Can be continuous or binary.}

\item{study}{A vector indicating which study the patient belongs to. Please change the study names into numbers (i.e. 1, 2, 3, etc)}

\item{treat}{A vector indicating which treatment the patient was assigned to (i.e. 1 for treatment, 0 for placebo)}

\item{X}{A matrix of covariate values for each patient. Dimension would be number of patients x number of covariates.}

\item{response}{Specification of the outcome type. Must specify either "normal" or "binomial".}

\item{type}{Assumption on the treatment effect: either "random" for random effects model or "fixed" for fixed effects model. Default is "random".}

\item{mean.alpha}{Prior mean for the study intercept}

\item{prec.alpha}{Prior precision for the study intercept}

\item{mean.beta}{Prior mean for the regression coefficients of the main effects of the covariates; main effects are assumed to have common effect.}

\item{prec.beta}{Prior precision for the regression coefficients of the main effects of the covariates}

\item{mean.gamma.within}{Prior mean for effect modifiers of within study information.}

\item{prec.gamma.within}{Prior precision for the effect modifiers of within study information.}

\item{mean.gamma.across}{Prior mean for the effect modifiers of across study information; effect modification is assumed to have common effect.}

\item{prec.gamma.across}{Prior precision for the effect modifiers of across study information}

\item{mean.delta}{Prior mean for the average treatment effect}

\item{prec.delta}{Prior precision for the average treatment effect}

\item{hy.prior}{Prior for the heterogeneity parameter. Supports uniform, gamma, and half normal for normal and binomial response
It should be a list of length 3, where first element should be the distribution (one of dunif, dgamma, dhnorm) and the next two are the parameters associated with the distribution. For example, list("dunif", 0, 5) give uniform prior with lower bound 0 and upper bound 5 for the heterogeneity parameter.}
}
\value{
\item{data.JAGS}{Data organized in a list so that it can be used when running code in JAGS}
\item{code}{JAGS code that is used to run the model. Use cat(code) to see the code in a readable format}
\item{model.JAGS}{JAGS code in a function. This is used when running model in parallel}
\item{Xbar}{Study specific average of covariates}
}
\description{
This function sets up data and JAGS code that is needed to run one-stage IPD-MA models in JAGS.
}
\references{
Fisher DJ, Carpenter JR, Morris TP, et al. Meta-analytical methods to identify who benefits most from treatments: daft, deluded, or deft approach?. \emph{BMJ}. 2017;356:j573 \doi{10.1136/bmj.j573}
}