#' bipd: A package for individual patient data meta-analysis using 'JAGS'
#'
#' A package for individual patient data meta-analysis using 'JAGS'
#'
#' We use a Bayesian approach to run individual patient data meta-analysis and network meta-analysis using 'JAGS'. 
#' The methods incorporate shrinkage methods and calculate patient-specific treatment effects as described in Seo et al. (2021) <DOI:10.1002/sim.8859>.
#' This package also includes user-friendly functions that impute missing data in an individual patient data using mice-related packages.
#'
#' @docType package
#' @name bipd-package
#' @references Audigier V, White I, Jolani S, et al. Multiple Imputation for Multilevel Data with Continuous and Binary Variables. \emph{Statistical Science}. 2018;33(2):160-183. \doi{10.1214/18-STS646}
#' @references Debray TPA, Moons KGM, Valkenhoef G, et al. Get real in individual participant data (IPD) meta-analysis: a review of the methodology. \emph{Res Synth Methods}. 2015;6(4):293-309. \doi{10.1002/jrsm.1160}
#' @references Dias S, Sutton AJ, Ades AE, et al. A Generalized Linear Modeling Framework for Pairwise and Network Meta-analysis of Randomized Controlled Trials. \emph{Medical Decision Making}. 2013;33(5):607-617. \doi{10.1177/0272989X12458724}
#' @references Fisher DJ, Carpenter JR, Morris TP, et al. Meta-analytical methods to identify who benefits most from treatments: daft, deluded, or deft approach?. \emph{BMJ}. 2017;356:j573. \doi{10.1136/bmj.j573}
#' @references O'Hara RB, Sillanpaa MJ. A review of Bayesian variable selection methods: what, how and which. \emph{Bayesian Anal}. 2009;4(1):85-117. \doi{10.1214/09-BA403}
#' @references Riley RD, Debray TP, Fisher D, et al. Individual participant data meta-analysis to examine interactions between treatment effect and participant-level covariates: Statistical recommendations for conduct and planning. \emph{Stat Med}. 2020:39(15):2115-2137. \doi{10.1002/sim.8516} 
#' @references Seo M, White IR, Furukawa TA, et al. Comparing methods for estimating patient-specific treatment effects in individual patient data meta-analysis. \emph{Stat Med}. 2021;40(6):1553-1573. \doi{10.1002/sim.8859}
NULL

#' @import coda
NULL

#' @import mvtnorm
NULL

#' @import dplyr
NULL

#' @importFrom stats update quantile rbinom rnorm runif sd
NULL
