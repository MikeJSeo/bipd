#' bipd: A package for individual patient data (network) meta analysis using 'JAGS'
#'
#' A package for running Bayesian individual patient data (network) meta analysis
#'
#' We use a Bayesian approach to run individual patient data (network) meta analysis. The methods allow shrinkage of effect modifiers as described in Seo et al.
#' This package is a user friendly application that can run models without having to code a JAGS model.
#' The program takes the input data and transforms it to a suitable format of analysis, generates a JAGS model and runs the model through the rjags package.
#' Furthermore, package includes functions to calculate measures such as patient-specific treatment effect.
#'
#' @docType package
#' @name bipd-package
#' @references Debray TPA, Moons KGM, Valkenhoef G, et al. Get real in individual participant data (IPD) meta-analysis: a review of the methodology. \emph{Res Synth Methods}. 2015;6(4):293-309. \doi{10.1002/jrsm.1160}
#' @references Fisher DJ, Carpenter JR, Morris TP, et al. Meta-analytical methods to identify who benefits most from treatments: daft, deluded, or deft approach?. \emph{BMJ}. 2017;356:j573 \doi{10.1136/bmj.j573}
#' @references Riley RD, Debray TP, Fisher D, et al. Individual participant data meta-analysis to examine interactions between treatment effect and participant-level covariates: Statistical recommendations for conduct and planning. \emph{Stat Med}. 2020:39(15):2115-2137. \doi{10.1002/sim.8516} 
#' @references Seo M, White IR, Furukawa TA, et al. Comparing methods for estimating patient-specific treatment effects in individual patient data meta-analysis. \emph{Stat Med}. 2021;40(6):1553-1573. \doi{10.1002/sim.8859}
NULL


#' @importFrom stats "update"
NULL

#' @import coda
NULL

#' @importFrom dclone "jags.parfit"
NULL

#' @import R2WinBUGS
NULL

#' @import mvtnorm
NULL

#' @import tidyr
NULL

#' @import dplyr 
NULL
