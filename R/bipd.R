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
#' @references A.J. Franchini, S. Dias, A.E. Ades, J.P. Jansen, N.J. Welton (2012), \emph{Accounting for correlation in network meta-analysis with multi-arm trials}, Research Synthesis Methods 3(2):142-160. [\url{https://doi.org/10.1002/jrsm.1049}]
NULL


#' @import stats
NULL

#' @import coda
NULL

#' @import dclone
NULL

#' @import R2WinBUGS
NULL
