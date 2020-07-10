#' bnma: A package for network meta analysis using Bayesian methods
#'
#' A package for running Bayesian network meta analysis
#'
#' Network meta-analysis or mixed treatment comparison (MTC) is a method that allows simultaneous comparison of more than two treatments.
#' We use a Bayesian approach to combine both direct and indirect evidence as in Dias et al. 2013a.
#' This package is a user friendly application that can run network meta analysis models without having to code a JAGS model.
#' The program takes the input data and transforms it to a suitable format of analysis, generates a JAGS model and reasonable
#' initial values and runs the model through the rjags package.
#' The focus of this package was inclusion of multinomial response and various options for adding covariates and/or baseline risks effects.
#' Also, while sampling, the package uses Gelman-Rubin convergence criteria to decide whether to continue sampling or not.
#' Furthermore, package includes different models such as contrast based models and unrelated mean effects (UME) model and nodesplitting model to test for inconsistency.
#'
#' @docType package
#' @name bipd-package
#' @references A.J. Franchini, S. Dias, A.E. Ades, J.P. Jansen, N.J. Welton (2012), \emph{Accounting for correlation in network meta-analysis with multi-arm trials}, Research Synthesis Methods 3(2):142-160. [\url{https://doi.org/10.1002/jrsm.1049}]
#' @seealso \code{\link{ipd.data}}, \code{\link{ipd.run}}
NULL


#' @import coda
NULL


