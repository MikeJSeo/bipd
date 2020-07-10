#' Run the model using the ipd object
#' 
#' This is the core function that runs the model in our program. Before running this function, we need to specify data, prior, JAGS code, etc. using \code{\link{ipd.data}}.
#' 
#' @param ipd ipd object created from \code{\link{ipd.data}} function
#' @param pars.save parameters to save. "beta" - coefficients for main effects; "gamma" - coefficients for effect modifiers; "delta" - average treaetment effect
#' @param inits initial values specified for the parameters to save
#' @param n.chains Number of MCMC chains to 
#' @param n.adapt The number of iterations for adaptation (this corresponds to the number of burnin)
#' @param n.iter The number of iterations to run after the adaptation
#'
#' @export

ipd.run <- function(ipd, pars.save = c("beta", "gamma", "delta"), inits = NULL, n.chains = 3, n.adapt = 1000, n.iter = 10000){
  
  mod <- rjags::jags.model(textConnection(ipd$code), data = ipd$data, inits = inits, n.chains = n.chains, n.adapt = n.adapt)
  samples <- rjags::coda.samples(model = mod, variable.names = pars.save, n.iter = n.iter)   

  return(samples)
}

#' Run the model using the ipd object with parallel computation
#' 
#' This function runs the model through parallel computation from dclone package. Before running this function, we need to specify data, prior, JAGS code, etc. using \code{\link{ipd.data}}.
#' 
#' @param ipd ipd object created from \code{\link{ipd.data}} function
#' @param pars.save parameters to save. "beta" - coefficients for main effects; "gamma" - coefficients for effect modifiers; "delta" - average treaetment effect
#' @param inits initial values specified for the parameters to save
#' @param n.chains Number of MCMC chains to 
#' @param n.adapt The number of iterations for adaptation (this corresponds to the number of burnin)
#' @param n.iter The number of iterations to run after the adaptation
#' @param n.cores The number of cores used for parallel computation
#'
#' @export


ipd.run.parallel <- function(ipd, pars.save = c("beta", "gamma", "delta"), inits = NULL, n.chains = 3, n.adapt = 1000, n.iter = 10000, n.cores = 3){

  cl <- parallel::makePSOCKcluster(n.cores)
  samples <- dclone::jags.parfit(cl = cl, data = ipd$data, params = pars.save, model = textConnection(ipd$code), inits = inits, n.chains = n.chains, n.adapt = n.adapt, n.iter = n.iter)
  
  return(samples)
}  