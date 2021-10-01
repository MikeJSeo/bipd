#' Run the model using the ipd object
#' 
#' This is the core function that runs the model in our program. Before running this function, we need to specify data, prior, JAGS code, etc. using ipd.data function.
#' 
#' @param ipd ipd object created from ipd.data function
#' @param pars.save parameters to save. "beta" - coefficients for main effects; "gamma" - coefficients for effect modifiers; "delta" - average treaetment effect
#' @param inits initial values specified for the parameters to save
#' @param n.chains Number of MCMC chains to sample
#' @param n.adapt The number of iterations for adaptation (Note that the samples from adaptation phase is non-Markovian and do not constitute a Markov chain)
#' @param n.burnin The number of iterations for burn-in
#' @param n.iter The number of iterations to run after the adaptation
#'
#' @export

ipd.run <- function(ipd, pars.save = c("beta", "gamma", "delta"), inits = NULL, n.chains = 3, n.adapt = 1000, n.burnin = 1000, n.iter = 10000){
  
  mod <- rjags::jags.model(textConnection(ipd$code), data = ipd$data.JAGS, inits = inits, n.chains = n.chains, n.adapt = n.adapt)
  rjags::update(mod, n.burnin)
  samples <- rjags::coda.samples(model = mod, variable.names = pars.save, n.iter = n.iter)   

  return(samples)
}

#' Run the model using the ipd object with parallel computation
#' 
#' This function runs the model through parallel computation from dclone package. Before running this function, we need to specify data, prior, JAGS code, etc. using ipd.data function.
#' 
#' @param ipd ipd object created from ipd.data function
#' @param pars.save parameters to save. "beta" - coefficients for main effects; "gamma" - coefficients for effect modifiers; "delta" - average treatment effect
#' @param inits initial values specified for the parameters to save
#' @param n.chains Number of MCMC chains to run. This corresponds the number of cores used for parallel computation.
#' @param n.adapt The number of iterations for adaptation (Note that the samples from adaptation phase is non-Markovian and do not constitute a Markov chain)
#' @param n.burnin The number of iterations for burn-in
#' @param n.iter The number of iterations to run after the adaptation
#'
#' @export


ipd.run.parallel <- function(ipd, pars.save = c("beta", "gamma", "delta"), inits = NULL, n.chains = 3, n.adapt = 1000, n.burnin = 1000, n.iter = 10000){

  cl <- parallel::makePSOCKcluster(n.chains)
  samples <- dclone::jags.parfit(cl = cl, data = ipd$data.JAGS, params = pars.save, model = ipd$model.JAGS, inits = inits, n.chains = n.chains, n.adapt = n.adapt, n.update = n.burnin, n.iter = n.iter)
  parallel::stopCluster(cl)
  
  return(samples)
}  