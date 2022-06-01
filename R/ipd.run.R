#' Run the model using the ipd object
#' 
#' This is the core function that runs the model in our program. Before running this function, we need to specify data, prior, JAGS code, etc. using ipd.model type function.
#' 
#' @param ipd ipd object created from ipd.model type function
#' @param pars.save parameters to save. For instance, "beta" - coefficients for main effects; "gamma" - coefficients for effect modifiers; "delta" - average treatment effect
#' @param inits initial values specified for the parameters to save
#' @param n.chains number of MCMC chains to sample
#' @param n.adapt number of iterations for adaptation (Note that the samples from adaptation phase is non-Markovian and do not constitute a Markov chain)
#' @param n.burnin number of iterations for burn-in
#' @param n.iter number of iterations to run after the adaptation
#' @return MCMC samples stored using JAGS. The returned samples have the form of mcmc.list and coda functions can be directly applied.
#' @examples
#' ds <- generate_ipdma_example(type = "continuous")
#' ipd <- with(ds, ipdma.model.onestage(y = y, study = studyid, treat = treat, X = cbind(z1, z2), 
#' response = "normal", shrinkage = "none"))
#' \donttest{
#' samples <- ipd.run(ipd, n.chains = 3, n.burnin = 500, n.iter = 5000)
#' }
#' @export

ipd.run <- function(ipd, pars.save = NULL, inits = NULL, n.chains = 3, n.adapt = 1000, n.burnin = 1000, n.iter = 10000){
  
  if(is.null(pars.save)){
    # default save parameters
    if(class(ipd) %in% c("ipdma.onestage", "ipdnma.onestage", "ipdnma.twostage.second")){
      pars.save <- c("alpha", "beta", "gamma", "delta")
    } else if(class(ipd) %in% c("ipdma.onestage.deft")){
      pars.save <- c("alpha", "beta", "gamma.within", "gamma.across", "delta")
    } else if(class(ipd) %in% c("ipdnma.twostage.first")){
      pars.save <- c("a", "b", "c", "d")
    }
  }
 
  mod <- rjags::jags.model(textConnection(ipd$code), data = ipd$data.JAGS, inits = inits, n.chains = n.chains, n.adapt = n.adapt)
  stats::update(mod, n.burnin)
  samples <- rjags::coda.samples(model = mod, variable.names = pars.save, n.iter = n.iter)   

  return(samples)
}

#' Run the model using the ipd object with parallel computation
#' 
#' This function runs the model through parallel computation using dclone R package. Before running this function, we need to specify data, prior, JAGS code, etc. using ipd.model type function.
#' 
#' @param ipd ipd object created from ipd.model type function
#' @param pars.save parameters to save. For instance, "beta" - coefficients for main effects; "gamma" - coefficients for effect modifiers; "delta" - average treatment effect
#' @param inits initial values specified for the parameters to save
#' @param n.chains number of MCMC chains to sample
#' @param n.adapt number of iterations for adaptation (Note that the samples from adaptation phase is non-Markovian and do not constitute a Markov chain)
#' @param n.burnin number of iterations for burn-in
#' @param n.iter number of iterations to run after the adaptation
#' @return MCMC samples stored using JAGS. The returned samples have the form of mcmc.list and coda functions can be directly applied.
#' @examples
#' ds <- generate_ipdma_example(type = "continuous")
#' ipd <- with(ds, ipdma.model.onestage(y = y, study = studyid, treat = treat, X = cbind(z1, z2), 
#' response = "normal", shrinkage = "none"))
#' \donttest{
#' samples <- ipd.run.parallel(ipd, n.chains = 2, n.burnin = 500, n.iter = 5000)
#' }
#' @export


ipd.run.parallel <- function(ipd, pars.save = NULL, inits = NULL, n.chains = 2, n.adapt = 1000, n.burnin = 1000, n.iter = 10000){

  if(is.null(pars.save)){
    # default save parameters
    if(class(ipd) %in% c("ipdma.onestage", "ipdnma.onestage", "ipdnma.twostage.second")){
      pars.save <- c("alpha", "beta", "gamma", "delta")
    } else if(class(ipd) %in% c("ipdma.onestage.deft")){
      pars.save <- c("alpha", "beta", "gamma.within", "gamma.across", "delta")
    } else if(class(ipd) %in% c("ipdnma.twostage.first")){
      pars.save <- c("a", "b", "c", "d")
    }
  }
  
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  
  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    num_workers <- 2L
  } else {
    num_workers <- n.chains
  }
  
  cl <- parallel::makePSOCKcluster(num_workers)
  samples <- dclone::jags.parfit(cl = cl, data = ipd$data.JAGS, params = pars.save, model = ipd$model.JAGS, inits = inits, n.chains = n.chains, n.adapt = n.adapt, n.update = n.burnin, n.iter = n.iter)
  parallel::stopCluster(cl)
  
  return(samples)
}  