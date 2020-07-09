
ipd.rjags <- function(ipd){
  
  
  type <- ipd$type
  model <- ipd$model
  
  code <- if(model == "onestage"){
    model.onestage(ipd)  
  } 
  code <- paste0("model {\n", code, "\n}")
  
  return(code)
}


model.onestage <- function(ipd){
  
  with(ipd, {
    
    code <- paste0("\n########## IPD-MA model",
                   "\nfor (i in 1:Np) {")
    
    if(response == "binary"){
      code <- paste0(code, "\n\ty[i] ~ dbern(p[i])",
                     "\n\tlogit(p[i]) <- a[studyid[i]] + inprod(beta[], X[i,]) +",
                     "\n\t\t(1 - equals(treat[i],1)) * inprod(gamma[], X[i,])")
    } else if(response == "continuous"){
      code <- paste0(code, "\n\ty[i] ~ dnorm(mu[i], sigma)",
                     "\n\tmu[i] <- a[studyid[i]] + inprod(beta[], X[i,]) +",
                     "\n\t\t(1 - equals(treat[i],1)) * inprod(gamma[], X[i,]) +")
    }
    
    if(type == "random"){
      code <- paste0(code, " d[studyid[i],treat[i]]",
                           "\n}")
    } else if (type == "fixed"){
      code <- paste0(code, "delta[treat[i]]",
                           "\n}")
    }
    
    if(response == "normal"){
      code <- paste0(code, "\nsigma ~ dgamma(0.001, 0.001)")
    }
    
    if(type == "random"){
      code <- paste0(code, "\n#####treatment effect",
                     "\nfor(j in 1:Nstudies){",
                     "\n\td[j,1] <- 0",
                     "\n\td[j,2] ~ dnorm(delta[2], tau)",
                     "\n}")
    }
    
    if(type == "random"){
      code <- paste0(code, hy.prior.rjags(hy.prior))  
    }
    
    code <- paste0(code, "\n## prior distribution for baseline effect",
                   "\nfor (j in 1:Nstudies){",
                   "\n\talpha[j] ~ dnorm(", mean.alpha, ", ", prec.alpha, ")",
                   "\n}")
    
    code <- paste0(code, "\nfor(k in 1:Ncovariate){",
                   "\nfor(k in 1:Ncovariate){",
                   "\n\tbeta[k] ~ dnorm(", mean.beta, ", ", prec.beta, ")",
                   "\n}")

    code <- paste0(code, "## prior distribution for average treatment effect",
                   "delta[1] <- 0",
                   "delta[2] ~ dnorm(", mean.delta, ", ", prec.delta, ")")
                   

    

    return(code)
  })
}



hy.prior.rjags <- function(hy.prior){
  
  code <- ""
  distr <- hy.prior[[1]]
  if (distr == "dunif") {
    code <- paste0(code,
                   "\n\tsd ~ dunif(hy.prior.1, hy.prior.2)",
                   "\n\ttau <- pow(sd,-2)")
  } else if(distr == "dgamma"){
    code <- paste0(code,
                   "\n\tsd <- pow(prec, -0.5)",
                   "\n\ttau ~ dgamma(hy.prior.1, hy.prior.2)")
  } else if(distr == "dhnorm"){
    code <- paste0(code,
                   "\n\tsd ~ dnorm(hy.prior.1, hy.prior.2)T(0,)",
                   "\n\ttau <- pow(sd, -2)")
  }
}
