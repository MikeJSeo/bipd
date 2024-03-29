ipdma.onestage.deft.rjags <- function(ipd){
  
  with(ipd, {
    
    code <- paste0("\n########## IPD-MA model",
                   "\nfor (i in 1:Np) {")
    
    if(response == "binomial"){
      
      code <- paste0(code, "\n\ty[i] ~ dbern(p[i])",
                       "\n\tlogit(p[i]) <- alpha[studyid[i]] + inprod(beta[], X[i,]) + (1 - equals(treat[i],1)) * inprod(gamma.across[], Xbar[studyid[i],]) +",
                       "\n\t\t(1 - equals(treat[i],1)) * inprod(gamma.within[], X[i,] - Xbar[studyid[i],]) +")
    } else if(response == "normal"){
      
      code <- paste0(code, "\n\ty[i] ~ dnorm(mu[i], sigma)",
                       "\n\tmu[i] <- alpha[studyid[i]] + inprod(beta[], X[i,]) + (1 - equals(treat[i],1)) * inprod(gamma.across[], Xbar[studyid[i],]) +",
                       "\n\t\t(1 - equals(treat[i],1)) * inprod(gamma.within[], X[i,] - Xbar[studyid[i],]) +") 
    }
    
    if(type == "random"){
      code <- paste0(code, " d[studyid[i],treat[i]]",
                     "\n}")
    } else if (type == "fixed"){
      code <- paste0(code, " delta[treat[i]]",
                     "\n}")
    }
    
    if(response == "normal"){
      code <- paste0(code, "\nsigma ~ dgamma(0.001, 0.001)")
    }
    
    if(type == "random"){
      code <- paste0(code, "\n\n#####treatment effect",
                     "\nfor(j in 1:Nstudies){",
                     "\n\td[j,1] <- 0",
                     "\n\td[j,2] ~ dnorm(delta[2], tau)",
                     "\n}")
      code <- paste0(code, hy.prior.rjags(ipd))  
    }
    
    code <- paste0(code, "\n\n## prior distribution for the average treatment effect",
                   "\ndelta[1] <- 0",
                   "\ndelta[2] ~ dnorm(", mean.delta, ", ", prec.delta, ")\n")
    
    code <- paste0(code, "\n\n## prior distribution for the study intercept",
                   "\nfor (j in 1:Nstudies){",
                   "\n\talpha[j] ~ dnorm(", mean.alpha, ", ", prec.alpha, ")",
                   "\n}")
    
    code <- paste0(code, "\n\n## prior distribution for the main effect of the covariates",
                   "\nfor(k in 1:Ncovariate){",
                   "\n\tbeta[k] ~ dnorm(", mean.beta, ", ", prec.beta, ")",
                   "\n}")
    
    code <- paste0(code, "\n\n## prior distribution for the effect modifiers of across study information",
                     "\nfor(k in 1:Ncovariate){",
                     "\n\tgamma.across[k] ~ dnorm(", mean.gamma.across, ", ", prec.gamma.across, ")",
                     "\n}")

    
    code <- paste0(code, "\n## prior distribution for the effect modifiers of within study information",
                   "\nfor(k in 1:Ncovariate){",
                   "\n\tgamma.within[k] ~ dnorm(", mean.gamma.within, ", ", prec.gamma.within, ") ",
                   "\n}")
    
    code <- paste0("model {\n", code, "\n}")
    
    return(code)
  })
}



hy.prior.rjags <- function(ipd){
  
  code <- ""
  hy.prior <- ipd$hy.prior
  distr <- hy.prior[[1]]
  if (distr == "dunif") {
    code <- paste0(code,
                   "\nsd ~ dunif(", hy.prior[[2]], ", ", hy.prior[[3]], ")",
                   "\ntau <- pow(sd,-2)")
  } else if(distr == "dgamma"){
    code <- paste0(code,
                   "\nsd <- pow(prec, -0.5)",
                   "\ntau ~ dgamma(", hy.prior[[2]], ", ", hy.prior[[3]], ")")
  } else if(distr == "dhnorm"){
    code <- paste0(code,
                   "\nsd ~ dnorm(", hy.prior[[2]], ", ", hy.prior[[3]], ")T(0,)",
                   "\ntau <- pow(sd, -2)")
  }
}

