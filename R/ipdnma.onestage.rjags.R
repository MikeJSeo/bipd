

ipdnma.onestage.rjags <- function(ipd){
  
  with(ipd, {
    
    code <- paste0("\n########## IPD-NMA model",
                   "\nfor (i in 1:Np) {")
    
    if(response == "binomial"){
      
      code <- paste0(code, "\n\ty[i] ~ dbern(p[i])",
                     "\n\tlogit(p[i]) <- alpha[studyid[i]] + inprod(beta[], X[i,]) +",
                     "\n\t\t(1 - equals(treat[i],1)) * inprod(gamma[], X[i,]) +")  
      
    } else if(response == "normal"){
      
      code <- paste0(code, "\n\ty[i] ~ dnorm(mu[i], sigma)",
                     "\n\tmu[i] <- alpha[studyid[i]] + inprod(beta[], X[i,]) +",
                     "\n\t\t(1 - equals(treat[i],1)) * inprod(gamma[], X[i,]) +")  
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
    
    code <- paste0(code, shrinkage.prior.rjags(ipd))  
    code <- paste0("model {\n", code, "\n}")
    
    return(code)
  })
}
