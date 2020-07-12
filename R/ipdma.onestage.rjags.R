

ipdma.onestage.rjags <- function(ipd){
  
  with(ipd, {
    
    code <- paste0("\n########## IPD-MA model",
                   "\nfor (i in 1:Np) {")
    
    if(response == "binomial"){
      code <- paste0(code, "\n\ty[i] ~ dbern(p[i])",
                     "\n\tlogit(p[i]) <- a[studyid[i]] + inprod(beta[], X[i,]) +",
                     "\n\t\t(1 - equals(treat[i],1)) * inprod(gamma[], X[i,]) +")
    } else if(response == "normal"){
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
      code <- paste0(code, "\n\n#####treatment effect",
                     "\nfor(j in 1:Nstudies){",
                     "\n\td[j,1] <- 0",
                     "\n\td[j,2] ~ dnorm(delta[2], tau)",
                     "\n}")
    }
    
    if(type == "random"){
      code <- paste0(code, hy.prior.rjags(ipd))  
    }
    
    code <- paste0(code, "\n\n## prior distribution for the average treatment effect",
                   "\ndelta[1] <- 0",
                   "\ndelta[2] ~ dnorm(", mean.delta, ", ", prec.delta, ")\n")
    
    code <- paste0(code, "\n\n## prior distribution for the study intercept",
                   "\nfor (j in 1:Nstudies){",
                   "\n\ta[j] ~ dnorm(", mean.a, ", ", prec.a, ")",
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


shrinkage.prior.rjags <- function(ipd){
  
  code <- ""
  with(ipd, {
    
    if(shrinkage == "none"){
      code <- paste0(code, "\n## prior distribution for the effect modifiers under no shrinkage",
      "\nfor(k in 1:Ncovariate){",
      "\n\tgamma[k] ~ dnorm(", mean.gamma, ", ", prec.gamma, ") ",
      "\n}")
      
    } else if(shrinkage == "laplace"){
      
      code <- paste0(code, "\n## prior distribution for the effect modifiers under laplacian shrinkage")
      if(lambda.prior[[1]] == "dgamma"){
        if(response == "normal"){
          code <- paste0(code, "\ntt <- lambda * sigma")
        } else if(response == "binomial"){
          code <- paste0(code, "\ntt <- lambda")
        }  
        code <- paste0(code, "\nlambda ~ dgamma(", lambda.prior[[2]], ", ", lambda.prior[[3]], ")",
                       "\nfor(k in 1:Ncovariate){",
                       "\n\tgamma[k] ~ ddexp(0, tt)",
                       "\n}")
                       
      } else if (lambda.prior[[1]] == "dunif"){
        code <- paste0(code, "\nlambda <- pow(lambda.inv, -1)",
                       "\nlambda.inv ~ dunif(", lambda.prior[[2]], ", ", lambda.prior[[3]], ")",
                       "\nfor(k in 1:Ncovariate){",
                       "\n\tgamma[k] ~ ddexp(0, lambda)",
                       "\n}")
      }
    } else if (shrinkage == "SSVS"){
      code <- paste0(code, "\n## prior distribution for the effect modifiers under SSVS",
                    "\nfor(k in 1:Ncovariate){",
                    "\n\tIndA[k] ~ dcat(Pind[,k])",
                    "\n\tInd[k] <- IndA[k] - 1",
                    "\n\tgamma[k] ~ dnorm(0, tauCov[IndA[k]])",
                    "\n}",
                    "\n",
                    "\nzeta <- pow(eta, -2)",
                    "\neta ~ dunif(", hy.prior.eta[[2]], ", ", hy.prior.eta[[3]], ")",
                    "\ntauCov[1] <- zeta * ", g, " # precision of spike",
                    "\ntauCov[2] <- zeta # precision of slab",
                    "\n\nfor(k in 1:Ncovariate){",
                    "\n\tPind[2,k] <- p.ind[k]",
                    "\n\tPind[1,k] <- 1- p.ind[k]",
                    "\n}"
                    )
    }
    return(code)
  })
}