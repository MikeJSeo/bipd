

ipdnma.onestage.rjags <- function(ipd){
  
  with(ipd, {
    
    code <- paste0("\n########## IPD-NMA model",
                   "\nfor (i in 1:Np) {")
    
    if(response == "binomial"){
      
      code <- paste0(code, "\n\ty[i] ~ dbern(p[i])",
                     "\n\tlogit(p[i]) <- alpha[studyid[i]] + inprod(beta[], X[i,]) +",
                     "\n\t\tinprod(gamma[,treat[i]], X[i,]) ")  
      
    } else if(response == "normal"){
      
      code <- paste0(code, "\n\ty[i] ~ dnorm(mu[i], sigma)",
                     "\n\tmu[i] <- alpha[studyid[i]] + inprod(beta[], X[i,]) +",
                     "\n\t\tinprod(gamma[,treat[i]], X[i,]) ")  
    }
    
    code <- paste0(code, "\n}")
    
    # if(type == "random"){
    #   code <- paste0(code, " d[studyid[i],treat[i]]",
    #                  "\n}")
    # } else if (type == "fixed"){
    #   code <- paste0(code, " delta[treat[i]]",
    #                  "\n}")
    # }
    
    if(response == "normal"){
      code <- paste0(code, "\nsigma ~ dgamma(0.001, 0.001)")
    }
    
    # if(type == "random"){
    #   code <- paste0(code, "\n\n#####treatment effect",
    #                  "\nfor(j in 1:Nstudies){",
    #                  "\n\td[j,1] <- 0",
    #                  "\n\td[j,2] ~ dnorm(delta[2], tau)",
    #                  "\n}")
    #   code <- paste0(code, hy.prior.rjags(ipd))  
    # }
    # 
    # code <- paste0(code, "\n\n## prior distribution for the average treatment effect",
    #                "\ndelta[1] <- 0",
    #                "\ndelta[2] ~ dnorm(", mean.delta, ", ", prec.delta, ")\n")
    
    code <- paste0(code, "\n\n## prior distribution for the study intercept",
                   "\nfor (j in 1:Nstudies){",
                   "\n\talpha[j] ~ dnorm(", mean.alpha, ", ", prec.alpha, ")",
                   "\n}")
    
    code <- paste0(code, "\n\n## prior distribution for the main effect of the covariates",
                   "\nfor(k in 1:Ncovariate){",
                   "\n\tbeta[k] ~ dnorm(", mean.beta, ", ", prec.beta, ")",
                   "\n}")
    
    code <- paste0(code, nma.shrinkage.prior.rjags(ipd))  
    code <- paste0("model {\n", code, "\n}")
    
    return(code)
  })
}




nma.shrinkage.prior.rjags <- function(ipd){
  
  code <- ""
  with(ipd, {
    
    if(shrinkage == "none"){
      code <- paste0(code, "\n## prior distribution for the effect modifiers under no shrinkage",
                     "\nfor(k in 1:Ncovariate){",
                     "\n\tgamma[k,1] <- 0",
                     "\n\tfor(l in 2:Ntreat){",
                     "\n\t\tgamma[k,l] ~ dnorm(", mean.gamma, ", ", prec.gamma, ") ",
                     "\n\t}",
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
        if(response == "normal"){
          code <- paste0(code, "\ntt <- lambda * sigma")
        } else if(response == "binomial"){
          code <- paste0(code, "\ntt <- lambda")
        }
        code <- paste0(code, "\nlambda <- pow(lambda.inv, -1)",
                       "\nlambda.inv ~ dunif(", lambda.prior[[2]], ", ", lambda.prior[[3]], ")",
                       "\nfor(k in 1:Ncovariate){",
                       "\n\tgamma[k] ~ ddexp(0, tt)",
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

