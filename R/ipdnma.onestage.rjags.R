

ipdnma.onestage.rjags <- function(ipd){
  
  with(ipd, {
    
    code <- paste0("\n########## IPD-NMA model",
                   "\nfor (i in 1:Np) {")
    
    if(response == "binomial"){
      
      code <- paste0(code, "\n\ty[i] ~ dbern(p[i])",
                     "\n\tlogit(p[i]) <- alpha[studyid[i]] + inprod(beta[], X[i,]) +",
                     "\n\t\tinprod(gamma[treat[i],], X[i,]) +")  
      
    } else if(response == "normal"){
      
      code <- paste0(code, "\n\ty[i] ~ dnorm(mu[i], sigma)",
                     "\n\tmu[i] <- alpha[studyid[i]] + inprod(beta[], X[i,]) +",
                     "\n\t\tinprod(gamma[treat[i],], X[i,]) +")  
    }
    
    if(type == "random"){
      code <- paste0(code, " d[studyid[i],treatment.arm[i]]",
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
                     "\nfor(i in 1:Nstudies){",
                     "\n\tw[i,1] <- 0",
                     "\n\td[i,1] <- 0",
                     "\n\tfor(k in 2:na[i]){",
                     "\n\t\td[i,k] ~ dnorm(mdelta[i,k], taudelta[i,k])",
                     "\n\t\tmdelta[i,k] <-  delta[t[i,k]] - delta[t[i,1]] + sw[i,k]",
                     "\n\t\ttaudelta[i,k] <- tau * 2 * (k-1)/k",
                     "\n\t\tw[i,k] <- d[i,k] - delta[t[i,k]] + delta[t[i,1]]",
                     "\n\t\tsw[i,k] <- sum(w[i, 1:(k-1)]) / (k-1)",
                     "\n\t}",
                     "\n}")
      code <- paste0(code, hy.prior.rjags(ipd))
    }

    code <- paste0(code, "\n\n## prior distribution for the average treatment effect",
                   "\ndelta[1] <- 0",
                   "\nfor(k in 2:Ntreat){",
                   "\n\tdelta[k] ~ dnorm(", mean.delta, ", ", prec.delta, ")",
                   "\n}")

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
                     "\n\tgamma[1,k] <- 0",
                     "\n\tfor(m in 2:Ntreat){",
                     "\n\t\tgamma[m,k] ~ dnorm(", mean.gamma, ", ", prec.gamma, ") ",
                     "\n\t}",
                     "\n}")
      
    } else if(shrinkage == "laplace"){
      
      code <- paste0(code, "\n## prior distribution for the effect modifiers under laplacian shrinkage")
      if(lambda.prior[[1]] == "dgamma"){
        if(response == "normal"){
          code <- paste0(code, 
                         "\nlambda[1] <- 0",
                         "\nfor(m in 2:Ntreat){",
                         "\n\ttt[m] <- lambda[m] * sigma",
                         "\n\tlambda[m] ~ dgamma(", lambda.prior[[2]], ", ", lambda.prior[[3]], ")",
                         "\n}"
                         )
        } else if(response == "binomial"){
          code <- paste0(code,
                         "\nlambda[1] <- 0",
                         "\nfor(m in 2:Ntreat){",
                         "\n\ttt[m] <- lambda[m]",
                         "\n\tlambda[m] ~ dgamma(", lambda.prior[[2]], ", ", lambda.prior[[3]], ")",
                         "\n}"
                         )
        } 
        code <- paste0(code, 
                       "\nfor(k in 1:Ncovariate){",
                       "\n\tgamma[1,k] <- 0",
                       "\n\tfor(m in 2:Ntreat){",
                       "\n\t\tgamma[m,k] ~ ddexp(0, tt[m])",
                       "\n\t}",
                       "\n}")
      } else if (lambda.prior[[1]] == "dunif"){
        if(response == "normal"){
          code <- paste0(code, 
                         "\nlambda[1] <- 0",
                         "\nlambda.inv[1] <- 0",
                         "\nfor(m in 2:Ntreat){",
                         "\n\ttt[m] <- lambda[m] * sigma",
                         "\n\tlambda[m] <- pow(lambda.inv[m], -1)",
                         "\n\tlambda.inv[m] ~ dunif(", lambda.prior[[2]], ", ", lambda.prior[[3]], ")",
                         "\n}")
        } else if(response == "binomial"){
          code <- paste0(code, 
                         "\nlambda[1] <- 0",
                         "\nlambda.inv[1] <- 0",
                         "\nfor(m in 2:Ntreat){",
                         "\n\ttt[m] <- lambda[m]",
                         "\n\tlambda[m] <- pow(lambda.inv[m], -1)",
                         "\n\tlambda.inv[m] ~ dunif(", lambda.prior[[2]], ", ", lambda.prior[[3]], ")",
                         "\n}")
        }
        code <- paste0(code, 
                       "\nfor(k in 1:Ncovariate){",
                       "\n\tgamma[1,k] <- 0",
                       "\n\tfor(m in 2:Ntreat){",
                       "\n\t\tgamma[m,k] ~ ddexp(0, tt[m])",
                       "\n\t}",
                       "\n}")
      }
    } else if (shrinkage == "SSVS"){
      code <- paste0(code, "\n## prior distribution for the effect modifiers under SSVS",
                     "\nfor(k in 1:Ncovariate){",
                     "\n\tgamma[1,k] <- 0",
                     "\n\tfor(m in 2:Ntreat){",
                     "\n\t\tIndA[m,k] ~ dcat(Pind[,k])",
                     "\n\t\tInd[m,k] <- IndA[m,k] - 1",
                     "\n\t\tgamma[m,k] ~ dnorm(0, tauCov[IndA[m,k]])",
                     "\n\t}",
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

