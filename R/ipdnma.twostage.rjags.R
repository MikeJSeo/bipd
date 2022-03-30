
ipdma.twostage.first.rjags <- function(ipd){
  
  with(ipd, {
    
    code <- paste0("\n########## ipdnma two-stage model - first stage",
                   "\nfor (i in 1:Np) {")
    
    if(response == "binomial"){
      
      code <- paste0(code, "\n\ty[i] ~ dbern(p[i])",
                     "\n\tlogit(p[i]) <- a + inprod(b[], X[i,]) +",
                     "\n\t\tinprod(c[treat[i],], X[i,]) + d[treat[i]]") 
      
    } else if(response == "normal"){
      
      code <- paste0(code, "\n\ty[i] ~ dnorm(mu[i], sigma)",
                     "\n\tmu[i] <- a + inprod(b[], X[i,]) +",
                     "\n\t\tinprod(c[treat[i],], X[i,]) + d[treat[i]]")  
    }
    
    if(response == "normal"){
      code <- paste0(code, "\nsigma ~ dgamma(0.001, 0.001)")
    }
    
    code <- paste0(code, "\n\n## prior distribution for the treatment effect",
                   "\nd[1] <- 0",
                   "\nfor(k in 2:Ntreat){",
                   "\n\td[k] ~ dnorm(", mean.d, ", ", prec.d, ")",
                   "\n}")
    
    code <- paste0(code, "\n\n## prior distribution for the study intercept",
                   "\na ~ dnrom(", mean.a, ", ", prec.a, ")",
                   "\n}")
                   
    code <- paste0(code, "\n\n## prior distribution for the main effect of the covariates",
                   "\nfor(k in 1:Ncovariate){",
                   "\n\tb[k] ~ dnorm(", mean.b, ", ", prec.b, ")",
                   "\n}")
    
    code <- paste0(code, nma.shrinkage.prior.rjags(ipd))  
    code <- paste0("model {\n", code, "\n}")
    
    return(code)
  })
}