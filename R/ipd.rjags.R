
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
    
    # code <- paste0 (code, )               
    return(code)
  })
}

