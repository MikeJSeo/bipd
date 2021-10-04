
#' Impute missing data in individual participant data with two treatments (i.e. placebo and a treatment).
#'
#' Impute missing data in individual participant data with two treatments. Data is clustered by different studies.
#' In the presence of systematically missing variables, the function uses 2l.2stage.norm, 2l.2stage.bin, 
#' and 2l.2stage.pois methods in micemd package depending on the type of the covariates to impute. 
#' If there are no systematically missing variables, the function defaults to use 2l.pmm in miceadds package 
#' which generalizes predictive mean matching using linear mixed model.
#' If there is only one study available, the function defaults to use pmm in mice package.
#'
#' @param dataset Data which contains variables of interests
#' @param covariates Vector of variable names to find missing data pattern
#' @param typeofvar Type of covariate variables; should be a vector of these values: "continuous", "binary", or "count".
#' Order should follow that of covariates parameter specified. Covariates that are specified "binary" are automatically factored.
#' @param interaction Indicator denoting whether treatment-covariate interactions should be included. Default is set to true.
#' @param meth User can specify imputation method to be used in the mice package. If left unspecified, function picks a reasonable one.
#' @param pred User can specify correct prediction matrix to be used in the mice package. If left unspecified, function picks a reasonable one.
#' @param studyname Study name in the data specified.
#' @param treatmentname Treatment name in the data specified.
#' @param outcomename Outcome name in the data specified.
#' @param m Number of imputed datasets. Default is set to 5.
#' @return 
#' \item{missingPattern}{missing Pattern object returned by running \code{\link{findMissingPattern}}}
#' \item{meth}{imputation method used with the mice function}
#' \item{pred}{prediction matrix used with the mice function}
#' \item{imp}{imputed datasets that is returned from the mice function}
#' \item{imp.list}{imputed datasets in a list format}
#'
#' @export

ipdma.impute <- function(dataset = NULL, covariates = NULL, typeofvar = NULL, interaction = TRUE,
                         meth = NULL, pred = NULL, studyname = NULL, treatmentname = NULL, outcomename = NULL, 
                         m = 5
                         ){
  
  if(is.null(dataset) | is.null(covariates) | is.null(typeofvar)){
    stop("dataset, covariates, and typeofvar have to be specified.")
  }
  
  if(is.null(studyname) | is.null(treatmentname) | is.null(outcomename)){
    stop("studyname, treatmentname, and outcomename have to be specified.")
  }
  
  if(length(typeofvar) != length(covariates)){
    stop("length of covariates and typeofvar should match")
  } else{
    names(typeofvar) <- covariates  
  }
  
  if(treatmentname %in% colnames(dataset)){
    dataset <- dataset %>% select(all_of(c(studyname, treatmentname, outcomename, covariates)))  
  } else{
    dataset <- dataset %>% select(all_of(c(studyname, outcomename, covariates)))
  }
  
  dataset[,covariates] <- dataset[,covariates] %>% mutate_if(typeofvar == "binary", as.factor)
  
  if(interaction == TRUE){
    for(i in 1:length(covariates)){
      varname <- paste0(covariates[i], treatmentname)
      dataset[[varname]] <- NA
    }
  }

  missingPattern <- findMissingPattern(dataset = dataset, covariates = covariates, typeofvar = typeofvar, 
                                       studyname = studyname, treatmentname = treatmentname, outcomename = outcomename)
  
  if(is.null(meth)){
    meth <- getCorrectMeth(dataset = dataset, missingPattern = missingPattern, interaction = interaction)
  }
  
  if(is.null(pred)){
    pred <- getCorrectPred(dataset = dataset, missingPattern = missingPattern, interaction = interaction)
  }

  imp <- mice(dataset, pred = pred, meth = meth, m = m)
  
  impc <- complete(imp, "long", include = "TRUE")
  impc.store <- impc[, c(".imp", ".id", studyname, outcomename, covariates, grep(treatmentname, colnames(impc), value = TRUE))]
  imp.list <- mitools::imputationList(split(impc.store, impc.store[,1])[-1])$imputations

  list(missingPattern = missingPattern, meth = meth, pred = pred, imp = imp, imp.list = imp.list)
}
  

#' Find missing data pattern in a given data
#'
#' Find missing data pattern in a given data i.e. whether they are systematically missing or sporadically missing. Also calculates missing count and percentage for exploratory purposes.
#'
#' @param dataset Data which contains variables of interests
#' @param covariates Vector of variable names that the user is interested in finding a missing data pattern
#' @param typeofvar Type of covariate variables; should be a vector of these values: "continuous", "binary", or "count". Order should follow that of covariates parameter.
#' @param studyname Study name in the data specified.
#' @param treatmentname Treatment name in the data specified.
#' @param outcomename Outcome name in the data specified.
#' 
#' @export

findMissingPattern <- function(dataset = NULL, covariates = NULL, typeofvar = NULL, 
                               studyname = NULL, treatmentname = NULL, outcomename = NULL){
  
  if(is.null(dataset) | is.null(covariates) | is.null(typeofvar)){
    stop("All parameters dataset, covariates, and typeofvar is not given")
  }
  
  if(is.null(studyname) | is.null(treatmentname) | is.null(outcomename)){
    stop("studyname, treatmentname, and outcomename have to be specified.")
  }
  
  if(treatmentname %in% covariates){
    stop("Treatment name should not be included as covariates")
  }
  
  if(outcomename %in% covariates){
    stop("Outcome name should not be included as covariates")
  }
  
  totalstudies <- dataset %>% select(studyname) %>% n_distinct()

  if(totalstudies == 1){
    
    missingcount <- dataset %>% select(all_of(covariates))%>% summarize_all(~sum(is.na(.)))
    missingn <- dataset %>% select(all_of(covariates))  %>% summarise_all(~length(.))
    missingpercent <- dataset %>% select(all_of(covariates)) %>% summarise_all(~round(sum(is.na(.)/length(.))*100))
    
    numberofNAs <- dataset %>% select(all_of(covariates)) %>% summarize_all(~sum(is.na(.)))
    studysize = dim(dataset)[1]
    
    sys_missing <- numberofNAs == studysize
    spor_missing <- numberofNAs < studysize & numberofNAs > 0
    sys_covariates <- colnames(sys_missing)[which(sys_missing == TRUE)]
    spor_covariates <- colnames(spor_missing)[which(spor_missing == TRUE)]
    without_sys_covariates <- colnames(sys_missing)[which(sys_missing == FALSE)]
    
  } else{
    
    missingcount <- dataset %>% select(all_of(c(studyname, covariates))) %>% group_by(across(all_of(studyname))) %>% summarize_all(~sum(is.na(.)))
    missingn <- dataset %>% select(all_of(c(studyname, covariates))) %>% group_by(across(all_of(studyname))) %>% summarise_all(~length(.))
    missingpercent <- dataset %>% select(all_of(c(studyname, covariates))) %>% group_by(across(all_of(studyname))) %>% summarise_all(~round(sum(is.na(.)/length(.))*100))
    
    numberofNAs <- dataset %>% select(all_of(c(studyname, covariates))) %>% group_by(across(all_of(studyname))) %>% summarize_all(~sum(is.na(.)))
    studysize <- dataset %>% select(all_of(c(studyname, covariates))) %>% group_by(across(all_of(studyname))) %>% summarize_all(~length(.))
    
    sys_missing <- apply(numberofNAs[,covariates] == studysize[,covariates], 2, any)
    sys_covariates <- names(sys_missing)[which(sys_missing == TRUE)]
    spor_missing <- apply(numberofNAs[,covariates], 2, sum) > 0 & !sys_missing
    spor_covariates <- names(spor_missing)[which(spor_missing == TRUE)]
    without_sys_covariates <- names(sys_missing)[which(sys_missing == FALSE)]
  }
  
  return(list(missingcount = missingcount, missingpercent = missingpercent, sys_missing = sys_missing, spor_missing = spor_missing, sys_covariates = sys_covariates, spor_covariates = spor_covariates, without_sys_covariates = without_sys_covariates, 
              covariates = covariates, typeofvar = typeofvar, studyname = studyname, treatmentname = treatmentname, outcomename = outcomename))
}



#' Find correct imputation method to be used in the mice package
#'
#' Find correct imputation method for the mice package
#'
#' @param dataset data which contains variables of interest
#' @param missingPattern missing pattern object created using \code{\link{findMissingPattern}}
#' @param interaction Indicator denoting whether treatment-covariate interactions should be included. Default is set to true.#'
#' @export

#Find correct imputation method to be used in the mice package

getCorrectMeth <- function(dataset = NULL, missingPattern = NULL, interaction = TRUE){
  
  if(is.null(dataset) | is.null(missingPattern)){
    stop("dataset and missingPattern have to be specified.")
  }
  
  with(missingPattern, {
    
    totalstudies <- dataset %>% select(studyname) %>% n_distinct()
    
    meth <- make.method(dataset)
    if(totalstudies == 1){
      
      meth[paste0(without_sys_covariates, treatmentname)] <- paste0("~ I(as.numeric(as.character(", without_sys_covariates, ")) *", treatmentname, ")")
      
      if(meth[outcomename] != ""){
        meth[outcomename] <- "pmm"
      }
      
      if(length(spor_covariates) != 0){
        meth[spor_covariates] <- "pmm"
      }
      
    } else{
      
      if(meth[outcomename] != ""){
        meth[outcomename] <- "2l.pmm" #assume outcome data is not systematically missing
      }
      
      if(length(spor_covariates) != 0){
        meth[spor_covariates] <- "2l.pmm"
      }
      
      if(length(sys_covariates) != 0){
        
        for(i in 1:length(sys_covariates)){
          if(typeofvar[sys_covariates[i]] == "continuous"){
            meth[sys_covariates[i]] <- "2l.2stage.norm"
          } else if(typeofvar[sys_covariates[i]] == "binary"){
            meth[sys_covariates[i]] <- "2l.2stage.bin"
          } else if(typeofvar[sys_covariates[i]] == "count"){
            meth[sys_covariates[i]] <- "2l.2stage.pois"
          }
        }
      }
      
      if(interaction == TRUE){
        meth[paste0(covariates, treatmentname)] <- paste0("~ I(as.numeric(as.character(", covariates, ")) *", treatmentname,  ")")
      }
    }
    return(meth)
  })

}


#' Find correct imputation prediction matrix to be used in the mice package
#'
#' Find correct imputation prediction matrix for the mice package
#'
#' @param dataset data which contains variables of interest
#' @param missingPattern missing pattern object created using \code{\link{findMissingPattern}}
#' @param interaction Indicator denoting whether treatment-covariate interactions should be included. Default is set to true.#'
#' @export

getCorrectPred <- function(dataset = NULL, missingPattern = NULL, interaction = TRUE){
  
  
  if(is.null(dataset) | is.null(missingPattern)){
    stop("dataset and missingPattern have to be specified.")
  }
  
  with(missingPattern, {

  totalstudies <- dataset %>% select(studyname) %>% n_distinct()
    
  if(totalstudies == 1){
    
    # Case when there are only one cluster/study
    pred <- make.predictorMatrix(dataset)
    pred[,] <- 0
      
    #outcome imputation
    pred[outcomename,] <- 1
    pred[outcomename, studyname] <- 0
      
    # sporadically missing predictors imputation
    if(length(spor_covariates) != 0){
      pred[spor_covariates,] <- 1
        
      for(i in 1:length(spor_covariates)){
        pred[spor_covariates[i], paste0(spor_covariates[i], treatmentname)] <- 0
      }
      pred[spor_covariates, studyname] <- 0
    }
    diag(pred) <- 0
  } else{
    
    # Case when there are multiple clusters/studies
    pred <- make.predictorMatrix(dataset)
    pred[,] <- 0
      
    # outcome imputation
    pred[outcomename,] <- 1
      
    if(treatmentname %in% colnames(pred)){
      pred[outcomename, treatmentname] <- 2
    }
    pred[outcomename, studyname] <- -2
      
    # sporadically missing predictors imputation
    if(length(spor_covariates) != 0){
      pred[spor_covariates,] <- 1
        
    if(treatmentname %in% colnames(pred)){
      pred[spor_covariates, treatmentname] <- 2
    }
        
    if(interaction == TRUE){
      for(i in 1:length(spor_covariates)){
        pred[spor_covariates[i], paste0(spor_covariates[i], treatmentname)] <- 0
      }
    }
    pred[spor_covariates, studyname] <- -2
    }
      
    # systematically missing predictors imputation
    if(length(sys_covariates) != 0){
      pred[c(outcomename, sys_covariates),] <- 1
      pred[sys_covariates, outcomename] <- 1
        
      if(treatmentname %in% colnames(pred)){
        pred[c(outcomename, sys_covariates), treatmentname] <- 2
      }
        
      if(interaction == TRUE){
        for(i in 1:length(sys_covariates)){
          pred[sys_covariates[i], paste0(sys_covariates[i], treatmentname)] <- 0
        }
      }
      pred[c(outcomename, sys_covariates), studyname] <- -2
    }
    diag(pred) <- 0
  }
    
  return(pred)
  })
}

