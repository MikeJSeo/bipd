
#' Impute missing data in individual participant data with study-level cluster.
#'
#' Impute missing data in individual participant data where data is clustered by different studies.
#' In the presence of systematically missing variables, we use 2l.2stage.norm, 2l.2stage.bin, 
#' and 2l.2stage.pois methods depending on the type of the covariates to impute. 
#' If there are no systematically missing variables, the function defaults to use 2l.pmm in miceadds package 
#' which generalizes predictive mean matching using linear mixed model.
#'
#'
#' @param dataset Data which contains variables of interests
#' @param covariates Vector of variable names to find missing data pattern
#' @param typeofvar Type of variables; should be a vector of these values: "continuous", "binary", or "count". Index value of the vector should be predictor names.
#' @param interaction Indicator denoting whether treatment-covariate interactions should be included
#' @param meth User can specify imputation method to be used in the mice package. If left unspecified, function picks a reasonable one.
#' @param pred User can specify correct prediction matrix to be used in the mice package. If left unspecified, function picks a reasonable one.
#' @param studyname Study name in the data specified. Default is "study".
#' @param treatmentname Treatment name in the data specified. Default is "treat".
#' @param outcomename Outcome name in the data specified. Default is "y".
#' @return 
#' \item{missingPattern}{missing Pattern object returned by running \code{\link{findMissingPattern}}}
#' \item{meth}{imputation method used with the mice function}
#' \item{pred}{prediction matrix used with the mice function}
#'
#' @export

ipdma.impute <- function(dataset = NULL, covariates = NULL, typeofvar = NULL, interaction = TRUE,
                         meth = NULL, pred = NULL, studyname = NULL, treatmentname = NULL, outcomename = NULL 
                         ){
  
  if(is.null(studyname) | is.null(treatmentname) | is.null(outcomename)){
    stop("studyname, treatmentname, and outcomename have to be specified.")
  }
  
  if(is.null(dataset) | is.null(covariates) | is.null(typeofvar)){
    stop("dataset, covariates, and typeofvar have to be specified.")
  }
  
  dataset <- dataset[,c(studyname, treatmentname, outcomename, covariates)]
  
  missingPattern <- findMissingPattern(dataset = dataset, covariates = covariates, studyname = studyname)
  
  meth = getCorrectMeth(dataset = dataset, missingPattern = missingPattern, typeofvar = typeofvar, interaction = TRUE,
                        studyname = studyname, treatmentname = treatmentname, outcomename = outcomename)
  
  list(missingPattern = missingPattern, meth = meth)
  
}
  

#' Find missing data pattern in a given data
#'
#' Find missing data pattern in a given data i.e. whether they are systematically missing or sporadically missing. Also calculates missing count and percentage for exploratory analysis.
#'
#' @param dataset Data which contains variables of interests
#' @param covariates Vector of variable names that the user is interested in finding a missing data pattern
#' @param studyname Study name in the data specified. Default is "study"
#'
#' @export

findMissingPattern <- function(dataset = NULL, covariates = NULL, studyname = NULL){
  
  if(is.null(dataset) | is.null(covariates) | is.null(studyname)){
    stop("All parameters dataset, covariates, and studyname Dataset is not given")
  }

  if(length(unique(dataset[,studyname])) == 1){
    
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
  
  return(list(missingcount = missingcount, missingpercent = missingpercent, sys_missing = sys_missing, spor_missing = spor_missing, sys_covariates = sys_covariates, spor_covariates = spor_covariates, without_sys_covariates = without_sys_covariates, covariates = covariates))
}



#' Find correct imputation method to be used in the mice package
#'
#' Find correct imputation method for the mice package
#'
#' @param dataset data which contains variables of interest
#' @param missingPattern missing pattern object created using \code{\link{findMissingPattern}}
#' @param typeofvar type of variables; should be a vector of these values: "continuous", "binary", or "count". Index value of the vector should be predictor names.
#' @param interaction indicator for including covariate-treatment interactions
#' @param studyname study name
#' @param treatmentname treatment name
#' @param outcomename outcome name
#'
#' @export

#Find correct imputation method to be used in the mice package

getCorrectMeth <- function(dataset = NULL, missingPattern = NULL, typeofvar = NULL, interaction = TRUE, studyname = NULL, treatmentname = NULL, outcomename = NULL){
  
  if(is.null(studyname) | is.null(treatmentname) | is.null(outcomename)){
    stop("studyname, treatmentname, and outcomename have to be specified.")
  }
  
  if(is.null(dataset) | is.null(covariates) | is.null(typeofvar)){
    stop("dataset, missingPattern, and typeofvar have to be specified.")
  }
  
  
  meth <- make.method(dataset)
  if(length(unique(dataset[,studyname])) == 1){
    
    meth[paste0(missingPattern$without_sys_covariates, treatmentname)] <- paste0("~ I(as.numeric(as.character(", missingPattern$without_sys_covariates, ")) *", treatmentname, ")")
    
    if(meth[outcomename] != ""){
      meth[outcomename] <- "pmm"
    }
    
    if(length(missingPattern$spor_covariates) != 0){
      meth[missingPattern$spor_covariates] <- "pmm"
    }
    
  } else{
    
    if(length(missingPattern$sys_covariates) != 0 & is.null(typeofvar)){
      stop("typeofvar needs to be specified to denote the type of the systematically missing variables")
    } 
    
    if(meth[outcomename] != ""){
      meth[outcomename] <- "2l.pmm" #assume outcome data is not systematically missing
    }
    
    if(length(missingPattern$spor_covariates) != 0){
      meth[missingPattern$spor_covariates] <- "2l.pmm"
    }
    
    if(length(missingPattern$sys_covariates) != 0){
      
      for(i in 1:length(missingPattern$sys_covariates)){
        if(typeofvar[missingPattern$sys_covariates[i]] == "continuous"){
          meth[missingPattern$sys_covariates[i]] <- "2l.2stage.norm"
        } else if(typeofvar[missingPattern$sys_covariates[i]] == "binary"){
          meth[missingPattern$sys_covariates[i]] <- "2l.2stage.bin"
        } else if(typeofvar[missingPattern$sys_covariates[i]] == "count"){
          meth[missingPattern$sys_covariates[i]] <- "2l.2stage.pois"
        }
      }
    }
    
    if(interaction == TRUE){
      meth[paste0(missingPattern$covariates, treatmentname)] <- paste0("~ I(as.numeric(as.character(", missingPattern$covariates, ")) *", treatmentname,  ")")
    }
  }
  return(meth)
}


#' Find correct imputation prediction matrix to be used in the mice package
#'
#' Find correct imputation prediction matrix for the mice package
#'
#' @param dataset data which contains variables of interest
#' @param missingPattern missing pattern object created using \code{\link{findMissingPattern}}
#' @param studyname study name
#' @param treatmentname treatment name
#' @param outcomename outcome name
#' @param interaction indicator for including covariate-treatment interactions
#'
#' @export

getCorrectPred <- function(dataset = NULL, missingPattern = NULL, studyname = NULL, treatmentname = NULL, outcomename = NULL, interaction = TRUE){
  
  
  if(is.null(studyname) | is.null(treatmentname) | is.null(outcomename)){
    stop("studyname, treatmentname, and outcomename have to be specified.")
  }
  
  if(is.null(dataset) | is.null(missingPattern)){
    stop("dataset and missingPattern have to be specified.")
  }
  
  if(length(unique(dataset[,studyname])) == 1){
    
    # Case when there are only one cluster/study
    with(missingPattern, {
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
      pred
    })
  } else{
    
    # Case when there are multiple clusters/studies
    with(missingPattern, {
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
      pred
    })
  }
}



####################### some other helpful functions
createinteractions <- function(dataset, cov) {
  for(i in 1:length(cov)){
    varname <- paste0(cov[i], "treat")
    df[[varname]] <- NA  
  }
  df
}


findVarianceUsingRubinsRule <- function(prediction.dummy, variance.dummy){
  
  avg.prediction <- apply(prediction.dummy, 1, mean)
  summ <- 0
  for(iii in 1:5){
    summ <- summ + (prediction.dummy[,iii] - avg.prediction)^2
  }
  return(apply(variance.dummy, 1, mean) + (5 + 1)/ (5^2 - 5) * summ)
}


