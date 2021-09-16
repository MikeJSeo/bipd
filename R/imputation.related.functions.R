##########################################################################
#helpful functions that were used to impute missing values using MICE package
##########################################################################

#' Find missing data pattern in a given data
#'
#' Find missing data pattern in a given data i.e. whether they are systematically missing or sporadically missing
#' 
#' @param dataset data which contains variables of interest
#' @param covariates vector of variable names to find missing data pattern
#' @param studyname study name
#'
#' @export

findMissingPattern <- function(dataset, covariates, studyname = "study"){
  
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
#' @param studyname study name
#' @param treatmentname treatment name
#' @param outcomename outcome name
#' @param interaction whether to include interaction or not
#' @param typeofvar type of variables; should be a vector of these values: "continuous", "binary", or "count". Index value of the vector should be predictor names.
#'
#' @export

getCorrectMeth <- function(dataset, missingPattern, studyname = "study", treatmentname = "treat", outcomename = "y", interaction = FALSE, typeofvar = NULL){
  
  if(length(missingPattern$sys_covariates) != 0 & is.null(typeofvar)){
    stop("typeofvar needs to be specified to denote the type of the systematically missing variables")
  } 
  
  meth <- make.method(dataset)
  if(length(unique(dataset[,studyname])) == 1){
    
    meth[paste0(missingPattern$without_sys_covariates, treatmentname)] <- paste0("~ I(as.numeric(as.character(", missingPattern$without_sys_covariates, ")) *", treatname, ")")
    
    if(meth[outcomename] != ""){
      meth[outcomename] <- "pmm"
    }
    
    if(length(missingPattern$spor_covariates) != 0){
      meth[missingPattern$spor_covariates] <- "pmm"
    }
    
  } else{
    
    if(meth[outcomename] != ""){
      meth[outcomename] <- "2l.pmm" #assume outcome data is not systematically missing
    }
    
    if(length(missingPattern$spor_covariates) != 0){
      
      for(i in 1:length(missingPattern$spor_covariates)){
        meth[missingPattern$spor_covariates[i]] <- "2l.pmm"
      }
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
    
    if(length(missingPattern$spor_covariates) != 0){
      meth[missingPattern$spor_covariates] <- "2l.pmm"
    }
    
    if(interaction == TRUE){
      meth[paste0(missingPattern$covariates, "treat")] <- paste0("~ I(as.numeric(as.character(", missingPattern$covariates, ")) *", treatmentname,  ")")
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
#' @param interaction whether to include interaction or not
#'
#' @export

getCorrectPred <- function(dataset, missingPattern, studyname = "study", treatmentname = "treat", outcomename = "y", interaction = FALSE){
  
  if(length(unique(dataset[,studyname])) == 1){
    
    with(missingPattern, {
      pred <- make.predictorMatrix(dataset)
      pred[,] <- 0
      
      if(length(spor_covariates) != 0){
        pred[c(outcomename, spor_covariates),] <- 1
        diag(pred) <- 0
        
        for(i in 1:length(spor_covariates)){
          pred[spor_covariates[i], paste0(spor_covariates[i], treatmentname)] <- 0
        }
        pred[c(outcomename, spor_covariates), studyname] <- 0
      } else{
        
        pred[outcomename,] <- 1
        diag(pred) <- 0
        pred[outcomename, studyname] <- 0
      }
      pred
    })
  } else{
    
    with(missingPattern, {
      pred <- make.predictorMatrix(dataset)
      pred[,] <- 0
      
      if(length(missingPattern$sys_covariates) == 0){
        if(length(spor_covariates) != 0){
          pred[c(outcomename, spor_covariates),] <- 1
          pred[c(outcomename, spor_covariates), treatmentname] <- 2 
          diag(pred) <- 0
          
          for(i in 1:length(spor_covariates)){
            pred[spor_covariates[i], paste0(spor_covariates[i], treatmentname)] <- 0
          }
          pred[c(outcomename, spor_covariates), studyname] <- -2
        } else{
          
          pred[outcomename,] <- 1
          pred[outcomename, treatmentname] <- 2
          diag(pred) <- 0
          pred[outcomename, studyname] <- -2
        }  
      } else if(length(missingPattern$sys_covariates) != 0){
        
        if(length(spor_covariates) != 0){
          pred[c(outcomename, spor_covariates),] <- 2
          pred[spor_covariates, outcomename] <- 1
          
          if(interaction == TRUE){
            pred[c(outcomename, spor_covariates), treatmentname] <- 2   
            for(i in 1:length(spor_covariates)){
              pred[spor_covariates[i], paste0(spor_covariates[i], treatmentname)] <- 0
            }
          }
          
          pred[c(outcomename, spor_covariates), studyname] <- -2
        }
        
        if(length(sys_covariates) != 0){
          pred[c(outcomename, sys_covariates),] <- 2
          pred[sys_covariates, outcomename] <- 1
          
          if(interaction == TRUE){
            pred[c(outcomename, sys_covariates), treatmentname] <- 2 
            for(i in 1:length(sys_covariates)){
              pred[sys_covariates[i], paste0(sys_covariates[i], treatmentname)] <- 0
            }
          }
          pred[c(outcomename, sys_covariates), studyname] <- -2
        }
        diag(pred) <- 0
      }
      pred
    })
  }
}

