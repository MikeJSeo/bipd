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
#' @param missingPattern missing pattern object created using \code{\link{findMissingPattern}}.
#' @param typeofvar type of variable: continuous, binary, or count
#' @param studyname study name
#' @param treatmentname treatment name
#' @param interaction whether to include interaction or not
#'
#' @export


getCorrectMeth <- function(dataset, missingPattern, typeofvar = NULL, studyname = "study", treatmentname = "treat", interaction = FALSE){
  
  meth <- make.method(dataset)
  if(length(unique(dataset[,studyname])) == 1){
    
    meth[paste0(missingP$without_sys_covariates, treatmentname)] <- paste0("~ I(as.numeric(as.character(", missingP$without_sys_covariates, ")) *", treatname, ")")
    
    if(meth["y"] != ""){
      meth["y"] <- "pmm"
    }
    
    if(length(missingP$spor_covariates) != 0){
      meth[missingP$spor_covariates] <- "pmm"
    }
    
  } else{
    
    if(meth["y"] != ""){
      meth["y"] <- "2l.pmm" #assume outcome data is not systematically missing
    }
    
    if(length(missingP$spor_covariates) != 0){
      
      for(i in 1:length(missingP$spor_covariates)){
        meth[missingP$spor_covariates[i]] <- "2l.pmm"
      }
    }
    
    if(length(missingP$sys_covariates) != 0){
      
      for(i in 1:length(missingP$sys_covariates)){
        if(typeofvar[missingP$sys_covariates[i]] == "continuous"){
          meth[missingP$sys_covariates[i]] <- "2l.2stage.norm"
          #meth[missingP$sys_covariates[i]] <- "2l.glm.norm"
        } else if(typeofvar[missingP$sys_covariates[i]] == "binary"){
          meth[missingP$sys_covariates[i]] <- "2l.2stage.bin" 
          #meth[missingP$sys_covariates[i]] <- "2l.glm.bin"
        } else if(typeofvar[missingP$sys_covariates[i]] == "count"){
          meth[missingP$sys_covariates[i]] <- "2l.2stage.pois"
          #meth[missingP$sys_covariates[i]] <- "2l.glm.pois"
        }
      }
    }
    
    if(length(missingP$spor_covariates) != 0){
      meth[missingP$spor_covariates] <- "2l.pmm"
    }
    
    if(interaction == TRUE){
      meth[paste0(missingP$covariates, "treat")] <- paste0("~ I(as.numeric(as.character(", missingP$covariates, ")) *", treatmentname,  ")")
    }
  }
  return(meth)
}
