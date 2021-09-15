##########################################################################
#helpful functions that were used to impute missing values using MICE package
##########################################################################

#' Find missing data pattern, i.e. whether they are systematically missing or sporadically missing
#'
#' Find missing data pattern in a given dataset
#' 
#' @param dummydata Dataset which contains covariates of interest
#' @param covariates A vector of covariate names to find missing data pattern
#'
#' @export

findMissingPattern <- function(dummydata, covariates, studyname = "study"){
  
  if(length(unique(dummydata[,studyname])) == 1){
    
    missingcount <- dummydata %>% select(all_of(covariates))%>% summarize_all(~sum(is.na(.)))
    missingn <- dummydata %>% select(all_of(covariates))  %>% summarise_all(~length(.))
    missingpercent <- dummydata %>% select(all_of(covariates)) %>% summarise_all(~round(sum(is.na(.)/length(.))*100))
    
    numberofNAs <- dummydata %>% select(all_of(covariates)) %>% summarize_all(~sum(is.na(.)))
    studysize = dim(dummydata)[1]
    sys_missing <- numberofNAs == studysize
    spor_missing <- numberofNAs < studysize & numberofNAs > 0
    sys_covariates <- colnames(sys_missing)[which(sys_missing == TRUE)]
    spor_covariates <- colnames(spor_missing)[which(spor_missing == TRUE)]
    without_sys_covariates <- colnames(sys_missing)[which(sys_missing == FALSE)]
    
  } else{

    missingcount <- dummydata %>% select(all_of(c(studyname, covariates))) %>% group_by(across(all_of(studyname))) %>% summarize_all(~sum(is.na(.)))
    missingn <- dummydata %>% select(all_of(c(studyname, covariates))) %>% group_by(across(all_of(studyname))) %>% summarise_all(~length(.))
    missingpercent <- dummydata %>% select(all_of(c(studyname, covariates))) %>% group_by(across(all_of(studyname))) %>% summarise_all(~round(sum(is.na(.)/length(.))*100))

    numberofNAs <- dummydata %>% select(all_of(c(studyname, covariates))) %>% group_by(across(all_of(studyname))) %>% summarize_all(~sum(is.na(.)))
    studysize <- dummydata %>% select(all_of(c(studyname, covariates))) %>% group_by(across(all_of(studyname))) %>% summarize_all(~length(.))

    sys_missing <- apply(numberofNAs[,covariates] == studysize[,covariates], 2, any)
    sys_covariates <- names(sys_missing)[which(sys_missing == TRUE)]

    spor_missing <- apply(numberofNAs[,covariates], 2, sum) > 0 & !sys_missing
    spor_covariates <- names(spor_missing)[which(spor_missing == TRUE)]
    without_sys_covariates <- names(sys_missing)[which(sys_missing == FALSE)]
  }

  return(list(missingcount = missingcount, missingpercent = missingpercent, sys_missing = sys_missing, spor_missing = spor_missing, sys_covariates = sys_covariates, spor_covariates = spor_covariates, without_sys_covariates = without_sys_covariates, covariates = covariates))
}

