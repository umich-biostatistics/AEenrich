
#' Covid Vaccine Adverse Event Data
#'
#' Adverse event data in the long format. Each row is a single adverse 
#' event, along with covariates.
#' 
#' \itemize{
#'   \item VAERS_ID Event ID
#'   \item VAX_LABEL Vaccine type
#'   \item AE_NAME Adverse event name
#'   \item AGE covariate
#'   \item SEX covariate
#' }
#'
"covid1"

#' Covid Vaccine Adverse Event Data
#'
#' Adverse event data in the short format. Each row is a count of adverse 
#' events with the given name.
#' 
#' \itemize{
#'   \item DRUG_TYPE Vaccine type
#'   \item AE_NAME Adverse event name
#'   \item AEYes Number of observations that have this AE
#'   \item AENo Number of observations that do not have this AE
#'   \item AGE covariate
#'   \item SEX covariate
#' }
#'
"covid2"

#' Group Structure Data
#'
#' Identifies which group each set of adverse events belongs. 
#' 
#' \itemize{
#'   \item AE_NAME Adverse event name
#'   \item GROUP_NAME Group name
#' }
#'
"group"
