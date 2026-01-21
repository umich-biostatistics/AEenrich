
#' Covid Vaccine Adverse Event Data
#'
#' Adverse event data in the long format. Each row is a single adverse 
#' event, along with covariates.
#' 
#' @return A data frame with the following columns:
#' \item{VAERS_ID}{Event ID}
#' \item{VAX_LABEL}{Vaccine type}
#' \item{AE_NAME}{Adverse event name}
#' \item{AGE}{Covariate}
#' \item{SEX}{Covariate}
#'
"covid1"

#' Covid Vaccine Adverse Event Data
#'
#' Adverse event data in the short format. Each row is a count of adverse 
#' events with the given name.
#' 
#' @return A data frame with the following columns:
#' \item{DRUG_TYPE}{Vaccine type}
#' \item{AE_NAME}{Adverse event name}
#' \item{AEYes}{Number of observations that have this adverse event}
#' \item{AENo}{Number of observations that do not have this adverse event}
#' \item{AGE}{Covariate}
#' \item{SEX}{Covariate}
#'
"covid2"

#' Group Structure Data
#'
#' Identifies which group each set of adverse events belongs. 
#' 
#' @return A data frame with the following columns:
#' \item{AE_NAME}{Adverse event name}
#' \item{GROUP_NAME}{Group name}
#'
"group"
