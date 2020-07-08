
#' Flu Vaccine Adverse Event Data
#'
#' Adverse event data in the long format. Each row is a single adverse 
#' event.
#' 
#' \itemize{
#'   \item VAERS_ID Event ID
#'   \item VAX_TYPE Vaccine type
#'   \item AE_NAME Adverse event name
#' }
#'
"flu1"

#' Flu Vaccine Adverse Event Data
#'
#' Adverse event data in the short format. Each row is a count of adverse 
#' events with the given name.
#' 
#' \itemize{
#'   \item VAX_TYPE Vaccine type
#'   \item AE_NAME Adverse event name
#'   \item Count Frequency of adverse event
#' }
#'
"flu2"

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
