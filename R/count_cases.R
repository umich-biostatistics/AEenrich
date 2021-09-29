
#' The count_cases function is used to convert data on the report level to
#' aggregated data, grouping by specified covariates.
#' 
#' @param data a data.frame with at least 3 columns, consisting data on the report
#' level, having ID, Drug type and AE name as the first 3 columns with
#' covariates(optional) followed. The order of columns is not interchangeable.
#' @param drug.case a character string for the target drug of interest.
#' @param drug.control a character string for the reference drug. If NULL(default),
#' all other drugs combined are the reference.
#' @param covar_cont a character vector of continuous covariates.
#' @param covar_disc a character vector of categorical covariates.
#' @param breaks a list consists of vectors used for creating specific bins to
#' transform continuous covariates into categorical. Breaks Should have the same
#' length as covar_cont. Given a vector of non-decreasing breakpoints in `breaks[i]`,
#' find the interval containing each element of `covar_cont[i]`; i.e., for each index
#' j in `breaks[i]`, value j is assigned to `covar_cont[i]` if and only if `breaks[i][j]`
#' `\leq covar_cont[i] < breaks[i][j+1]`.
#' @param min_AE the minimum number of cases required to start counting
#' for a specific AE. Default 10.
#' @param cores the number of cores to use for parallel execution.
#' 
#' @return A **data.frame** consists of aggregated data.
#' 
#' The returned data.frame contains the following columns:
#' \itemize{
#'   \item{DRUG_TYPE: }{type of the drug, DrugYes for target drug and DrugNo for referenced drug}
#'   \item{AE_NAME: }{the name of the adverse event}
#'   \item{AEyes: }{number of observations that have this AE}
#'   \item{AEno: }{number of observations that do not have this AE}
#'   \item{covariates: }{covariates specifed by user}
#' }
#' 
#' @examples
#' 
#' 
#' # count_cases(data = covid1, drug.case = "COVID19", drug.control = "OTHER",
#' #             covar_cont = c("AGE"), covar_disc = c("SEX"),
#' #             breaks = list(c(16,30,50,65,120)))
#'             

count_cases = function(data, drug.case = drug.case, drug.control = NULL,
                       covar_disc = NULL, covar_cont = NULL, breaks = NULL,
                       cores = detectCores(), min_AE = 10){
  . <- "Muted"
  data = as_tibble(data) 
  ## check the breaks-covar_cont pairs
  if(!length(breaks) == length(covar_cont)){
    stop("The length of breaks does not match that of continuous covariates")
  }
  ## check the basic first three columns
  if(length(names(data)) < 3){
    stop("Unexpected data type")
  }
  ## check the existence of continuous covariates
  if(!is.null(covar_cont)){
    if(!all(covar_cont %in% names(data))){
      stop("nonexistent column") } else {
        if(!all(sapply(data, is.numeric)[covar_cont])){
          stop("Discrete covariates misclassified as continuous")
        }
      }
  }
  ## check the existence of discrete covariates
  if(!is.null(covar_disc)){
    if(!all(covar_disc %in% names(data))){
      stop("nonexistent column") }else {
        if(any(sapply(data, is.numeric)[covar_disc])){
          stop("Continuous covariates misclassified as discrete")
        }
      }
  }
  ## rename the first three columns
  names(data)[1:3] = c('ID', 'DRUG_TYPE', 'AE_NAME')
  ## filter by drug.case and drug.control
  if(!is.null(drug.control)){
    drug_list = c(drug.case, drug.control)
    data = data[data$DRUG_TYPE %in% drug_list, ]
    }
  ## remove NAs rename drug type based on drug.case
  data = data[complete.cases(data), ] %>%
    mutate(DRUG_TYPE = ifelse(DRUG_TYPE %in% drug.case,
                              "DrugYes",
                              "DrugNo") )
  
  ## filter out reports with both case and control vaccines
  ID_No = data %>%
    filter(DRUG_TYPE == "DrugNo") %>%
    .$ID %>%
    unique()
  ID_Yes = data %>%
    filter(DRUG_TYPE == "DrugYes") %>%
    .$ID %>%
    unique()
  Confused_ID = intersect(ID_No, ID_Yes)
  
  ## convert all the character variable to factor
  data_temp = data %>%
    filter(!ID %in% Confused_ID)
  
  ## for continuous covariates, classifying each obs by argument 'breaks'
  if(!is.null(covar_cont)){
    for(i in 1:length(breaks)){
      breaks_temp = breaks[[i]]
      var_temp = covar_cont[i]
      covar_group = data_temp %>%
        dplyr::select(all_of(var_temp)) %>%
        .[[1]] %>%
        findInterval(breaks_temp)
      data_temp = data_temp %>%
        dplyr::select(-all_of(var_temp)) %>%
        tibble({{var_temp}} := covar_group) %>%
        filter(!!as.symbol(var_temp) != 0) %>%
        filter(!!as.symbol(var_temp) != length(breaks_temp)) %>%
        mutate(!!as.symbol(var_temp) := as.factor(!!as.symbol(var_temp)))
      }
  }
  ## filter out AE with less than 10 observations
  AE_list = data_temp %>%
    group_by(AE_NAME) %>%
    summarise(count = n()) %>%
    filter(count >= as.integer(min_AE)) %>%
    .$AE_NAME %>%
    as.character()
  data_comp = data_temp %>%
    filter(AE_NAME %in% AE_list) %>%
    mutate_if(sapply(data, is.character), as.factor)

  cl = makeCluster(cores)
  registerDoParallel(cl)
  results = foreach(i = 1:length(AE_list),
                    .packages = c("tidyverse"),
                    .combine = bind_rows 
                    ) %dopar% {
    AE = AE_list[i]
    AE_yes = data_comp %>%
      filter(AE_NAME == AE) %>%
      mutate(AE_NAME = "AEYes") %>%
      distinct(ID, .keep_all = TRUE)
    ## filter by ID
    AE_no = data_comp %>%
      filter(! ID %in% AE_yes$ID) %>%
      mutate(AE_NAME = "AENo") %>%
      distinct(ID, .keep_all = TRUE)
    
    ## combine together
    data_AE = AE_yes %>%
      bind_rows(AE_no) %>%
      mutate(AE_NAME = as.factor(AE_NAME))
    ## grouping variables
    grp_cols = c("DRUG_TYPE", "AE_NAME", covar_cont, covar_disc)
    ## convertnto count
    data_count = data_AE %>%
      dplyr::select(-ID) %>%
      group_by_at(grp_cols) %>%
      summarise(count = n(), .groups = "drop") %>%
      pivot_wider(names_from = AE_NAME, values_from = count) %>%
      mutate(AEYes = replace_na(AEYes, 0),
             AENo = replace_na(AENo, 0)) %>%
      mutate(AE_NAME = {{AE}}) %>%
      relocate(DRUG_TYPE, AE_NAME, AEYes, AENo)
    data_count
    }
  stopCluster(cl)
  return(results)
}

#' @description The Count_cases function is used to convert data on the report
#' level to aggregated data, grouping by specified covariates.
#' 
#' Use the function `count_cases` to convert report level data into aggregated data.
#' 
#' See our \href{https://github.com/umich-biostatistics/AEenrich}{Github home page} 
#' or run ?count_cases for examples.
"_PACKAGE"
