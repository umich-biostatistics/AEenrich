## Function odds_ratio: -------------------------------------------------------
## Estimate log odds ratio for each AE, perform Logistic regression on data, 
## use Logit as link function.
###
###  Input: 
###     1. data: 
###           Type I:  ID   DRUG_TYPE   AE_NAME    AGE   SEX
###                   201     FLUN      Insomnia   69     F
###                   ...     ...       ...        ...   ...
###                   299     FLU       Chills     68     M
###
###           Type II: DRUG_TYPE   AE_NAME   COUNT(YES)  COUNT(NO)  AGE   SEX
###                     FLUN      Insomnia   640          6544      69     F
###                     ...        ...       ...          ...       ...   ...
###                     FLU       Chills     586          3720      68     M
###
###       For Type II data, the 3rd and 4th Columns give the numbers of
###         successes(have AE) and failures(Do not have AE) respectively.       
###
###     2. drug.case: Target vaccine type
###     3. drug.control: Reference vaccine type
###     4. covar: Covariates for logistic regression
###     5. min_AE: The minimum number of cases required to start counting
###        for a specific AE.
###     6. cores: The number of cores to use for parallel execution.
###  Output: 
###       a tibble(data.frame) with six columns:
### 
###         AE_NAME                   BETA   p_value   `95Lower`
###         <chr>                     <dbl>   <dbl>     <dbl>
###         Injection site joint pain -1.58  4.75e- 18  ...
###         Apathy                    -0.166 7.02e-  1  ...
###         Arthralgia                 0.225 2.07e- 24  ...
###         
###         `95Upper`   `se(logOR)`
###          <dbl>        <dbl>
###           ...         ...
###           ...         ...
###           ...         ...
###     AE_NAME: The name of the adverse event.
###     BETA: Log odds ratio of AEs.
###     p_value: p value of log odds.
###     `95Lower` and `95Upper`: The lower/upper bound of confidence
###     interval of odds ratio.
# 79: -------------------------------------------------------------------------
odds_ratio = function(data, drug.case, drug.control = NULL, covar = NULL,
                      min_AE = 10, cores = detectCores()){
  i <- "Muted"
  . <- "Muted"
  data = as_tibble(data)
  if(!is.null(covar)){
    if(!all(covar %in% names(data))){
      stop("covariates not found") }
  }
  
  ## Check data type
  if (!sapply(data, is.numeric)[3]){
    ## change the names of columns 
    names(data)[1:3] = c('ID', 'DRUG_TYPE', 'AE_NAME')
    if(!is.null(drug.control)){
      drug_list = c(drug.case, drug.control)
      data = data[data$DRUG_TYPE %in% drug_list, ]
      }
    ## remove NA  
    data = data[complete.cases(data), ]
    data = data %>%
      mutate(DRUG_TYPE = ifelse(DRUG_TYPE %in% drug.case, "DrugYes", "DrugNo") )
    ## filter out AE with less than 10 observations
    AE_list = data %>%
      group_by(AE_NAME) %>%
      summarise(count = n()) %>%
      filter(count >= as.integer(min_AE)) %>%
      .$AE_NAME %>%
      as.character()
    ## Convert character columns to factor
    data_temp = data %>%
      filter(AE_NAME %in% AE_list) %>%
      mutate_if(sapply(data, is.character), as.factor)
    
    data_comp = data_temp
    ## A set consists of unique AE names
    AE_SET = unique(as.character(data_comp$AE_NAME))
    
    cl = makeCluster(cores)
    registerDoParallel(cl)
    results = foreach(i = 1:length(AE_SET),
                      .packages = c("tidyverse"),
                      .combine = bind_rows 
                      ) %dopar% {
        AE = AE_SET[i]
        ## those have AE
        AE_yes = data_comp %>%
          filter(AE_NAME == AE) %>%
          mutate(AE_NAME = "AEYes") %>%
          distinct(ID, .keep_all = TRUE)
        ## those who don't
        ID_list = AE_yes$ID
        ## filter by ID
        AE_no = data_comp %>%
          filter(! ID %in% ID_list) %>%
          mutate(AE_NAME = "AENo") %>%
          distinct(ID, .keep_all = TRUE)
        ## combine data together
        data_AE = AE_yes %>%
          bind_rows(AE_no) %>%
          mutate(AE_NAME = as.factor(AE_NAME))
      
        covar_formula = ifelse(is.null(covar), "",
                              paste("+", paste(covar, collapse = " + ")))
        str_formula = as.formula(paste("AE_NAME ~ DRUG_TYPE", covar_formula))
        ## Logistic regression
        mod = glm(formula = str_formula, family = binomial(link = logit),
                  data = data_AE)
        ## return the log odds ratio
        tibble(AE_NAME = AE,
               BETA = coef(mod)["DRUG_TYPEDrugYes"],
               p_value = coef(summary(mod))[,4]["DRUG_TYPEDrugYes"],
               `95Lower` = exp(confint.default(mod)[2,])[1],
               `95Upper` = exp(confint.default(mod)[2,])[2],
               `se(logOR)` = coef(summary(mod))[,2]["DRUG_TYPEDrugYes"])
                
        }
    stopCluster(cl)
    } else {
      if (!sapply(data, is.numeric)[4]){
        stop("Invalid data type")
        }
      ## change the names of columns 
      names(data)[1:4] = c('DRUG_TYPE', 'AE_NAME', 'YES', 'NO')
      if(!is.null(drug.control)){
        drug_list = c(drug.case, drug.control)
        data = data[data$DRUG_TYPE %in% drug_list, ]
      }
      ## remove NA  
      data = data[complete.cases(data), ]
      data = data %>%
        mutate(DRUG_TYPE = ifelse(DRUG_TYPE %in% drug.case,
                                  "DrugYes",
                                  "DrugNo") )
      ## Type two data, so every covariate should be factor
      if(length(names(data)) == 4){
        index = 1:2
      } else{
          index = c(1:2, 5:length(names(data))) }
      data_temp = data %>%
        mutate_at(.vars = index, as.factor)
      # filter out AEs with less than 10 observations
      AE_list = data_temp %>%
        group_by(AE_NAME) %>%
        summarize(count = sum(YES)) %>%
        filter(count >= as.integer(min_AE)) %>%
        .$AE_NAME %>%
        as.character()
      data_comp = data_temp %>%
        filter(AE_NAME %in% AE_list)
      ## A set consists of AE names
      AE_SET = unique(as.character(data_comp$AE_NAME))
      
      cl = makeCluster(cores)
      registerDoParallel(cl)
      results = foreach(i = 1:length(AE_SET),
                        .packages = c("tidyverse"),
                        .combine = bind_rows 
      ) %dopar% {
        AE = AE_SET[i]
        data_count = data_comp %>%
          filter(AE_NAME == {{AE}}) %>%
          mutate(YES = YES + 1,
                 NO = NO + 1)
      
        covar_formula = ifelse(is.null(covar), "",
                               paste("+", paste(covar, collapse = " + ")))
      
        str_formula = as.formula(paste("cbind(YES, NO) ~ DRUG_TYPE", covar_formula))
        ## Logistic regression
        mod = glm(formula = str_formula, family = binomial(link = logit),
                  data = data_count)
        ## return the log odds ratio
        tibble(AE_NAME = AE,
               BETA = coef(mod)["DRUG_TYPEDrugYes"],
               p_value = coef(summary(mod))[,4]["DRUG_TYPEDrugYes"],
               `95Lower` = exp(confint.default(mod)[2,])[1],
               `95Upper` = exp(confint.default(mod)[2,])[2],
               `se(logOR)` = coef(summary(mod))[,2]["DRUG_TYPEDrugYes"])
        }
      stopCluster(cl)
      }
  return(list(res = results, ae = AE_list))
}
