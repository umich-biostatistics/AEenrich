## Function KS_enrichment: ----------------------------------------------------
## perform Enrichment analysis using proposed KS score and permutation test
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
###     2. dd.group (data.frame): dd.meddra (AE_NAME GROUP_NAME)
###     3. drug.case: Name of the target vaccine. 
###     4. drug.control: Name of the Reference vaccine. 
###     5. n_perms: Number of permutations to get null distribution of ES.
###     6. p: An exponent p to control the weight of the step.
###     7. zero: logical, if TRUE, calculate zero inflated KS score. If FALSE,
###         calculate KS score without zero indicator.
###     8. min_size: The minimum size of group required for enrichment analysis.
###     9. min_AE: The minimum number of cases required to start counting
###        for a specific AE.
###     10. cores: The number of cores to use for parallel execution.
###  Output:
###       a list with 2 data.frames:
###            1. Final_result (data.frame):
###                                     GROUP_NAME        ES p_value GROUP_SIZE
### 1                          Acid-base disorders 0.0000000   1.000          2
### 2                          Allergic conditions 0.0000000   1.000         41
### 3 Anaemias nonhaemolytic and marrow depression 0.6874658   0.300          1
### 4                  Ancillary infectious topics 0.4867987   0.206          3
### 5                     Angioedema and urticaria 0.0000000   1.000         17
### 6               Anxiety disorders and symptoms 0.0000000   1.000         11
### 
###            2. AE_info (data.frame):
###                                   AE_NAME               RR    p_value L U
###                                 1 Abdomen scan normal  0.25    0.428
###                                 2 Abdominal discomfort 0.137   0.517
###                                 3 Abdominal distension 0.269   0.065
###                                 4 Abdominal mass       1       0.137
###                                 5 Abdominal pain       0.282   0
###                                 6 Abdominal pain lower 0.25    0.216
# 79: -------------------------------------------------------------------------
KS_enrichment = function(data, drug.case, drug.control = NULL, covar = NULL,
                         dd.group, n_perms = 1000, p = 0, zero = FALSE,
                         min_size = 5, min_AE = 10, cores = detectCores()){
  i <- "Muted"
  ## check p is between 0 and 1
  if(p > 1 | p < 0){
    stop("P should take any value between 0 and 1")
  }
  n_perms = as.integer(n_perms)
  ## Calculate the log odds ratio.
  odds_out = odds_ratio(data, drug.case, drug.control, covar = covar,
                        min_AE = min_AE, cores = cores)
  ## Convert log odds ratio to odds ratio.
  OR_data = odds_out$res %>%
    mutate(OR = exp(BETA)) %>%
    select(AE_NAME, OR, p_value, `95Lower`, `95Upper`, `se(logOR)`)
  # Get the group size for groups of interest
  dd.group = dd.group[!duplicated(dd.group), ] %>%
    filter(!is.na(GROUP_NAME))
  df_size = dd.group %>%
    filter(AE_NAME %in% odds_out$ae) %>%
    group_by(GROUP_NAME) %>%
    summarise(GROUP_SIZE = n(), .groups = 'drop_last') %>%
    filter(GROUP_SIZE >= as.integer(min_size))
  
  ## filter out groups in which the number of AEs is less than the minimum size
  dd.group = dd.group %>%
    filter(GROUP_NAME %in% df_size$GROUP_NAME) %>%
    filter(AE_NAME %in% odds_out$ae)
  ## calculate the enrichment score(generalized KS statistic).
  ks_raw = get_ES(dd.group, OR_data[, 1:2], p, zero)
  ## Change the ES from long to wide.
  ks_true = ks_raw %>%
    as_tibble() %>%
    pivot_wider(names_from = group,
                values_from = ES)
  ## get ready for permutation test.
  perms = modelr::permute(OR_data, n_perms, OR)
  ## calculate ES for each permutation.
  perm_object = perms$perm
  cl = makeCluster(cores)
  registerDoParallel(cl)
  models = foreach(i = 1:length(perm_object),
                 .packages = c("tidyverse"),
                 .export = c("get_ES", "HitMiss_Curve")
                 ) %dopar% {
    dat = perm_object[i] %>%
      as.data.frame()
    get_ES(dd.group, dat[,1:2], p, zero)
    }
  stopCluster(cl)
  ## combine all the results of permutations together into a data frame.
  cl = makeCluster(cores)
  registerDoParallel(cl)
  ks_null = foreach(i = 1:length(models),
                  .packages = c("tidyverse"),
                  .combine = bind_rows 
                  ) %dopar% {
                    ES_i = models[[i]] %>%
                      as_tibble() %>%
                      pivot_wider(names_from = group,
                                  values_from = ES)
                    ES_i
                    }
  stopCluster(cl)
  ## add true ES to the end of the data frame.
  ks_all = ks_true %>%
    bind_rows(ks_null)
  ## calculate the p value based on permutation results.
  p_value_ks = sapply(ks_all, function(x) mean(x[2:n_perms+1] >= x[1]))
  Final_result_ks = tibble(GROUP_NAME = ks_raw$group,
                           ES = ks_raw$ES,
                           p_value = p_value_ks)
  Final_result_ks = Final_result_ks %>%
    left_join(df_size, by = "GROUP_NAME")
  return(list(Final_result = Final_result_ks, AE_info = OR_data))
}
