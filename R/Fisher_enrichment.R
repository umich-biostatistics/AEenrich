## Function Fisher_enrichment: ------------------------------------------------
## performe Enrichment test using Fisher test
###  Input: 
###     1. data (data.frame):
###           Type I: ID   DRUG_TYPE   AE_NAME
###                   201    FLUN     Insomnia
###                   ...    ...       ...
###                   299    FLU      Chills
###
###     2. dd.group (data.frame): dd.meddra (AE_NAME GROUP_NAME)
###
###     3. drug.case: Name of the Target vaccine 
###     4. drug.control: Name of the Reference vaccine 
###     5. n_perms: Number of simulations to get null distribution of ES
###     6. q.cut (numerical value): q value cut deciding the significance of
###        each AE
###     7. or.cut (numerical value): odds ratio cut deciding the significance
###        of each AE
###     8. zero: Default False, perform classic fisher exact test. If
###        True, add zero indicator to the Enrichment score.
###     9. min_size: The minimum size of group required for enrichment analysis.
###     10. min_AE: The minimum number of cases required to start counting
###        for a specific AE.
###     11. cores: The number of cores to use for parallel execution.
###  Output:
###     1. a list with 2 data.frames:
###         1. Final_result (data.frame):
###           GROUP_NAME                      ES p_value GROUP_SIZE 
###           Acid-base disorders             0       1          2    
###           Allergic conditions             0       1         41    
###           ...                            ...     ...        ...
###           Angioedema and urticaria        0       1         17    
###           Anxiety disorders and symptoms  0       1         11    
### 
###            2. AE_info (data.frame):
###                                         AE_NAME       OR    p_value   
###                                    Abdominal mass    2.13   0.175    
###                                         ...           ...     ...
###                                    Abdominal pain     0        1       
###                                         `95Lower`   `95Upper`   `se(logOR)`
###                                           2.0           2.2       ... 
###                                           ...           ...       ...
###                                            0              0       ...
#------------------------------------------------------------------------------

Fisher_enrichment = function(data, dd.group, drug.case, drug.control = NULL,
                             n_perms = 1000, q.cut = 0.05, or.cut = 1.5,
                             zero = FALSE, min_size = 5, covar = NULL,
                             min_AE = 10, cores = detectCores()){
  # Get odds ratio estimate and related info for null simulation
  odds_out = odds_ratio(data, drug.case, drug.control, covar = covar,
                        min_AE = min_AE, cores = cores)
  fisher_res =  odds_out$res %>%
    mutate(OR = exp(BETA)) %>%
    select(-BETA) %>%
    mutate(isRatio0 = as.logical(ifelse(OR < 1e-5, "TRUE", "FALSE"))) %>%
    relocate(AE_NAME, OR, p_value, isRatio0, `95Lower`, `95Upper`, `se(logOR)`)
  
  fisher_res_withQ = fisher_res %>% 
    select(-`95Lower`, -`95Upper`, -`se(logOR)`) %>%
    mutate(qval = (qvalue(p_value))$qvalue)

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
  # Calculate the Enrichment Score for each group
  fisher_result = fisher_test(dd.group, fisher_res_withQ, q.cut, or.cut,
                              n_perms, zero)
  
  Final_result_fisher = fisher_result %>%
    merge(df_size, by = 'GROUP_NAME') %>%
    as_tibble()
  return(list(Final_result = Final_result_fisher, AE_info = fisher_res[,-4]))
}
