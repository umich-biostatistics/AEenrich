## Function get_ES: -----------------------------------------------------------
## Calculate Enrichment Score over AE Groups, called by KS_enrichment function.
###
###  Input: 
###     1. dd.group (data.frame): dd.meddra (AE_NAME GROUP_NAME)
###     2. data_Ratio (data.frame):  AE_NAME                           OR
###                                 <chr>                            <dbl>
###                                 Injected limb mobility decreased 0.0914
###                                 Injection site joint pain        0.206 
###                                 ...                               ... 
###                                 Arthralgia                       1.25  
###     3. p: An exponent p to control the weight of the step. Default 0, which
###          corresponds to standard Kolmogorov-Smirnov statistic.
###     4. zero: logical, if TRUE, calculate zero inflated KS score. If FALSE,
###         calculate KS score without zero indicator.
###
###  Output: 
###     a data.frame:  group            ES
###                   Respiratory      0.83
###                     ...             ...
###                   Infections       0.46
# 79: -------------------------------------------------------------------------
get_ES = function(dd.group, data_Ratio, p, zero){
  . <- "Muted"
  # add zero to AE if it was not mentioned with the target vaccine
  ddF = data_Ratio %>%
    right_join(dd.group, by = "AE_NAME") %>%
    mutate(OR = coalesce(OR, 0)) %>%
    arrange(desc(OR)) # order by odds ratio
  
  # check if there are 0's
  flag_0 = any(ddF$OR == 0)

  # get interesting group names
  group.enrich = ddF %>%
    .$GROUP_NAME %>%
    unique()
  ng = length(group.enrich)
  
  get_score = function(ddF, j, p, flag_0, zero){
    # check which AE in this group
    hit_ind = ddF$GROUP_NAME == group.enrich[j]
    AE_vec = (ddF$AE_NAME)[hit_ind]
    # get the miss index (handle one AE with multiple groups) 
    miss_ind = (hit_ind == FALSE & !(ddF$AE_NAME %in% AE_vec))
    
    h_m_lst = HitMiss_Curve(ddF, miss_ind = miss_ind, p)
    h_vec = h_m_lst$hit
    m_vec = h_m_lst$miss
    position = h_m_lst$pos
    n_pos = length(position)
    
    if ((flag_0) & (zero == TRUE)){
      n_pos = length(position)
      ES = max((h_vec-m_vec)[1:(n_pos-1)])*((h_vec-m_vec)[n_pos-1]>=0)
    }else{
      ES = max(h_vec-m_vec)
    }
    ES
  }
  
  ES_vec = sapply(1:ng, function(x) get_score(ddF, x, p, flag_0, zero))
  return(data.frame(group = group.enrich, ES = ES_vec))
}
