## Function HitMiss_Curve: ----------------------------------------------------
## Calculate miss and hit value in the KS method(considering ties),  called by
## get_ES function.
###
###  Input: 
###     1. ddF: a data frame
###         AE_NAME                 OR      GROUP_NAME                     
###           <chr>               <dbl>     <chr>                          
###         Dysgeusia              5.30     Tongue conditions              
###         ...                    ...      ...
###         Tongue disorder        4.15     Tongue conditions              
###         Paraesthesia oral      4.09     Oral soft tissue conditions    
###         Paraesthesia oral      4.09     Neurological disorders NEC  
###
###     2. miss_ind (vector): if miss, 1; Otherwise, 0
###     3. p: An exponent p to control the weight of the step.
###
###  Output:
###     a list with three vectors: P_hit, P_miss and position
# 79: -------------------------------------------------------------------------
HitMiss_Curve = function(ddF, miss_ind, p){
  
  # create hit index based on the miss index.
  # This is to avoid multiple selections of a certain AE which corresponds to
  # multiple groups.
  ddF_temp = ddF %>%
    mutate(miss = miss_ind) %>%
    mutate(hit = ifelse(miss_ind == TRUE, FALSE, TRUE)) %>%
    select(AE_NAME, OR, hit, miss) %>%
    distinct()
  # found distinct positions 
  position = sapply(unique(ddF_temp$OR),
                    function(x) tail(which(ddF_temp$OR == x), n = 1),
                    simplify = T)
  # number of different positions 
  n_pos = length(position)
  
  # The sum of modified correlation metric(odds ratio here)
  # reduces to the standard Kolmogorov-Smirnov statistic when p = 0
  N_R = ddF_temp %>%
    filter(hit == TRUE) %>%
    mutate(OR_p = abs(OR)^p) %>%
    summarize(Nr = sum(OR_p)) %>%
    as.numeric()
  ## Number of miss hits
  N_miss = ddF_temp %>%
    summarize(N_M = sum(miss_ind)) %>%
    as.numeric()
  
  if (N_R == 0){
    hit_value = rep(0, n_pos)
  }else{
    if(p == 0){
      hit_value = cumsum(ddF_temp$hit / N_R)
      hit_value = hit_value[position]
    } else{
    OR_hit = ddF_temp %>%
      mutate(OR_p = abs(OR)^p * hit)
    hit_value = cumsum(OR_hit$OR_p / N_R)
    hit_value = hit_value[position]
    } }
  miss_value = cumsum(ddF_temp$miss / N_miss)
  miss_value = miss_value[position]
  return(list(hit = hit_value, miss = miss_value, pos = position))
}
