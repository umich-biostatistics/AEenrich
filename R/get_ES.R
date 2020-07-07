
#---------------------------------------------------------------------------------------------------------------
##------ get_ES function: calculate Enrichment Score over AE Groups
##       used in KS_enrichment function
###
###  Input: 
###       dd.group (data.frame): dd.meddra (AE_NAME GROUP_NAME)
###       data_Ratio (data.frame):  AE_NAME     Ratio
###                                 Insomnia     0.35
###                                  ...         ...
###                                 Chills       0.26
###
###  Output: 
###       a data.frame:  group            ES
###                      Respiratory      0.83
###                       ...             ...
###                      Infections       0.46
#----------------------------------------------------------------------------------------------------------------

get_ES<-function(dd.group,data_Ratio){
  
  
  # add zero to AE if it was not mentioned with the target vaccine
  ddF=merge(data_Ratio,dd.group,by = 'AE_NAME',all.y = TRUE) %>% mutate(Ratio = coalesce(Ratio, 0))
  ddF = ddF %>% arrange(desc(Ratio)) #order the observation by Ratio (reporting rate)
  
  # check if there are 0's
  flag_0 = any(ddF$Ratio == 0)
  # found distinct positions 
  position = sapply(unique(ddF$Ratio),function(x) tail(which(ddF$Ratio==x),n=1), simplify = T)
  # get interesting group name
  group.enrich = ddF %>% arrange(GROUP_NAME) %>% distinct(GROUP_NAME)
  group.enrich = group.enrich$GROUP_NAME
  ng=length(group.enrich)
  
  ES_vec = c()
  #dist.0_vec = c()
  for (j in 1:ng){
    # check which AE in this group
    hit_ind=ddF$GROUP_NAME==group.enrich[j] 
    AE_vec = (ddF$AE_NAME)[hit_ind]
    # get the miss index (handle one AE with multiple groups) 
    miss_ind = (hit_ind==F & !(ddF$AE_NAME %in% AE_vec))
    # get hit vector and miss vector
    value = ddF$Ratio
    h_m_lst = HitMiss_Curve(position,hit_ind = hit_ind, miss_ind = miss_ind)
    h_vec = h_m_lst$hit
    m_vec = h_m_lst$miss
    if (flag_0){
      n_pos = length(position)
      ES = max((h_vec-m_vec)[1:(n_pos-1)])*((h_vec-m_vec)[n_pos-1]>=0)
    }else{
      ES = max(h_vec-m_vec)
    }
    ES_vec = c(ES_vec,ES)
  }
  
  return(data.frame(group=group.enrich,ES = ES_vec))
}
