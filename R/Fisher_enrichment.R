
#---------------------------------------------------------------------------------------------------------------
##--Fisher_enrichment: performe Enrichment test using Fisher test
###  Input: 
###       df (data.frame): two types
###             Type I: ID   DRUG_TYPE   AE_NAME
###                     201    FLUN     Insomnia
###                     ...    ...       ...
###                     299    FLU      Chills
###
###             Type II: DRUG_TYPE   AE_NAME   COUNT(FREQ)
###                       FLUN      Insomnia   39
###                       ...        ...       ...
###                       FLU       Chills     13
###
###       dd.group (data.frame): dd.meddra (AE_NAME GROUP_NAME)
###
###       drug.case: Name of the Target vaccine 
###       drug.control: Name of the Reference vaccine 
###       n_iter: Number of simulations to get null distribution of ES
###       q.cut (numerical value): q value cut deciding the significance of each AE
###       or.cut (numerical value): odds ratio cut deciding the significance of each AE
###
###  Output:
###       a list with 2 data.frames:
###            1. Final_result (data.frame):
###                                             GROUP_NAME ES p_value GROUP_SIZE 
###         1                          Acid-base disorders  0       1          2    
###         2                          Allergic conditions  0       1         41    
###         3 Anaemias nonhaemolytic and marrow depression  0       1          1  
###         4                  Ancillary infectious topics  0       1          3  
###         5                     Angioedema and urticaria  0       1         17    
###         6               Anxiety disorders and symptoms  0       1         11    
### 
###            2. AE_info (data.frame):
###                                         AE_NAME       OR    p_value   
###                                    Abdominal mass    2.13   0.175     
###                                    Abdominal pain     0        1       
#---------------------------------------------------------------------------------------------------------------

Fisher_enrichment = function(df,dd.group,drug.case,drug.control=NULL,n_iter = 1000,q.cut=0.05,or.cut=1.5){
  # Get reporting rate estimate and related info for null simulation
  fisher_res = fisher_eachAE(df,drug.case,drug.control)
  fisher_res_withQ = fisher_res %>% mutate(qval = (qvalue(p_value))$qvalue)
  # Get the group size for groups of interest
  dd.group=dd.group[!duplicated(dd.group),] %>% filter(!is.na(GROUP_NAME))
  df_size = dd.group %>%  group_by(GROUP_NAME) %>% summarise(GROUP_SIZE=n()) %>% ungroup()
  # Calculate the Enrichment Score for each group
  fisher_result = fisher_test(dd.group,fisher_res_withQ,q.cut,or.cut,n_iter)
  
  Final_result_fisher = fisher_result %>% merge(df_size, by = 'GROUP_NAME') 
  return(list(Final_result=Final_result_fisher, AE_info=fisher_res[,1:3]))
}
