
#---------------------------------------------------------------------------------------------------------------
##-KS_enrichment: perform Enrichment using proposed KS
###  Input: 
###       df (data.frame): two types
###
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
###       drug.case: Name of the target vaccine 
###       drug.control: Name of the Reference vaccine 
###       n_iter: Number of simulations to get null distribution of ES
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
###                                   AE_NAME              Ratio p_value
###                                 1 Abdomen scan normal  0.25    0.428
###                                 2 Abdominal discomfort 0.137   0.517
###                                 3 Abdominal distension 0.269   0.065
###                                 4 Abdominal mass       1       0.137
###                                 5 Abdominal pain       0.282   0
###                                 6 Abdominal pain lower 0.25    0.216
#---------------------------------------------------------------------------------------------------------------
#example: KS_enrichment(df = flu1, dd.group = group, drug.case = 'FLU')
KS_enrichment = function(df, dd.group, drug.case, drug.control=NULL, n_iter = 1000){
  # Get reporting rate estimate and related info for null simulation
  Ratio_result = get_ratio(df,drug.case,drug.control)
  Ratio = Ratio_result$Ratio
  n_DRUG = Ratio_result$n_DRUG
  data_info = Ratio_result$data_info
  # Get the group size for groups of interest
  dd.group=dd.group[!duplicated(dd.group),] %>% filter(!is.na(GROUP_NAME))
  df_size = dd.group %>%  group_by(GROUP_NAME) %>% summarise(GROUP_SIZE=n()) %>% ungroup()
  # Calculate the Enrichment Score for each group
  ks_true = get_ES(dd.group,Ratio[,1:2])  
  
  # Start simulations to get Null distribution of ES for each group of interest
  # Ratio_whole = Ratio
  # an initial data.frame with n_iter rows and n_Groups columns (each column corresponds to
  # the Null Distribution of a certain group)
  ks_null = data.frame(matrix(NA,ncol=nrow(ks_true), nrow = n_iter))
  Freq = rmultinom(n_iter, n_DRUG, data_info$n_AE/data_info$n)
  for ( i in 1:n_iter){
    # display the processing 
    if ( (i %% floor(n_iter/10)) == 0){
      cat(sprintf('Iteration Processing: %4.2f Percent \n', i/n_iter*100))
    }
    # simulate the data based on the null reporting rate and Poisson structure  
    data_null = data_info %>% mutate(Freq=Freq[,i])
    Ratio_null = data_null %>% mutate(Ratio = Freq/n_AE) %>% 
      dplyr::select(AE_NAME,Ratio)
    # get ES for a simulated data
    ks_null_temp = get_ES(dd.group, Ratio_null)
    ks_null[i,] = ks_null_temp$ES
    #Ratio_whole = merge(Ratio_whole, Ratio_null, by = 'AE_NAME', all.x = T, 
    #suffixes = c(as.character(i-1),as.character(i)))
  }
  # calculate p value for each AE
  # Ratio_whole[is.na(Ratio_whole)] = 0 # replace NA with 0
  # Ratio_mat = as.matrix(Ratio_whole[,3:ncol(Ratio_whole)])
  # p_AE = rowMeans((Ratio_mat-Ratio_whole$Ratio0)>=0)
  # Ratio_final = Ratio %>% mutate(p_value = p_AE)
  # p value for ks
  ks_true_df = data.frame(matrix(NA, nrow = 1, ncol = nrow(ks_true)))
  ks_true_df[1,] = ks_true$ES
  # combines true ES and permuted ES
  ks_all = rbind(ks_null, ks_true_df)
  p_value_ks = sapply(ks_all, function(x) mean(x[1:n_iter]>=x[n_iter+1]))
  
  Final_result_ks = data.frame(GROUP_NAME = ks_true$group, ES = ks_true$ES, p_value = p_value_ks)
  Final_result_ks = Final_result_ks %>% merge(df_size, by = 'GROUP_NAME') 
  
  return(list(Final_result=Final_result_ks, AE_info=Ratio))
}
