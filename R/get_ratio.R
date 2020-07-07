
#--------------------------------------------------------------------------------------------------------
##------- get_ratio function: estimate reporting rate for each AE, which is the MLE of lambda in Poisson in formula (1）
### Input: 
###       data: two types of input (pairs or counts)
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
###       drug.case: Target vaccine type
###       drug.control: Reference vaccine type
###
### Output: 
###       a list with three elements: 
###           1. Ratio (data.frame): AE_NAME     Ratio
###                                  Insomnia     0.35
###                                   ...         ...
###                                  Chills       0.26
###
###           2. n_DRUG (single value): the number of target vaccine reports
###
###           3. data_info (data.frame): AE_NAME    n_AE
###                                      Insomnia    26
###                                       ...        ...
###                                      Chills      118
###
###           Note: 1 is the main result we will use ; 2 and 3 are used to generate the null
#--------------------------------------------------------------------------------------------------------
get_ratio = function(data,drug.case=drug.case, drug.control=NULL){
  data = as.data.frame(data)
  if (!is.numeric(data[,3])){
    #-- Data Type I case
    ## change the name of columns 
    names(data) = c('ID', 'DRUG_TYPE', 'AE_NAME')
    if(!is.null(drug.control)){
      drug.list=c(drug.case,drug.control)
      data=data[data$DRUG_TYPE %in% drug.list,]
    }
    data = data[complete.cases(data),] # in case there are NA's
    data = data %>% mutate(DRUG_TYPE = ifelse(DRUG_TYPE %in% drug.case, "DrugYes", "DrugNo")) 
    
    data = data %>% group_by(DRUG_TYPE, AE_NAME) %>% summarize(Freq = n())%>% ungroup()
    ## calculate useful info: n_DRUG, n_AE and n
    data = data %>%  group_by(DRUG_TYPE) %>% mutate(n_DRUG = round(sum(Freq))) %>% ungroup() %>%  
      group_by(AE_NAME) %>% mutate(n_AE = round(sum(Freq))) %>% ungroup() %>%
      mutate(n = round(sum(Freq))) %>% ungroup()
    ## calculate the estimate of reporting rate which I named Ratio
    data = data %>% mutate(Ratio = Freq/n_AE)
    ## save useful ratio into data_info (used later in the NULL simulation process)
    data_info = data %>% dplyr::select(AE_NAME,n_AE,n) %>% distinct() %>% arrange(AE_NAME)
    ## filter by target vaccine
    data = filter(data, DRUG_TYPE=="DrugYes") 
    ## get the reporting number of target vaccine
    n_DRUG = data$n_DRUG[1]
    
    data_Ratio = data %>% transmute(AE_NAME,Ratio,p_value=1-ppois(Freq-1,n_DRUG*n_AE/n))
    
  }else{
    #-- Data Type II case
    names(data) = c('DRUG_TYPE', 'AE_NAME', 'Freq')
    if(!is.null(drug.control)){
      drug.list=c(drug.case,drug.control)
      data=data[data$DRUG_TYPE %in% drug.list,]
    }
    data = data %>% mutate(DRUG_TYPE = ifelse(DRUG_TYPE %in% drug.case, "DrugYes", "DrugNo")) 
    ## calculate useful info: n_DRUG, n_AE and n
    data = data %>%  group_by(DRUG_TYPE) %>% mutate(n_DRUG = round(sum(Freq))) %>% ungroup() %>%  
      group_by(AE_NAME) %>% mutate(n_AE = round(sum(Freq))) %>% ungroup() %>%
      mutate(n = round(sum(Freq))) %>% ungroup()
    ## calculate the estimate of reporting rate which I named Ratio
    data = data %>% mutate(Ratio = Freq/n_AE)
    ## save useful ratio into data_info (used later in the NULL simulation process)
    data_info = data %>% dplyr::select(AE_NAME,n_AE,n) %>% distinct() %>% arrange(AE_NAME)
    ## filter by target vaccine
    data = filter(data, DRUG_TYPE=="DrugYes") 
    ## get the reporting number of target vaccine
    n_DRUG = data$n_DRUG[1]
    
    data_Ratio = data %>% transmute(AE_NAME,Ratio,p_value=1-ppois(Freq-1,n_DRUG*n_AE/n))
  }
  return(list(Ratio = data_Ratio,n_DRUG=n_DRUG,data_info=data_info))
}

