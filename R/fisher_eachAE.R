
#--------------------------------------------------------------------------------------------------------
##------- fisher_eachAE function: fisher's exact test for each AE (called by Fisher_enrichment)
###
### Input: 
###       data: two types
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
###       a data.frame:
###                   AE_NAME       OR    p_value   isRatio0
###               Abdominal mass    2.13   0.175      FALSE
###               Abdominal pain     0        1       TRUE   
#--------------------------------------------------------------------------------------------------------

fisher_eachAE = function(data,drug.case=drug.case, drug.control=NULL){
  data = as.data.frame(data)
  if (!is.numeric(data[,3])){
    #-- Input Data Type I 
    ## change the name of columns 
    names(data) = c('ID', 'DRUG_TYPE', 'AE_NAME')
    if(!is.null(drug.control)){
      drug.list=c(drug.case,drug.control)
      data=data[data$DRUG_TYPE %in% drug.list,]
    }
    data = data[complete.cases(data),] # in case there are NA's
    data = data %>% mutate(DRUG_TYPE = ifelse(DRUG_TYPE %in% drug.case, "DrugYes", "Z_DrugNo")) 
    
    data = data %>% group_by(DRUG_TYPE, AE_NAME) %>% summarize(Freq = n())%>% ungroup()
    ## calculate useful info: n_DRUG, n_AE and n
    data = data %>%  group_by(DRUG_TYPE) %>% mutate(n_DRUG = round(sum(Freq))) %>% ungroup() %>%  
      group_by(AE_NAME) %>% mutate(n_AE = round(sum(Freq))) %>% ungroup() %>%
      mutate(n = round(sum(Freq))) %>% ungroup()
    ## extract info of each AE
    data$DRUG_TYPE=factor(data$DRUG_TYPE)
    dataList=split(data,data$AE_NAME)
    pval=OR=isRatio0=c()
    ## construct a loop to do fisher's test for each AE
    for(l in 1:length(dataList)){
      s=dataList[[l]]
      ## calculate the column of no_AE
      if (nrow(s)==1){
        s_c = data.frame(DRUG_TYPE=levels(s$DRUG_TYPE)[which(levels(s$DRUG_TYPE)!=s$DRUG_TYPE)],
                         AE_NAME=s$AE_NAME,Freq=0,n_DRUG = s$n-s$n_DRUG, n_AE=s$n_AE,
                         n=s$n)
        s = rbind(s,s_c)
      }
      s.noAE=cbind.data.frame(DRUG_TYPE=s[,1],AE_NAME=rep("Z_AE",nrow(s)),Freq=(s$n_DRUG-s$Freq))
      ## extract the columns of AE
      s.AE=s[,c("DRUG_TYPE","AE_NAME","Freq")]
      ## combine AE column and no_AE column together
      ss=rbind.data.frame(s.AE,s.noAE)
      dd.array = xtabs(Freq ~ DRUG_TYPE + AE_NAME, data=ss)
      ## do fisher test
      isRatio0[l] = !dd.array[1,1]
      foo=fisher.test(dd.array)
      pval[l]=foo$p.value
      OR[l]=foo$estimate
    }
    pval = ifelse(pval>1,1,pval)
    res=cbind.data.frame(AE_NAME=names(dataList),OR=round(OR,2),p_value=pval,isRatio0=isRatio0)
    
  }else{
    #-- Input Data Type II 
    names(data) = c('DRUG_TYPE', 'AE_NAME', 'Freq')
    if(!is.null(drug.control)){
      drug.list=c(drug.case,drug.control)
      data=data[data$DRUG_TYPE %in% drug.list,]
    }
    data = data %>% mutate(DRUG_TYPE = ifelse(DRUG_TYPE %in% drug.case, "DrugYes", "Z_DrugNo")) 
    ## calculate useful info: n_DRUG, n_AE and n
    data = data %>%  group_by(DRUG_TYPE) %>% mutate(n_DRUG = round(sum(Freq))) %>% ungroup() %>%  
      group_by(AE_NAME) %>% mutate(n_AE = round(sum(Freq))) %>% ungroup() %>%
      mutate(n = round(sum(Freq))) %>% ungroup()
    ## extract info of each AE
    data$DRUG_TYPE=factor(data$DRUG_TYPE)
    dataList=split(data,data$AE_NAME)
    pval=OR=isRatio0=c()
    ## construct a loop to do fisher's test for each AE
    for(l in 1:length(dataList)){
      s=dataList[[l]]
      ## calculate the column of no_AE
      if (nrow(s)==1){
        s_c = data.frame(DRUG_TYPE=levels(s$DRUG_TYPE)[which(levels(s$DRUG_TYPE)!=s$DRUG_TYPE)],
                         AE_NAME=s$AE_NAME,Freq=0,n_DRUG = s$n-s$n_DRUG, n_AE=s$n_AE,
                         n=s$n)
        s = rbind(s,s_c)
      }
      s.noAE=cbind.data.frame(DRUG_TYPE=s[,1],AE_NAME=rep("Z_AE",nrow(s)),Freq=(s$n_DRUG-s$Freq))
      ## extrac the column of AE
      s.AE=s[,c("DRUG_TYPE","AE_NAME","Freq")]
      ## combine AE column and no_AE column together
      ss=rbind.data.frame(s.AE,s.noAE)
      dd.array = xtabs(Freq ~ DRUG_TYPE + AE_NAME, data=ss)
      ## do fisher test
      isRatio0[l] = !dd.array[1,1]
      foo=fisher.test(dd.array)
      pval[l]=foo$p.value
      OR[l]=foo$estimate
    }
    pval = ifelse(pval>1,1,pval)
    res=cbind.data.frame(AE_NAME=names(dataList),OR=round(OR,2),p_value=pval,isRatio0=isRatio0)
  }
  return(res)
}
