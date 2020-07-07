
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
##------- fisher_test: permute significant AE label; for each permutated data, construct a contingency table for each Group and compute OR*Zero indicator  
#                                                                         grp
#                                                                sig    0    1
#                                                                  0 1614   17
#                                                                  1  197    0

###       Used in Fisher_enrichment function
###  Input: 
###       dd.group (data.frame): dd.meddra (AE_NAME GROUP_NAME)
###
###       fisher_res (data.frame): 
###                                  AE_NAME              OR   qval    isRatio0
###
###                                1 Abdomen scan normal  0.25  0.611   FALSE
###                                2 Abdominal discomfort 0.137 0.684   FALSE
###                                3 Abdominal distension 0.269 0.273   FALSE
###
###       q.cut (numerical value): q value cut deciding the significance of each AE
###       or.cut (numerical value): odds ratio cut deciding the significance of each AE
###       n_iter: number of permutations 
###
###  Output: 
###            1. Final_result (data.frame):
###                                             GROUP_NAME ES p_value 
###         1                          Acid-base disorders  0       1 
###         2                          Allergic conditions  0       1   
###         3 Anaemias nonhaemolytic and marrow depression  0       1       
###         4                  Ancillary infectious topics  0       1         
###         5                     Angioedema and urticaria  0       1          
###         6               Anxiety disorders and symptoms  0       1         
#---------------------------------------------------------------------------------------------------------------
fisher_test<-function(dd.group, fisher_res, q.cut, or.cut, n_iter){
  
  
  # add statistics for AE not mentioned with the target vaccine. (OR=0,q=1, isRatio0=TRUE)
  ddF=merge(fisher_res,dd.group,by = 'AE_NAME',all.y = TRUE) %>% mutate(OR = coalesce(OR, 0),
                                                                        qval = coalesce(qval, 1),
                                                                        isRatio0 = coalesce(isRatio0,TRUE))
  # determine the significance of each AE ??? 1.5 should be replaced by OR.cut 
  ddF$sig=factor(ifelse((ddF$qval<q.cut)&(ddF$OR>or.cut),1,0), levels = c(0,1))
  
  # get interesting group name
  group.enrich = ddF %>% arrange(GROUP_NAME) %>% distinct(GROUP_NAME)
  group.enrich = group.enrich$GROUP_NAME
  ng=length(group.enrich)
  
  # put group name together for a single AE.
  # For example, if AE "pain" belongs to two groups "G1" and "G2", the original data.frame 
  # looks like 
  #           AE_NAME    OR   isRatio0  sig     GROUP_NAME
  #           pain       1.3   FALSE     0          G1
  #           pain       1.3   FALSE     0          G2
  # But after we summarise, it will be like
  #           AE_NAME    OR   isRatio0  sig     GROUP_NAME
  #           pain       1.3   FALSE     0      list(G1,G2)
  ddF_new = ddF %>% group_by(AE_NAME) %>%  summarise(OR = OR[1],
                                                     isRatio0 = isRatio0[1],
                                                     sig = sig[1],
                                                     GROUP_NAME = list(GROUP_NAME))
  # Calculate true ES
  # initialize values 
  ES=c()
  grp_info = list()
  for(j in 1:ng){
    # This group indicator function can handle one AE with multiple groups
    grp=sapply(ddF_new$GROUP_NAME, function(x) group.enrich[j] %in% x)
    # store the grp info so that can speed up the program
    grp_info[[j]] = grp
    # calculate 0 proportion for each group
    in_0 = sum(ddF_new$isRatio0[grp==T])/sum(grp)
    out_0 = sum(ddF_new$isRatio0[grp==F])/sum(grp==F)
    # do one-sided fisher's exact test 
    table2=table(sig=ddF_new$sig,grp)
    or = table2[1,1]*table2[2,2]/(table2[1,2]*table2[2,1])
    # in_0 is the ratio of zero (#zero AE/  total #AE in target group) for target drug; 
    # out_0 is the ratio of zero (#zero AE/  total #AE in other groups) for target drug; 
    ES[j]=or*(in_0<=out_0)
    if (is.na(ES[j])){
      ES[j] = 0 # Inf*0 = NAN
    }
  }
  # get the index of ES not 0
  index = which(ES!=0)
  
  # Calculate the null distribution of ES
  ES_null_df = data.frame(matrix(NA,ncol=length(ES), nrow = n_iter))
  for (perm in 1:n_iter){
    # display the processing 
    if ( (perm %% floor(n_iter/10)) == 0){
      cat(sprintf('Iteration Processing: %4.2f Percent \n', perm/n_iter*100))
    }
    # use permutation to construct null data.frame
    ind = sample(nrow(ddF_new))
    ddF_null = ddF_new %>% mutate(sig = sig[ind], isRatio0=isRatio0[ind])
    ES_null= rep(0, ng)
    for(j in index){
      grp=grp_info[[j]]
      # calculate 0 proportion 
      in_0 = sum(ddF_null$isRatio0[grp==T])/sum(grp)
      out_0 = sum(ddF_null$isRatio0[grp==F])/sum(grp==F)
      # do one-sided fisher's exact test
      table2=table(sig=ddF_null$sig,grp)
      or = table2[1,1]*table2[2,2]/(table2[1,2]*table2[2,1])
      # in_0 is the ratio of zero (#zero AE/  total #AE in target group) for target drug; 
      # out_0 is the ratio of zero (#zero AE/  total #AE in other groups) for target drug; 
      ES_null[j]=or*(in_0<=out_0)
      if (is.na(ES_null[j])){
        ES_null[j] = 0 # Inf*0 = NAN
      }
    }
    ES_null_df[perm,] = ES_null
  }
  # calculate tail probability (p value)
  ES_true_df = data.frame(matrix(NA, nrow = 1, ncol = length(ES)))
  ES_true_df[1,] = ES
  # combines true ES and permuted ES
  ES_all = rbind(ES_null_df, ES_true_df)
  p_value = sapply(ES_all, function(x) mean(x[1:n_iter]>=x[n_iter+1]))
  
  
  res=cbind.data.frame(GROUP_NAME=group.enrich,ES=ES,p_value=p_value)
  
  return(res)
} 
