# rm(list=ls())
# options(stringsAsFactors = FALSE)
# library(dplyr)
# library(qvalue)
# 
# 
# ##-- import data
# ### data type 1
# # flu1 = read.csv('data/flu1.csv')
# # 
# # ### data type 2
# # 
# # flu2 = read.csv('data/flu2.csv')
# # 
# # 
# # ##-- import AE hierarchy
# # group = read.csv('data/GroupStructure.csv')
# 
# load('data/flu1.RData')
# load('data/flu2.RData')
# load('data/group.RData')
# 
# source('R/Enrichment_Ratio_V4.R')
# 
# ##---Experiment 1 
# drug.case = 'FLUN'
# drug.control = 'FLU'
# 
# # AEKS
# ## Data Type 1
# KS_result1 = ae.enrich(df=flu1, dd.group=group, 
#                                    drug.case=drug.case, drug.control=drug.control, method = 'aeks',
#                                    n_iter=1000)
# ## Data Type 2
# KS_result2 = ae.enrich(df=flu2, dd.group=group, 
#                                    drug.case=drug.case, drug.control=drug.control, method = 'aeks',
#                                    n_iter=1000)
# 
# 
# # AEFisher
# ## Data Type 1
# fisher_result1 = ae.enrich(df=flu1, dd.group=group, 
#                                        drug.case=drug.case, drug.control=drug.control, method = 'aefisher',
#                                        n_iter=1000, q.cut = 0.1, or.cut=1.5)
# 
# ## Data Type 2
# fisher_result2 = ae.enrich(df=flu2, dd.group=group, 
#                                   drug.case=drug.case, drug.control=drug.control, method = 'aefisher',
#                                   n_iter=1000, q.cut = 0.1, or.cut=1.5)
# 
# 
# 




