
#' Perform Adverse Event Enrichment Tests
#' 
#' The enrich function is used to perform Adverse event (AE) enrichment analysis.
#' Unlike the continuous gene expression data, AE data are counts. Therefore,
#' AE data has many zeros and ties. We propose two enrichment tests. AEFisher is
#' a modified Fisher's exact test based on pre-selected significant AEs, while
#' AEKS is based on a modified Kolmogorov-Smirnov statistic.
#' 
#' @param data a data.frame. Two data types are allowed. Type I data consisting
#' data on the report level, having ID, Drug type and AE name as the first 3
#' columns with covariates(optional) followed. Type II data have drug type and
#' AE name as the first two columns, with the 3rd and 4th Columns giving the
#' numbers of successes(have AE) and failures(Do not have AE) respectively, then 
#' followed by covariates. See example data for details. 
#' @param dd.group a data.frame with AE name and Group name. This data.frame have
#' the group information for each individual AE.
#' @param drug.case a character string for the target drug of interest.
#' @param drug.control a character string for the reference drug. If NULL(default),
#' all other drugs combined are the reference.
#' @param method a character string specifying the method for the enrichment test. 
#' It must take "aeks" (default) or "aefisher"; "aeks" is the rank-based 
#' enrichment test, and "aefisher" is the Fisher enrichment test. See details
#' described in the paper (see reference section of this document).
#' @param n_perms an integer value specifying the number of permutations in
#' permutation test.
#' @param covar a character vector specifying the columns of covariates, default
#' NULL.
#' @param p a numerical value to control the weight of the step, can take any
#' value between 0 and 1. If 0(default), reduces to the standard Kolmogorov-Smirnov
#' statistics.
#' @param q.cut a numerical value specifying the significance cut for q value 
#' of AEs in aefisher.
#' @param or.cut a numerical value specifying the significance cut for odds ratio 
#' of AEs in aefisher.
#' @param zero logical, default FALSE.If TRUE, add zero indicator to enrichment score.
#' @param min_size the minimum size of group required for enrichment analysis.
#' @param min_AE the minimum number of cases required to start counting
#' for a specific AE.
#' @param cores the number of cores to use for parallel execution.
#' @references Li, S. and Zhao, L. (2020). Adverse event enrichment tests using 
#' VAERS. \href{https://arxiv.org/abs/2007.02266}{arXiv:2007.02266}.
#' 
#' Subramanian, A.e.a. (2005). Gene set enrichment analysis: a knowledge-based
#' approach for interpreting genome-wide expression profiles. Proc Natl Acad
#' Sci U S A. Proceedings of the National Academy of Sciences. 102. 15545-15550. 
#' 
#' Tian, Lu & Greenberg, Steven & Kong, Sek Won & Altschuler, Josiah & Kohane, Isaac & Park,
#' Peter. (2005). Discovering statistically significant pathways in expression profiling studies.
#' Proceedings of the National Academy of Sciences of the United States of America.
#' 102. 13544-9. 10.1073/pnas.0506577102. 
#' 
#' @return A list containing 2 data.frames named **Final_result** and **AE_info**.
#' 
#' The **Final_result** data.frame contains the following columns: 
#' \itemize{
#'   \item{GROUP_NAME: }{AE group names}
#'   \item{ES: }{enrichment score}
#'   \item{p_value: }{p value of the enrichment test}
#'   \item{GROUP_SIZE: }{number of AEs per group} 
#' }
#' 
#' The **AE_info** contains the following columns:
#' \itemize{
#'   \item{AE_NAME: }{AE names}
#'   \item{OR: }{odds ratio for each individual AE}
#'   \item{p_value: }{p value for AE-drug association}
#'   \item{95Lower: }{lower bound of 95 percent confidence interval of odds ratio}
#'   \item{95Lower: }{upper bound of 95 percent confidence interval of odds ratio}
#'   \item{se(logOR): }{standard error of log odds ratio}
#' }
#' 
#' @examples
#' 
#' \donttest{
#'# AEKS
#'
#'### Type I data: data on report level
#'# enrich(data = covid1, covar = c("SEX", "AGE"), p = 0, method = "aeks",
#'#        n_perms = 1000, drug.case = "COVID19", dd.group = group, cores = 2,
#'#        drug.control = "OTHER", min_size = 5, min_AE = 10, zero = FALSE)
#'       
#'## Type II data: aggregated data
#'# enrich(data = covid2, covar = c("SEX", "AGE"), p = 0, method = "aeks",
#'#        n_perms = 1000, drug.case = "DrugYes", dd.group = group, cores = 2,
#'#        drug.control = "DrugNo", min_size = 5, min_AE = 10)
#'       
#'# AEFISHER
#'## Type I data: data on report level
#'# enrich(data = covid1, covar = c("SEX", "AGE"), p = 0, method = "aefisher",
#'#        n_perms = 1000, drug.case = "COVID19", dd.group = group,
#'#        drug.control = "OTHER", min_size = 5, min_AE = 10, q.cut = 0.05, 
#'#        or.cut = 1.5, cores = 2)
#'       
#'## Type II data: aggregated data
#'# enrich(data = covid2, covar = c("SEX", "AGE"), p = 0, method = "aefisher",
#'#        n_perms = 1000, drug.case = "DrugYes", dd.group = group,
#'#        drug.control = "DrugNo", min_size = 5, min_AE = 10, cores = 2)
#'       }



enrich = function(data, dd.group, drug.case, drug.control = NULL,
                  method = 'aeks', n_perms = 1000, covar = NULL, p = 0,
                  q.cut = 0.1, or.cut = 1.5, zero = FALSE, min_size = 5,
                  min_AE = 10, cores = detectCores()
                  ) {
  names(dd.group) = c('AE_NAME', 'GROUP_NAME')
  if (method == 'aeks'){
    KS_result = KS_enrichment(data, drug.case, drug.control, covar = covar,
                              dd.group = dd.group, n_perms = n_perms, p = p,
                              zero = zero, min_size = min_size, min_AE = min_AE,
                              cores = cores)
    return(KS_result)
  }else if (method == 'aefisher'){
    fisher_result = Fisher_enrichment(data, dd.group, drug.case, drug.control, 
                                      n_perms = n_perms, q.cut = q.cut,
                                      or.cut = or.cut, zero = zero,
                                      min_size = min_size, covar = covar,
                                      min_AE = min_AE, cores = cores)
    return(fisher_result)
  }else{
    stop('Please choose one of two methods: aeks or fisher')
  }
}

#' @description Perform Adverse Event Enrichment Tests
#' The enrich function is used to perform Adverse event (AE) enrichment analysis.
#' Unlike the continuous gene expression data, AE data are counts. Therefore,
#' AE data has many zeros and ties. We propose two enrichment tests. AEFisher is
#' a modified Fisher's exact test based on pre-selected significant AEs, while
#' AEKS is based on a modified Kolmogorov-Smirnov statistic.
#' 
#' Use the function `enrich` to fit models and inspect results.
#' 
#' See our \href{https://github.com/umich-biostatistics/AEenrich}{Github home page} 
#' or run ?enrich for examples.
"_PACKAGE"
