
#' Perform Adverse Event Enrichment Tests
#' 
#' Adverse event (AE) enrichment analysis. Unlike the continuous gene expression data, AE 
#' data are counts. Therefore, AE data has many zeros and ties. We propose two 
#' enrichment tests. AEFisher is a modified Fisher's exact test based on pre-selected 
#' significant AEs, while AEKS is based on a modified Kolmogorov-Smirnov 
#' statistic.
#' 
#' @param df a data.frame with 3 columns. The function allows two data types. One type (data type I) consists data on the report level, including
#' ID, Drug type and AE name. The other type (data type II) consists of aggregated data, including
#' drug type, AE name and Count. Data should be ordered as 
#' ID, Drug type, AE name in data type I and Drug type, AE name, Count in data type II.
#' @param dd.group a data.frame with AE name and Group 
#' name. This data.frame have the group information for each individual AE.
#' @param drug.case a character string for the target drug of interest.
#' @param drug.control a character string for the reference drug. If NULL(default),
#' all other drugs combined are the reference.
#' @param method a character string specifying the method for the enrichment test. 
#' It must take "aeks" (default) or "aefisher"; "aeks" is the rank-based 
#' enrichment test, and "aefisher" is the modified Fisher enrichment test. See
#' details described in the paper (see reference section of this document). 
#' @param n_iter an integer value specifying the number of iterations in aeks 
#' method or the number of permutations in aefisher.
#' @param q.cut a numerical value specifying the significance cut for q value 
#' of AEs in aefisher.
#' @param or.cut a numerical value specifying the significance cut for odds ratio 
#' of AEs in aefisher.
#' @param seed a numeric seed for reproducible analysis.
#' @param verbose logical, if TRUE, print iterations. If FALSE, silence printing
#' to the console. Default is verbose = FALSE.
#' 
#' @references Li, S. and Zhao, L. (2020). Adverse event enrichment tests using 
#' VAERS. \href{https://arxiv.org/abs/2007.02266}{arXiv:2007.02266}.
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
#' The **AE_info** in 'aeks' contains the following columns:
#' \itemize{
#'   \item{AE_NAME: }{AE names}
#'   \item{Ratio: }{reporting rate}
#'   \item{p_value: }{p value for AE-drug association}
#' }
#' 
#' The **AE_info** in 'aefisher' contains the following columns: 
#' \itemize{
#'   \item{AE_NAME: }{AE names}
#'   \item{OR: }{Odds ratio for AE-drug association}
#'   \item{p_value: }{p value for AE-drug association}
#' }
#' 
#' @examples 
#' 
#' drug.case = 'FLUN'
#' drug.control = 'FLU'
#' 
#' # AEKS
#' ## Input data using data Type I
#' KS_result1 = enrich(df = flu1, dd.group = group, drug.case = drug.case, 
#'                     drug.control = drug.control, method = 'aeks', n_iter = 10)
#' ## Input data using data Type II
#' \donttest{
#' KS_result2 = enrich(df = flu2, dd.group = group, drug.case = drug.case, 
#'                     drug.control = drug.control, method = 'aeks', n_iter = 1000)
#' 
#' # AEFisher
#' fisher_result1 = enrich(df = flu1, dd.group = group, drug.case = drug.case, 
#'                         drug.control = drug.control, method = 'aefisher', 
#'                         n_iter = 1000, q.cut = 0.1, or.cut=1.5)
#' }



enrich = function(df, dd.group, drug.case, drug.control = NULL, method = 'aeks',  
                    n_iter = 1000, q.cut = 0.1, or.cut = 1.5, seed = NULL, 
                    verbose = FALSE) {
  if(!is.null(seed)) { set.seed(seed) }
  names(dd.group) = c('AE_NAME', 'GROUP_NAME')
  if (method == 'aeks'){
    KS_result = KS_enrichment(df, dd.group, drug.case, drug.control, n_iter = n_iter, 
                              verbose = verbose)
    return(KS_result)
  }else if (method == 'aefisher'){
    fisher_result = Fisher_enrichment(df, dd.group, drug.case, drug.control, 
                                      n_iter = n_iter, q.cut = q.cut, or.cut = or.cut, 
                                      verbose = verbose)
    return(fisher_result)
  }else{
    stop('Please choose one of two methods: aeks or fisher')
  }
}

#' @description Adverse event (AE) enrichment 
#' analysis. Unlike the continuous gene expression data, AE data are counts. 
#' Therefore, AE data has many zeros and ties. We propose two enrichment tests. 
#' AEFisher is a modified Fisher's exact test based on pre-selected significant AEs, 
#' while AEKS is based on a modified Kolmogorov-Smirnov statistic.
#' 
#' Use the function `enrich` to fit models and inspect results.
#' 
#' See our \href{https://github.com/umich-biostatistics/AEenrich}{Github home page} 
#' or run ?enrich for examples.
"_PACKAGE"
