
#' Perform Adverse Event Enrichment Tests
#' 
#' Gene enrichment tests to perform adverse 
#' event (AE) enrichment analysis. Unlike the continuous gene expression data, AE 
#' data are counts. Therefore, AE data has many zeros and ties. We propose two 
#' enrichment tests. One is a modified Fisher's exact test based on pre-selected 
#' significant AEs, while the other is based on a modified Kolmogorov-Smirnov 
#' statistic.
#' 
#' @param df a data.frame with 3 columns. Two types are acceptable. One is with 
#' three columns corresponding to ID, Drug type and AE name. The other one is 
#' three columns with Vaccine type, AE name and Count. The order matters. Either 
#' order the columns ID, Drug type, AE name or Drug type, AE name, Count.
#' @param dd.group a data.frame with first column AE name and second column Group 
#' name. This data.frame provides the group for which each AE belongs.
#' @param drug.case a character string of target vaccine type.
#' @param drug.control a character string of reference vaccine type. If NULL(default),
#' all other types shown in the data are reference.
#' @param method a character string specifying the method for the enrichment test. 
#' It must be "aeks" (default) or "aefisher"; "aeks" means the rank-based 
#' enrichment test, and "aefisher" means the Fisher enrichment test. These two 
#' methods are described in the paper (see reference section of this document). 
#' @param n_iter an integer value specifying the iteration time for the aeks 
#' method or the permutation time for the aefisher method.
#' @param q.cut a numerical value specifying the significance cut for q value 
#' of AEs when using aefisher enrichment (method="aefisher").
#' @param or.cut a numerical value specifying the significance cut for odds ratio 
#' of AEs when using aefisher enrichment (method="aefisher").
#' @param seed set a numeric seed for reproducible analysis.
#' 
#' @references Li, S. and Zhao, L. (2020). Adverse event enrichment tests using 
#' VAERS. \href{https://arxiv.org/abs/2007.02266}{arXiv:2007.02266}.
#' 
#' @return A list with 
#' \itemize{
#'   \item{Final_result: }{data.frame with GROUP_NAME, ES, p_value, and GROUP_SIZE}
#'   \item{AE_info: }{data.frame with AE_name, Ratio, and p_value}
#' }
#' 
#' @examples 
#' 
#' ##---Experiment 1 
#' drug.case = 'FLUN'
#' drug.control = 'FLU'
#' 
#' # AEKS
#' ## Data Type 1
#' KS_result1 = enrich(df = flu1, dd.group = group, drug.case = drug.case, 
#'                     drug.control = drug.control, method = 'aeks', n_iter = 1000)
#' ## Data Type 2
#' KS_result2 = enrich(df = flu2, dd.group = group, drug.case = drug.case, 
#'                     drug.control = drug.control, method = 'aeks', n_iter = 1000)
#' 
#' # AEFisher
#' ## Data Type 1
#' fisher_result1 = enrich(df = flu1, dd.group = group, drug.case = drug.case, 
#'                         drug.control = drug.control, method = 'aefisher', 
#'                         n_iter = 1000, q.cut = 0.1, or.cut=1.5)
#' ## Data Type 2
#' fisher_result2 = enrich(df = flu2, dd.group = group, drug.case = drug.case, 
#'                         drug.control = drug.control, method = 'aefisher', 
#'                         n_iter = 1000, q.cut = 0.1, or.cut=1.5)

enrich = function(df, dd.group, drug.case, drug.control = NULL, method = 'aeks',  
                    n_iter = 1000, q.cut = 0.1, or.cut = 1.5, seed = NULL) {
  if(!is.null(seed)) { set.seed(seed) }
  else { set.seed(floor(runif(1, 1, 1000))) }
  names(dd.group) = c('AE_NAME', 'GROUP_NAME')
  if (method == 'aeks'){
    KS_result = KS_enrichment(df, dd.group, drug.case, drug.control, n_iter = n_iter)
    return(KS_result)
  }else if (method == 'aefisher'){
    fisher_result = Fisher_enrichment(df, dd.group, drug.case, drug.control, 
                                      n_iter = n_iter, q.cut = q.cut, or.cut = or.cut)
    return(fisher_result)
  }else{
    stop('Please choose one of two methods: aeks or fisher')
  }
}

#' @description Gene enrichment tests to perform adverse event (AE) enrichment 
#' analysis. Unlike the continuous gene expression data, AE data are counts. 
#' Therefore, AE data has many zeros and ties. We propose two enrichment tests. 
#' One is a modified Fisher's exact test based on pre-selected significant AEs, 
#' while the other is based on a modified Kolmogorov-Smirnov statistic.
#' 
#' Use the function `enrich` to fit models and inspect results.
#' 
#' See our \href{https://github.com/umich-biostatistics/AEenrich}{Github home page} 
#' or run ?enrich for examples.
"_PACKAGE"
