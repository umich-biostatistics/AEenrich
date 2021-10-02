
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R package `AEenrich`

[![](https://img.shields.io/badge/devel%20version-1.1.0-blue.svg)](https://github.com/umich-biostatistics/AEenrich)
[![](https://img.shields.io/github/languages/code-size/umich-biostatistics/AEenrich.svg)](https://github.com/umich-biostatistics/AEenrich)
[![](http://cranlogs.r-pkg.org/badges/grand-total/AEenrich?color=blue)](https://cran.r-project.org/package=AEenrich)
[![](http://cranlogs.r-pkg.org/badges/last-month/AEenrich?color=green)](https://cran.r-project.org/package=AEenrich)
[![](https://img.shields.io/badge/Read%20on-arXiv-orange.svg)](https://arxiv.org/abs/2007.02266)
[![CRAN
checks](https://cranchecks.info/badges/summary/AEenrich)](https://cran.r-project.org/web/checks/check_results_AEenrich.html)

## Overview

We extend existing gene enrichment tests to perform adverse event
enrichment analysis. Unlike the continuous gene expression data, adverse
event data are counts. Therefore, adverse event data has many zeros and
ties. We propose two enrichment tests. One is a modified Fisher’s exact
test based on pre-selected significant adverse events, while the other
is based on a modified Kolmogorov-Smirnov statistic. We add covariate
adjustment to improve the analysis.

## Install from CRAN

``` r
install.packages("AEenrich")
```

Then, load the package with

``` r
library(AEenrich)
```

## Install from Github

If the devtools package is not yet installed, install it first:

``` r
install.packages('devtools')
```

Then run:

``` r
# install AEenrich from Github:
devtools::install_github('umich-biostatistics/AEenrich') 
library(AEenrich)
```

## Example usage

For documentation pages:

``` r
?AEenrich
?enrich
?count_cases
```

### Quick example:

``` r
# AEKS
## Type I data: data on report level
enrich(data = covid1, covar = c("SEX", "AGE"), p = 0, method = "aeks",
       n_perms = 1000, drug.case = "COVID19", dd.group = group_info,
       drug.control = "OTHER", min_size = 5, min_AE = 10, zero = FALSE)
## Type II data: aggregated data
enrich(data = covid2, covar = c("SEX", "AGE"), p = 0, method = "aeks",
       n_perms = 1000, drug.case = "DrugYes", dd.group = group_info,
       drug.control = "DrugNo", min_size = 5, min_AE = 10)
# AEFISHER
## Type I data: data on report level
enrich(data = covid1, covar = c("SEX", "AGE"), p = 0, method = "aefisher",
       n_perms = 1000, drug.case = "COVID19", dd.group = group_info,
       drug.control = "OTHER", min_size = 5, min_AE = 10, q.cut = 0.05, 
       or.cut = 1.5, cores = 8)
## Type II data: aggregated data
enrich(data = covid2, covar = c("SEX", "AGE"), p = 0, method = "aefisher",
       n_perms = 1000, drug.case = "DrugYes", dd.group = group_info,
       drug.control = "DrugNo", min_size = 5, min_AE = 10)

## Convert type I data to type II data
count_cases (data = covid1, drug.case = "COVID19", drug.control = "OTHER",
             covar_cont = c("AGE"), covar_disc = c("SEX"),
             breaks = list(c(16,30,50,65,120)))
```

### Current Suggested Citation

1.  Li, S. and Zhao, L. (2020). Adverse event enrichment tests using
    VAERS. [arXiv:2007.02266](https://arxiv.org/abs/2007.02266).

2.  Subramanian, A.e.a. (2005). Gene set enrichment analysis: a
    knowledge-based approach for interpreting genome-wide expression
    profiles. Proc Natl Acad Sci U S A. Proceedings of the National
    Academy of Sciences. 102. 15545-15550.

3.  Tian, Lu & Greenberg, Steven & Kong, Sek Won & Altschuler, Josiah &
    Kohane, Isaac & Park, Peter. (2005). Discovering statistically
    significant pathways in expression profiling studies. Proceedings of
    the National Academy of Sciences of the United States of America.
    102. 13544-9. 10.1073/pnas.0506577102.
