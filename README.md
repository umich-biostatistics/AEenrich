
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R package `AEenrich`

[![](https://img.shields.io/badge/devel%20version-1.0.0-blue.svg)](https://github.com/umich-biostatistics/AEenrich)
[![](https://img.shields.io/github/languages/code-size/umich-biostatistics/AEenrich.svg)](https://github.com/umich-biostatistics/AEenrich)
[![](http://cranlogs.r-pkg.org/badges/grand-total/AEenrich?color=blue)](https://cran.r-project.org/package=AEenrich)
[![](http://cranlogs.r-pkg.org/badges/last-month/AEenrich?color=green)](https://cran.r-project.org/package=AEenrich)
[![](https://img.shields.io/badge/Read%20on-arXiv-orange.svg)](https://arxiv.org/abs/2007.02266)
[![CRAN
checks](https://cranchecks.info/badges/summary/AEenrich)](https://cran.r-project.org/web/checks/check_results_AEenrich.html)

## Overview

Adverse event (AE) enrichment analysis. Unlike continuous gene
expression data, AE data are counts. Therefore, AE data has many zeros
and ties. We propose two enrichment tests, AEFisher is a modified
Fisher’s exact test based on pre-selected significant AEs, while AEKS
is based on a modified Kolmogorov-Smirnov statistic.

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
```

### Quick example:

``` r
##---Experiment 1 
drug.case = 'FLUN'
drug.control = 'FLU'

# AEKS
## Data Type 1
KS_result1 = enrich(df = flu1, dd.group = group, drug.case = drug.case, 
                    drug.control = drug.control, method = 'aeks', n_iter = 1000)
## Data Type 2
KS_result2 = enrich(df = flu2, dd.group = group, drug.case = drug.case, 
                    drug.control = drug.control, method = 'aeks', n_iter = 1000)

# AEFisher
## Data Type 1
fisher_result1 = enrich(df = flu1, dd.group = group, drug.case = drug.case, 
                        drug.control = drug.control, method = 'aefisher', 
                        n_iter = 1000, q.cut = 0.1, or.cut=1.5)
## Data Type 2
fisher_result2 = enrich(df = flu2, dd.group = group, drug.case = drug.case, 
                        drug.control = drug.control, method = 'aefisher', 
                        n_iter = 1000, q.cut = 0.1, or.cut=1.5)
```

### Current Suggested Citation

Li, S. and Zhao, L. (2020). Adverse event enrichment tests using VAERS.
[arXiv:2007.02266](https://arxiv.org/abs/2007.02266).
