
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ProtE

<!-- badges: start -->
<!-- badges: end -->

The goal of ProtE is to get proteomics data from softwares like DIA-NN,
and Proteome Discoverer and apply a variety of proccesses like
Filtering, Imputation and Normalization and afterwards it checks the
quality of the data. Then it performs a statistical analysis between the
different Groups of samples in your data (specifically Mann-Whitney,
t-test, Kruskal-Wallis, ANOVA),and proceeds to create a heatmap, PCA
plots, Boxplots and Violin plots and a protein abundance Rank plot.

## Installation

You can install the development version of PACKAGE from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("theomargel/PACKAGE")
```
