## Project Description

This repository contains the code necessary to reproduce my thesis. In a Bayesian setting, Gaussian Graphical Models parametrize the precision matrix of a multivariate normal distribution to capture the conditional dependence structure between variables, represented as a Markov Random Field. Simulations assessing the sensitivity of spike and slab priors to various non-Gaussian distributions were conducted, and differences between the accuracy of posterior probabilities at various thresholds were noted. A non-parametric transformation was applied to non-Gaussian data, and differences in predictive accuracy between models fit on Gaussian data and transformed non-Gaussian data were compared. This thesis was conducted under the supervision of Jack Jewson and Ruben Loaiza Maya and expanded upon the findings of the reasearch paper [Bayesian computation for high-dimensional Gaussian Graphical Models with spike-and-slab priors](https://arxiv.org/abs/2511.01875).

## Guide to Repository

R/ contains all the R scripts and R markdown scripts

Plots/ contains all the main figures

data/ contains empirical data used in this project

## Reproduce Results

1.  Download all files from this repository, click "<> Code" and then click "Download Zip"

2.  Open your R IDE (for example R Studio or Visual Studio Code)

3.  In your R IDE, set your working directory to the folder you downloaded by running `setwd()` with the path of the folder you download.

4.  Run the code chunk below to use all packages used within this repository

``` r
install.packages("renv")
renv::restore()
```
## Alternative to step 4: 

5. This code chunk below can be run if step 4 doesn't work

``` r
install.packages(c("metRology","MASS","mvtnorm","mombf","fredr","tidyr","tibble","tidyverse","","igraph","grid","png","gridExtra","dbplyr","sn","ggplot2","pROC","huge","forecast","purrr","ggnetwork","maps","ggraph","tseries"))
```

To render any R markdown scripts, run `rmarkdown::render()` in the console.

**Citation:** Sulem, D., Jewson, J., & Rossell, D. (2025). Bayesian computation for high-dimensional Gaussian Graphical Models with spike-and-slab priors. arXiv. https://doi.org/10.48550/ARXIV.2511.01875
