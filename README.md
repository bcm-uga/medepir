# medepir: MEthylation DEconvoluation PIpeline in R

Package {medepir} (MEthylation DEconvoluation PIpeline in R) provides practical implementation of recommended guidelines for inference of cell-type proportions with confounders, based on non-negative matrix factorization of DNA methylation data. 

## Installation

To get the current version from GitHub:

```R
remotes::install_github("bcm-uga/medepir")
```

To also build vignettes, use:

```R
# Take 1-5 min to build
remotes::install_github("bcm-uga/medepir", 
                        build_opts = c("--no-resave-data", "--no-manual"))
                        
# To see the two vignettes:
vignette("simulations", package = "medepir")
vignette("deconvolution", package = "medepir")
```
