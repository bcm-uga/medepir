# medepir
Data challenge methods 

The `medepir` package (MEthylation DEconvoluation PIpeline in R) provides practical implementation of recommended guidelines for inference of cell-type proportions with confounders, based on non-negative matrix factorization of DNA methylation data. 

## Installation

To get the current development version from github:

```
git clone https://github.com/bcm-uga/medepir.git
cd medepir 
R
```

## Build package

```R
install.packages("devtools")
# bigstatr ggplot2 RefFreeEWAS doParallel
devtools::load_all(); devtools::document(); devtools::install()
```

## Build vignettes

```R
setwd("vignette")
rmarkdown::render("vignette_deconv.Rmd")
rmarkdown::render("vignette_simu.Rmd")
```

