# Powerful and robust non-parametric association testing for microbiome data via a zero-inflated quantile approach (ZINQ) - Version 2

## Overview

The R package developed for Ling, W. et al. (2021). Powerful and robust non-parametric association testing for microbiome data via a zero-inflated quantile approach (ZINQ). Microbiome 9, 181.

## Instructions for use

From an `R` session, install `ZINQ` by:
```
devtools::install_github("wdl2459/ZINQ-v2")
```
And find the vignettes at: https://wdl2459.github.io/ZINQ-v2/ZINQ.Vignette.html.


To include the vignettes during installing the package:
```
devtools::install_github("wdl2459/ZINQ-v2", build_vignettes = TRUE, force=TRUE)
```
To view the vignettes, from the `R` terminal, type: 
```
browseVignettes("ZINQ")
```

From an `R` session, library the package by:
```
library(ZINQ)
```

v2.0 updates the logistic component to Firth logistic regression, in order to mitigate the bias when the presence-absence status of taxa is imbalanced across the different groups of the interested clinical variable(s).

Details can be found in the manual: https://github.com/wdl2459/ZINQ-v2/blob/dev/ZINQ_2.0.pdf.
