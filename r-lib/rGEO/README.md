# rGEO


## Overview

This package contains functions to read GSE raw data and returns a tidy format. It's especially useful for chip annotation file, the format of which is really in a messy. You just need to provide the raw data file and it will give you a nice table containing probes and [HUGO](https://www.genenames.org/) gene symbols.


## Installation

```r
if (!('devtools' %in% .packages(T))) install.packages('devtools');
remotes::install_github('dongzhuoer/thesis/r-lib/rGEO');
```



## Usage

Inpatient people please refer to `vignette('rGEO')`.

For why this package, read this `vignette('probe2symbol')`.



# Develop

> Wired as it seems, it's best way I can think about so far to maintain reproducibility. 

Basically, you make change in `R/platform.R`,  `Load All` (C-S-L, or ), then Knit `R-raw/platform.Rmd` to check whether the rule returns exactly what you want. (you can also execute all chunk before, and then repeat modify rule - load- run a specific chunk to save time)


