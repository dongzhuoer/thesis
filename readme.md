# undergraduate thesis online version
[![Build Status](https://travis-ci.com/dongzhuoer/thesis.svg?branch=bookdown)](https://travis-ci.com/dongzhuoer/thesis)



## usgae

- install WenQuanYi Zen Hei font (`fonts-wqy-zenhei` for Ubuntu)

- open R

```r
install.packages(c('magrittr', 'tidyverse', 'cowplot', 'bookdown'))

bookdown::render_book('thesis.Rmd')
browseURL('gitbook/index.html')
```


