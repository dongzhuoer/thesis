# Overview
[![Build Status](https://travis-ci.com/dongzhuoer/thesis.svg?branch=master)](https://travis-ci.com/dongzhuoer/thesis)

Graduation Project at Nankai University

There is a brief [introduction](https://bookdown.dongzhuoer.com/dongzhuoer/thesis/abstract.html) about the scientific significance, but I think the most valuable part is the R packages:

1. **[rGEO](https://dongzhuoer.github.io/rGEO)** is the core package, it converts GEO microarray data to a standard format, i.e, map probes to HUGO gene symbol.
1. **[qGSEA](https://dongzhuoer.github.io/qGSEA)** is developed based on rGEO, it can quickly prepare GSEA input files from GEO. The biggest feature is a GUI in browser, user just need to specify an accession in the simplest case.
1. **[hgnc](https://dongzhuoer.github.io/hgnc)** is a dependency of rGEO. In short, rGEO find where information lies and hgnc convert that information to standard format. hgnc may be also useful in other situations so I make it a separated package.
1. **[rGEO.data](https://dongzhuoer.github.io/rGEO)** is mainly for prevening rGEO from getting too large, it contains several big data used in rGEO testing  and functions to create them. 



# Usage

## figures

You can produce the [figures](figure.md) used in paper and presentation by:

```r
install.packages(c("rmarkdown", "prettydoc"))
rmarkdown::render("figure.Rmd", output_format = "prettydoc::html_pretty")
browseURL("figure.html")
```

## full text

You can read full text [online](https://bookdown.dongzhuoer.com/dongzhuoer/thesis/), or build it on your own:

```r
setwd("bookdown")
install.packages(c("bookdown"))
# install "WenQuanYi Zen Hei" font for proper rendering of Chinese character in plot
bookdown::render_book('')
browseURL("gitbook/index.html")
```

## others

You can explore other code,

- [整个研究的代码](notebook/)
- [在工作站用到的代码](workstation/)
- [总结出的文件格式，数据库架构等](format.md)

under the following conventions,

- `lem4` is alias for ssh workstation
- output is rendered on 2019-04-03, and masked by `#> `
- all code containing output can run, except for `gsea_output_full.rds` (too big)



## progress outline

Understand the principle of GSEA. Then its input file format, parameter meaning, and how to interpret its output.

Explore the architecture of GEO, the kinds of data it provides, and their file format. The most onerous part is to extract meaningful things, which leads to development of **rGEO**

The next step is quite nature, convert that meaningful thing to GSEA input and run GSEA (thousands of, thanks to Prof. Xie's workstation). In the process I develop **qGSEA**, I think its significance fall far behind **rGEO**, though many user might love the former.

Finally, the most tough part, synthesis the result and draw conclusion. There is no precedent to follow, I have to make every decision by myself. Part of the process is presented in the thesis, but that is just the tip of the iceburg.




## Useful resource

- LEM4 is also know as ANKLE2 ([PDB link](https://www.uniprot.org/uniprot/Q86XL3#interaction))
- GSEA [User Guide](http://software.broadinstitute.org/gsea/doc/GSEAUserGuideFrame.html?_Preparing_Data_Files), [Data Format](http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats)
- [SOFT format](https://www.ncbi.nlm.nih.gov/geo/info/soft.html), [GEO search](https://www.ncbi.nlm.nih.gov/geo/browse/?view=series&display=20)



-----------------------
[![Creative Commons License](https://i.creativecommons.org/l/by-nc/4.0/88x31.png)](http://creativecommons.org/licenses/by-nc/4.0/)  
This work is licensed under a [Creative Commons Attribution-NonCommercial 4.0 International License](http://creativecommons.org/licenses/by-nc/4.0/)
