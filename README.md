# IdeasCustom
The URL to the GitHub (i.e., the source code) is: 
https://github.com/TatiZhang/IdeasCustom

The URL to the Pkgdown webpage is: 
https://tatizhang.github.io/IdeasCustom/articles/IdeasCustom.html

## Purpose

The IdeasCustom package provides a tool for analyzing gene expression data on individual-level using IDEAS mehtod and decomposed components of Wasserstein-2 distance. It includes functions to arrange gene expression data by donors, initialize distance array lists, compute divergence metrics using the Wasserstein distance, and more. This package is designed to facilitate the exploration and analysis of complex gene expression datasets.

## Installation

You can install the IdeasCustom package from GitHub using the following commands:

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install IdeasCustom from GitHub
devtools::install_github("TatiZhang/IdeasCustom")
```
## Dependencies

## Session Info
```r
> devtools::session_info()
─ Session info ──────────────────────────
 setting  value
 version  R version 4.3.1 (2023-06-16)
 os       macOS Ventura 13.1
 system   x86_64, darwin20
 ui       RStudio
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       America/Los_Angeles
 date     2024-05-27
 rstudio  2023.12.0+369 Ocean Storm (desktop)
 pandoc   3.1.1 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/ (via rmarkdown)

─ Packages ──────────────────────────────
 ! package      * version    date (UTC) lib source
   askpass        1.2.0      2023-09-03 [1] CRAN (R 4.3.0)
   brio           1.1.4      2023-12-10 [1] CRAN (R 4.3.0)
   cachem         1.0.8      2023-05-01 [1] CRAN (R 4.3.0)
   callr          3.7.5      2024-02-19 [1] CRAN (R 4.3.2)
   cli            3.6.2      2023-12-11 [1] CRAN (R 4.3.0)
   codetools      0.2-19     2023-02-01 [1] CRAN (R 4.3.1)
   colorspace     2.1-0      2023-01-23 [1] CRAN (R 4.3.0)
   commonmark     1.9.0      2023-03-17 [1] CRAN (R 4.3.0)
   crayon         1.5.2      2022-09-29 [1] CRAN (R 4.3.0)
   credentials    2.0.1      2023-09-06 [1] CRAN (R 4.3.0)
   curl           5.2.0      2023-12-08 [1] CRAN (R 4.3.0)
   data.table   * 1.15.4     2024-03-30 [1] CRAN (R 4.3.2)
   desc           1.4.3      2023-12-10 [1] CRAN (R 4.3.0)
   devtools       2.4.5      2022-10-11 [1] CRAN (R 4.3.0)
   digest         0.6.34     2024-01-11 [1] CRAN (R 4.3.0)
   dlstats        0.1.7      2023-05-24 [1] CRAN (R 4.3.0)
   doRNG        * 1.8.6      2023-01-16 [1] CRAN (R 4.3.0)
   downlit        0.4.3      2023-06-29 [1] CRAN (R 4.3.0)
   dplyr          1.1.4      2023-11-17 [1] CRAN (R 4.3.0)
   ellipsis       0.3.2      2021-04-29 [1] CRAN (R 4.3.0)
   evaluate       0.23       2023-11-01 [1] CRAN (R 4.3.0)
   fansi          1.0.6      2023-12-08 [1] CRAN (R 4.3.0)
   fastmap        1.1.1      2023-02-24 [1] CRAN (R 4.3.0)
   foreach      * 1.5.2      2022-02-02 [1] CRAN (R 4.3.0)
   fs             1.6.3      2023-07-20 [1] CRAN (R 4.3.0)
   generics       0.1.3      2022-07-05 [1] CRAN (R 4.3.0)
   gert           2.0.1      2023-12-04 [1] CRAN (R 4.3.0)
   ggplot2        3.4.4      2023-10-12 [1] CRAN (R 4.3.0)
   gh             1.4.0      2023-02-22 [1] CRAN (R 4.3.0)
   gitcreds       0.1.2      2022-09-08 [1] CRAN (R 4.3.0)
   glue           1.7.0      2024-01-09 [1] CRAN (R 4.3.0)
   gtable         0.3.4      2023-08-21 [1] CRAN (R 4.3.0)
   htmltools      0.5.7      2023-11-03 [1] CRAN (R 4.3.0)
   htmlwidgets    1.6.4      2023-12-06 [1] CRAN (R 4.3.0)
   httpuv         1.6.13     2023-12-06 [1] CRAN (R 4.3.0)
   httr           1.4.7      2023-08-15 [1] CRAN (R 4.3.0)
   httr2          1.0.0      2023-11-14 [1] CRAN (R 4.3.0)
 P IdeasCustom  * 0.0.0.9002 2024-05-28 [?] load_all()
   iterators      1.0.14     2022-02-05 [1] CRAN (R 4.3.0)
   jsonlite       1.8.8      2023-12-04 [1] CRAN (R 4.3.0)
   knitr          1.45       2023-10-30 [1] CRAN (R 4.3.0)
   later          1.3.2      2023-12-06 [1] CRAN (R 4.3.0)
   lifecycle      1.0.4      2023-11-07 [1] CRAN (R 4.3.0)
   magrittr       2.0.3      2022-03-30 [1] CRAN (R 4.3.0)
   memoise        2.0.1      2021-11-26 [1] CRAN (R 4.3.0)
   mime           0.12       2021-09-28 [1] CRAN (R 4.3.0)
   miniUI         0.1.1.1    2018-05-18 [1] CRAN (R 4.3.0)
   munsell        0.5.0      2018-06-12 [1] CRAN (R 4.3.0)
   openssl        2.1.1      2023-09-25 [1] CRAN (R 4.3.0)
   pillar         1.9.0      2023-03-22 [1] CRAN (R 4.3.0)
   pkgbuild       1.4.3      2023-12-10 [1] CRAN (R 4.3.0)
   pkgconfig      2.0.3      2019-09-22 [1] CRAN (R 4.3.0)
   pkgdown        2.0.7      2022-12-14 [1] CRAN (R 4.3.0)
   pkgload        1.3.4      2024-01-16 [1] CRAN (R 4.3.0)
   prettyunits    1.2.0      2023-09-24 [1] CRAN (R 4.3.0)
   processx       3.8.3      2023-12-10 [1] CRAN (R 4.3.0)
   profvis        0.3.8      2023-05-02 [1] CRAN (R 4.3.0)
   promises       1.2.1      2023-08-10 [1] CRAN (R 4.3.0)
   ps             1.7.6      2024-01-18 [1] CRAN (R 4.3.0)
   purrr          1.0.2      2023-08-10 [1] CRAN (R 4.3.0)
   R6             2.5.1      2021-08-19 [1] CRAN (R 4.3.0)
   rappdirs       0.3.3      2021-01-31 [1] CRAN (R 4.3.0)
   rcmdcheck      1.4.0      2021-09-27 [1] CRAN (R 4.3.0)
   RColorBrewer   1.1-3      2022-04-03 [1] CRAN (R 4.3.0)
   Rcpp           1.0.12     2024-01-09 [1] CRAN (R 4.3.0)
   remotes        2.4.2.1    2023-07-18 [1] CRAN (R 4.3.0)
   rlang          1.1.3      2024-01-10 [1] CRAN (R 4.3.0)
   rmarkdown      2.25       2023-09-18 [1] CRAN (R 4.3.1)
   rngtools     * 1.5.2      2021-09-20 [1] CRAN (R 4.3.0)
   roxygen2       7.3.1      2024-01-22 [1] CRAN (R 4.3.2)
   rprojroot      2.0.4      2023-11-05 [1] CRAN (R 4.3.0)
   rstudioapi     0.15.0     2023-07-07 [1] CRAN (R 4.3.0)
   scales         1.3.0      2023-11-28 [1] CRAN (R 4.3.0)
   sessioninfo    1.2.2      2021-12-06 [1] CRAN (R 4.3.0)
   shiny          1.8.0      2023-11-17 [1] CRAN (R 4.3.0)
   stringi        1.8.3      2023-12-11 [1] CRAN (R 4.3.0)
   stringr        1.5.1      2023-11-14 [1] CRAN (R 4.3.0)
   sys            3.4.2      2023-05-23 [1] CRAN (R 4.3.0)
   testthat       3.2.1.1    2024-04-14 [1] CRAN (R 4.3.2)
   tibble         3.2.1      2023-03-20 [1] CRAN (R 4.3.0)
   tidyselect     1.2.0      2022-10-10 [1] CRAN (R 4.3.0)
   transport      0.15-2     2024-05-10 [1] CRAN (R 4.3.3)
   urlchecker     1.0.1      2021-11-30 [1] CRAN (R 4.3.0)
   usethis        2.2.2      2023-07-06 [1] CRAN (R 4.3.0)
   utf8           1.2.4      2023-10-22 [1] CRAN (R 4.3.0)
   vctrs          0.6.5      2023-12-01 [1] CRAN (R 4.3.0)
   whisker        0.4.1      2022-12-05 [1] CRAN (R 4.3.0)
   withr          3.0.0      2024-01-16 [1] CRAN (R 4.3.0)
   xfun           0.42       2024-02-08 [1] CRAN (R 4.3.2)
   xml2           1.3.6      2023-12-04 [1] CRAN (R 4.3.0)
   xopen          1.0.0      2018-09-17 [1] CRAN (R 4.3.0)
   xtable         1.8-4      2019-04-21 [1] CRAN (R 4.3.0)
   yaml           2.3.8      2023-12-11 [1] CRAN (R 4.3.0)

 [1] /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/library