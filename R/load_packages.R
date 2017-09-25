devtools::load_all(file.path('~/antigenic_space', 'libs', 'fasanalysis'))
theme_set(theme_fas(base_size = 8))
pacman::p_load('data.table')
pacman::p_load('ggplot2')
pacman::p_load('rlm')
pacman::p_load('dplyr')
pacman::p_load('ggrepel')
pacman::p_load('plyr')
pacman::p_load('NMF')
devtools::load_all(file.path('~/libs', 'serializer'))
sf <- gen_serializer('rds')

# load_bioconductor <- function(pkg_name = 'DESeq2') {
#   if (!require(pkg_name)) {
#     source('https://bioconductor.org/biocLite.R')
#     biocLite(pkg_name)
#   }
#   library(pkg_name)
# }
# load_bioconductor('DESeq2')
# load_bioconductor('mygene')
if (!require(mygene)) {
  source('https://bioconductor.org/biocLite.R')
  biocLite(mygene)
}
library(mygene)
if (!require(DESeq2)) {
  source('https://bioconductor.org/biocLite.R')
  biocLite(DESeq2)
}
library(DESeq2)

if (!require('ggsea')) {
  devtools::install_github("NKI-CCB/ggsea")
  # devtools::install("NKI-CCB/ggsea")
  # devtools::install('/home/m.slagter/libs/ggsea')
  library("ggsea")
}
library('DESeq2')
