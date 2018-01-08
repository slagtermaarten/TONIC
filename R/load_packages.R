CRAN_packs <- c('praise', 'plyr', 'dtplyr', 'data.table', 'ggplot2', 'rlm',
                'ggrepel', 'NMF', 'testthat', 'rex', 'covr')
for (f in CRAN_packs) {
  if (!f %in% installed.packages()[, 1]) {
    install.packages(f, dependencies = c('Depends', 'Suggests'))
  }
  library(package = f, character.only = T)
}


github_packages <- c('decisionpatterns/operator.tools',
                     'decisionpatterns/formula.tools',
                     'slagtermaarten/serializer',
                     'slagtermaarten/maartenutils',
                     'NKI-CCB/ggsea')
devtools::install_github(github_packages, 
                         dependencies = c('Depends', 'Suggests'))
# devtools::install_github('NKI-CCB/ggsea', force = T)
library(serializer)
sf <- gen_serializer('rds')
## Load most bleeding edge version of maartenutils
devtools::load_all('~/libs/maartenutils')
library(ggsea)
# devtools::load_all('~/libs/ggsea')

load_bioconductor <- function(pkg_name = 'DESeq2') {
  if (!require(pkg_name, character.only = T)) {
    source('https://bioconductor.org/biocLite.R')
    biocLite(pkg_name)
  }
  library(pkg_name, character.only = T)
}
# load_bioconductor('DESeq2')
# load_bioconductor('mygene')
