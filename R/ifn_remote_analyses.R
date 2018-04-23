setwd('/DATA/users/m.slagter/TONIC')
source('R/init.R')

if (T && !maartenutils::local_run) {
  pacman::p_load('doParallel')
  doParallel::registerDoParallel(cores=16)
}

if (F) {
  library(knitr)
  knitr::knit(file.path(root_dir, 'ML.Rmd'))
}

if (T) {
  pacman::p_load(shinystan)
  pacman::p_load("rstan")
  if (!require(brms)) {
    devtools::install_github("paul-buerkner/brms")
  } 
  library(brms)
  plyr::l_ply(test_genes, run_brm_gene, .parallel = !local_run)
  plyr::l_ply(genesets, run_brm_geneset, .parallel = !local_run)
}

