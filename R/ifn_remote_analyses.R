root_dir <- '/DATA/users/m.slagter/TONIC'
setwd(root_dir)
source(file.path(root_dir, 'R', 'init.R'))

if (T && !maartenutils::local_run) {
  pacman::p_load('doParallel')
  doParallel::registerDoParallel(cores=16)
}

if (F) {
  library(knitr)
  knitr::knit(file.path(root_dir, 'ML.Rmd'))
}

if (F) {
  pacman::p_load(shinystan)
  pacman::p_load("rstan")
  if (!require(brms)) {
    devtools::install_github('paul-buerkner/brms')
  } 
  library(brms)
  plyr::l_ply(test_genes, run_brm_gene, .parallel = !local_run)
  plyr::l_ply(genesets, run_brm_geneset, .parallel = !local_run)
}

if (T) {
  # devtools::install_github('karthik/wesanderson')
  compute_tp_comp_FCs(tp1 = 'Baseline', tp2 = 'Post-induction', 
                      tp3 = '-2', tp4 = '0')
  compute_tp_comp_FCs(tp1 = 'Baseline', tp2 = 'On nivo', 
                      tp3 = '-2', tp4 = '6')
  compute_tp_comp_FCs(tp1 = 'Baseline', tp2 = 'On nivo', 
                      tp3 = '-2', tp4 = '10')
  compute_tp_comp_FCs(tp1 = 'Baseline', tp2 = 'On nivo', 
                      tp3 = '-2', tp4 = '12')
}

