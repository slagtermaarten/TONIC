hostname <- Sys.info()[["nodename"]]
servers <- c('paranoid', 'steroid', 'medoid', 'void', 'coley')
local_run <- !(hostname %in% servers)
dir.create('rds2', showWarnings = F)
source('~/antigenic_space/bin/install_packages.R')


if (!local_run) {
  root_dir <- '/DATA/users/m.slagter/TONIC'
} else {
  root_dir <- '~/Projects/TONIC'
}
setwd(root_dir)
source(file.path(root_dir, 'R', 'init.R'))


if (T && !maartenutils::local_run) {
  library('doParallel')
  doParallel::registerDoParallel(cores = 8)
}


if (F) {
  library(knitr)
  knitr::knit(file.path(root_dir, 'ML.Rmd'))
}


if (F) {
  library(shinystan)
  library("rstan")
  if (!require(brms)) {
    devtools::install_github('paul-buerkner/brms')
  } 
  library(brms)
  plyr::l_ply(test_genes, run_brm_gene, .parallel = !local_run)
  plyr::l_ply(genesets, run_brm_geneset, .parallel = !local_run)
}


if (F) {
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


if (F) {
  ## GSEA analyses
  source(file.path(root_dir, 'R', 'init.R'))
  source('R/load_rna_dat.R')
  source('R/GSEA_funcs_paired.R')
  # res <- gsea_all_arms(resp_exp = 'timepoint', nperm = 1e3,
  #                      abs = F, fn_extra = '')
  # debugonce(gsea_all_arms)
  res <- gsea_all_arms(gene_sets = HALLMARK_pathways,
                       patients = rna_sample_annotation[, unique(patient)],
                       allowed_timepoints = c('Baseline', 'Post-induction'),
                       gene_score_fn = my_paired_WC_test,
                       resp_exp = 'timepoint',
                       paired_test = T, 
                       nperm = 1e3,
                       abs = F,
                       fn_extra = '')
  saveRDS(res, file = 'rds2/GSEA_paired.rds')
}


if (F) {
  source(file.path(root_dir, 'R', 'init.R'))
  source('R/load_rna_dat.R')
  source('R/GSEA_funcs_paired.R')
  res <- gsea_all_arms(gene_sets = HALLMARK_pathways,
                       patients = rna_sample_annotation[!is.na(clinical_response), 
                                                        unique(patient)],
                       allowed_timepoints = 'Baseline',
                       gene_score_fn = my_unpaired_WC_test,
                       resp_exp = 'clinical_response',
                       paired_test = F, 
                       nperm = 1e3,
                       abs = F,
                       fn_extra = '')
  saveRDS(res, file = 'rds2/GSEA_unpaired_baseline.rds')
}


if (F) {
  source(file.path(root_dir, 'R', 'init.R'))
  source('R/load_rna_dat.R')
  source('R/GSEA_funcs_paired.R')
  res <- gsea_all_arms(gene_sets = HALLMARK_pathways,
                       patients = rna_sample_annotation[!is.na(clinical_response), 
                                                        unique(patient)],
                       allowed_timepoints = 'Post-induction',
                       gene_score_fn = my_unpaired_WC_test,
                       resp_exp = 'clinical_response',
                       paired_test = F, 
                       nperm = 1e3,
                       abs = F,
                       fn_extra = '')
  saveRDS(res, file = 'rds2/GSEA_unpaired_postinduction.rds')
}


if (T) {
  source(file.path(root_dir, 'R', 'init.R'))
  source('R/load_rna_dat.R')
  source('R/GSEA_funcs_paired.R')
  res <- gsea_wrapper(gene_sets = STING_pathways,
                      patients = rna_sample_annotation[!is.na(clinical_response), 
                                                       unique(patient)],
                      allowed_timepoints = 'Baseline',
                      gene_score_fn = my_unpaired_WC_test,
                      resp_exp = 'clinical_response',
                      paired_test = F, 
                      nperm = 1e3,
                      abs = F,
                      fn_extra = '')
  saveRDS(res, file = 'rds2/GSEA_unpaired_baseline_STING.rds')
  file.exists('rds2/GSEA_unpaired_baseline_STING.rds')
}


if (T) {
  source(file.path(root_dir, 'R', 'init.R'))
  source('R/load_rna_dat.R')
  source('R/GSEA_funcs_paired.R')
  res <- gsea_wrapper(gene_sets = STING_pathways,
                      patients = rna_sample_annotation[!is.na(clinical_response), 
                                                       unique(patient)],
                      allowed_timepoints = 'Post-induction',
                      gene_score_fn = my_unpaired_WC_test,
                      resp_exp = 'clinical_response',
                      paired_test = F, 
                      nperm = 1e3,
                      abs = F,
                      fn_extra = '')
  saveRDS(res, file = 'rds2/GSEA_unpaired_postinduction_STING.rds')
  file.exists('rds2/GSEA_unpaired_postinduction_STING.rds')
}
