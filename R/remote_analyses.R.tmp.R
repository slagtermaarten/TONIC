hostname <- Sys.info()[["nodename"]]
servers <- c('paranoid', 'steroid', 'medoid', 'void', 'coley')
local_run <- !(hostname %in% servers)

if (!local_run) {
  root_dir <- '/DATA/users/m.slagter/TONIC'
} else {
  root_dir <- '~/Projects/TONIC'
}
setwd(root_dir)
source(file.path(root_dir, 'R', 'init.R'))

if (T && !maartenutils::local_run) {
  pacman::p_load('doParallel')
  doParallel::registerDoParallel(cores = 24)
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
  source('R/load_rna_dat.R')
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

  # source(file.path(root_dir, 'R', 'init.R'))
  # source('R/load_rna_dat.R')
  # source('R/GSEA_funcs_paired.R')
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

  # source(file.path(root_dir, 'R', 'init.R'))
  # source('R/load_rna_dat.R')
  # source('R/GSEA_funcs_paired.R')
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
  NPERM <- 1000
  if (maartenutils::local_run) {
    NPERM <- 3
  }
  # rownames(tpms_salmon)
  # rownames(rna_read_counts_salmon)
  # rownames(rna_read_counts_salmon_tmm)
  # colnames(tpms_salmon)
  # colnames(rna_read_counts_salmon_tmm)
  # class(rna_read_counts_salmon_tmm)
  # class(tpms_salmon)
  # dim(rna_read_counts_salmon_tmm)
  # dim(tpms_salmon)
  ## GSEA analyses
  if (!maartenutils::local_run) {
    source('R/load_rna_dat.R')
  }
  # debugonce(my_paired_WC_test)
  # which(is.na(colnames(rna_read_counts_salmon_tmm)))
  if (F) {
    res <- gsea_all_arms(gene_sets = HALLMARK_pathways,
                         patients = rna_sample_annotation[, unique(patient)],
                         allowed_timepoints = c('Baseline', 'Post-induction'),
                         gene_score_fn = my_paired_WC_test,
                         exp_mat = rna_read_counts_salmon_tmm_M,
                         # exp_mat = tpms_salmon,
                         resp_exp = 'timepoint',
                         paired_test = T, 
                         nperm = NPERM,
                         abs = F,
                         fn_extra = '')
    if (!local_run) {
      saveRDS(res, file = 'rds2/GSEA_paired_TMM.rds')
    }
  }

  source('R/GSEA_funcs_paired.R')
  res <- gsea_wrapper(gene_sets = HALLMARK_pathways,
                      patients = rna_sample_annotation[!is.na(clinical_response), 
                                                       unique(patient)],
                      allowed_timepoints = 'Baseline',
                      exp_mat = rna_read_counts_salmon_tmm_M,
                      gene_score_fn = my_unpaired_WC_test,
                      resp_exp = 'clinical_response',
                      paired_test = F, 
                      nperm = NPERM,
                      abs = F,
                      fn_extra = '')
  if (!local_run) {
    saveRDS(res, file = 'rds2/GSEA_unpaired_baseline_TMM.rds')
  }

  res <- gsea_wrapper(gene_sets = HALLMARK_pathways,
                      patients = rna_sample_annotation[!is.na(clinical_response), 
                                                       unique(patient)],
                      allowed_timepoints = 'Post-induction',
                      gene_score_fn = my_unpaired_WC_test,
                      resp_exp = 'clinical_response',
                      exp_mat = rna_read_counts_salmon_tmm_M,
                      paired_test = F, 
                      nperm = NPERM,
                      abs = F,
                      fn_extra = '')
  if (!local_run) {
    saveRDS(res, file = 'rds2/GSEA_unpaired_postinduction_TMM.rds')
  }
}
