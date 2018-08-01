# library(maartenutils)
devtools::load_all('~/libs/maartenutils')
# pacman::p_load('extrafont')
options(error = traceback)
options(warn=1)
force_installation <- F
if (!local_run) {
  source(file.path('~/antigenic_space', 'bin', 'install_packages.R'))
  ## No idea why this is necessary
  # devtools::install('~/libs/maartenutils')
  devtools::install('~/libs/GSEAgenesets')
}
pacman::p_load(data.table)
pacman::p_load(dtplyr)
pacman::p_load(dplyr)
pacman::p_load(tidyverse)
if (!require(ggsea)) {
  devtools::install_github('NKI-CCB/ggsea')
} else {
  library(ggsea)
}
library(limma)
library(edgeR)
# devtools::load_all('~/libs/GSEAgenesets')
# devtools::load_all('~/libs/serializer')
library(maartenutils)
library(GSEAgenesets)
library(serializer)
if (local_run) {
  p_root <- file.path('/Users/maartenslagter/surfdrive/Projects/TONIC')
} else {
  p_root <- file.path('/DATA/users/m.slagter/TONIC')
}
setwd(p_root)
sr <- serializer::gen_serializer(file.path(p_root, 'rds'))
img_dir <- file.path(p_root, 'plots')
rds_dir <- file.path(p_root, 'rds')
dir.create('rds2', showWarnings = F)
plot_dir <- file.path(p_root, 'plots')
data_dir <- file.path(p_root, 'data-raw')
forge_mirror <- file.path(data_dir, 'forge', 'userdata', 
                          'sLpLDXtXQNrqruxTWwnV2GSuNVSSGAdgP60YoQzR', 
                          'TONIC_stage1_WES')
# source('R/read_xlsx.R')
source('R/constants.R')
sr(blood_adaptive)
sr(patient_labels)
source('R/load_data.R')
source('R/load_cache.R')
source('R/read_DNASeq.R')
source('R/seq_stat_overview.R')

source('R/plotting_nanostring.R')
source('R/plotting_rna.R')
source('R/plotting_adaptive.R')

source('R/ML.R')
source('R/bayes.R')
source('R/test_associations.R')
source('R/ton_percentile_comp.R')
source('R/rna.R')
source('R/ML_RNASeq.R')

source('R/GSEA_plotting.R')
source('R/GSEA_funcs_paired.R')
# library('extrafont')
# font_import(pattern="[A/a]rial", prompt=FALSE)
library(cowplot)
# ggplot2::theme_set(theme_ms(base_size = 8,
#                             text=element_text(size=8, family='Arial')))
ggplot2::theme_set(theme_ms(base_size = 10))


#' Filter out patients that have dubious annotation
#'
#'
filter_patients <- name <- function(p_dat, ...) {
  comb_vars <- as.character(...)
  clinical_params <- c('response', 'clinical_response', 'comb_time_resp')
  if (is.null(comb_vars)) return(p_dat)
  if (any(comb_vars %in% clinical_params)) {
    # resp_var <- comb_vars[which(comb_vars %in% clinical_params)]
    p_dat <- p_dat[patient %nin% c('pat_63', 'pat_64')]
    for (varn in clinical_params) {
      if (varn %in% colnames(p_dat)) {
        N <- p_dat[is.na(get(varn)), .N]
        if (N > 0) {
          messagef('Removing %d donors due to absence of %s', N, varn)
          p_dat <- p_dat[!is.na(get(varn))]
        }
      }
    }
  }
  if (F) {
    pres_cols <-
      intersect(colnames(p_dat), c('patient', 'arm', 'timepoint',
                                   'blood_timepoint',
                                   'filename', 'adaptive_sample_name'))
    p_dat <- unique(p_dat, by = pres_cols)
  }
  return(p_dat)
}

