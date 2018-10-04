# library(maartenutils)
devtools::load_all('~/libs/maartenutils')
# pacman::p_load('extrafont')
options(error = traceback)
options(warn = 1)
force_installation <- F
if (F && !local_run) {
  source(file.path('~/antigenic_space', 'bin', 'install_packages.R'))
  ## No idea why this is necessary
  # devtools::install('~/libs/maartenutils')
  # devtools::install('~/libs/GSEAgenesets')
}
suppressMessages(pacman::p_load(data.table))
library(magrittr)
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
rds_dir <- file.path(p_root, 'rds')
dir.create(rds_dir, showWarnings = F)
plot_dir <- file.path(p_root, 'plots2')
img_dir <- plot_dir
data_dir <- file.path(p_root, 'data-raw')
forge_mirror <- file.path(data_dir, 'forge', 'userdata', 
                          'sLpLDXtXQNrqruxTWwnV2GSuNVSSGAdgP60YoQzR', 
                          'TONIC_stage1_WES')
# source('R/read_xlsx.R')
cond_rm('blood_adaptive')
sr(blood_adaptive)
sr(patient_labels)
source('R/constants.R')

source('R/misc.R')
source('R/load_data.R')
source('R/load_cache.R')
source('R/read_DNASeq.R')
source('R/load_adaptive.R')
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
source('R/plot_DNASeq.R')

source('R/GSEA_plotting.R')
source('R/GSEA_funcs_paired.R')
# library('extrafont')
# font_import(pattern="[A/a]rial", prompt=FALSE)
library(cowplot)
library(naturalsort)
# ggplot2::theme_set(theme_ms(base_size = 8,
#                             text=element_text(size=8, family='Arial')))
ggplot2::theme_set(theme_ms(base_size = 10, 
  panel.border = ggplot2::element_rect(colour = 'grey20', fill = NA, size = 1,
                                       linetype = 'solid')))

patient_labels %<>% 
  controlled_merge(wes_table[, .('cf_number' = tumor_cf, brca1_like)])
