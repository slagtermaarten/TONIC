# library(maartenutils)
packrat::extlib('devtools')
# devtools::install('~/libs/maartenutils')
if (!require(maartenutils)) {
  devtools::install('~/libs/maartenutils')
}
devtools::load_all('~/libs/maartenutils')
if (!require(genesets)) {
  devtools::install('~/libs/genesets')
}
devtools::load_all('~/libs/genesets')
if (!require(flexgsea)) {
  tryCatch(devtools::install_github('NKI-CCB/ggsea'), 
    error = function(e) { print(e); }) 
} 
tryCatch(library(flexgsea), error = function(e) { print(e) }) 
# library('extrafont')
options(error = traceback)
options(warn = 1)
force_installation <- F
if (F && !local_run) {
  source(file.path('~/antigenic_space', 'bin', 'install_packages.R'))
  ## No idea why this is necessary
  # devtools::install('~/libs/maartenutils')
}
# suppressMessages(library(data.table))
library(data.table)
library(magrittr)
library(dtplyr)
library(dplyr)
library('tidyverse')
library(magrittr)
library(scales)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!require(limma))
  BiocManager::install("limma")
if (!require(edgeR))
  BiocManager::install("edgeR")
library(limma)
library(edgeR)
p_root <- path.expand(file.path('~/TONIC'))
setwd(p_root)
rds_dir <- file.path(p_root, 'rds')
dir.create(rds_dir, showWarnings = F)
# if (!require(serializer)) {
#   devtools::install('~/libs/serializer')
# }
devtools::load_all('~/libs/serializer')
sr <- serializer::gen_serializer(rds_dir)
plot_dir <- file.path(p_root, 'plots2')
img_dir <- plot_dir
data_dir <- file.path(p_root, 'data-raw')
forge_mirror <- file.path(data_dir, 'forge', 'userdata', 
  'sLpLDXtXQNrqruxTWwnV2GSuNVSSGAdgP60YoQzR', 
  'TONIC_stage1_WES')
# source('R/read_xlsx.R')
cond_rm('blood_adaptive')
# sr(blood_adaptive)
blood_adaptive <- readRDS(file.path(rds_dir, 'blood_adaptive.rds'))
# sr(patient_labels)
patient_labels <- readRDS(file.path(rds_dir, 'patient_labels.rds'))
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
# library(cowplot)
# install.packages('naturalsort')
library(naturalsort)
# ggplot2::theme_set(theme_ms(base_size = 8,
#                             text=element_text(size=8, family='Arial')))
ggplot2::theme_set(theme_ms(base_size = 8, 
  panel.border = ggplot2::element_rect(colour = 'grey20', fill = NA, size = 1,
                                       linetype = 'solid')))
