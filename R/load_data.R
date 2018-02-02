library(data.table)
# devtools::install('~/libs/maartenutils')
library(maartenutils)
library(dtplyr)
library(tidyverse)
p_root <- file.path('~/Projects/TONIC')

patient_labels <- read.csv(file.path(p_root, 'data-raw/patient_labels.csv'),
                           dec = ',', sep = ';') %>% as.data.table
normalize_colnames(patient_labels)
setnames(patient_labels, 'x', 'filename')
setnames(patient_labels, gsub('\\.', '_', colnames(patient_labels))) 
maartenutils::set_dt_types(patient_labels, 
                           c('mean_log2_hk' = 'numeric', 
                             'tis_score' = 'numeric'))
patient_labels[, patient := paste0('pat_', patient)]
setkey(patient_labels, filename)

exp_levels <- read.csv(file.path(p_root, 'data-raw/normalized_gene_expression.csv'),
                    # verbose = T,
                    dec = ',', 
                    sep = ';', skip = 1L) %>% as.data.table

setnames(exp_levels, gsub('.RCC$', '', colnames(exp_levels))) 
setnames(exp_levels, gsub('^X', '', colnames(exp_levels))) 
setnames(exp_levels, 'File.Name', 'gene_symbol')
setdiff(colnames(exp_levels), 'gene_symbol') %>% 
  { setNames(rep('numeric', length(.)), .) } %>%
  set_dt_types(exp_levels, .)
maartenutils::set_dt_types(exp_levels, 
                           c('mean_log2_hk' = 'numeric', 
                             'tis_score' = 'numeric'))

# length(intersect(patient_labels[, filename], colnames(exp_levels)))
# length(colnames(exp_levels))

stopifnot(all(patient_labels[, filename] %in% colnames(exp_levels)))
