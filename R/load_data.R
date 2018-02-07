pacman::p_load(data.table)
# devtools::install('~/libs/maartenutils')
pacman::p_load(maartenutils)
pacman::p_load(dtplyr)
pacman::p_load(tidyverse)
p_root <- file.path('~/Projects/TONIC')
img_dir <- file.path(p_root, 'plots')

timepoints <- c('baseline' = 'Baseline', 'post.induction' = 'Post-induction',
                'on.nivo' = 'On nivo')

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

response_data <- read.csv(file.path(p_root, 'data-raw/response_data.csv'),
                          dec = ',', sep = ';') %>% as.data.table %>%
  normalize_colnames %>%
  ## Select last columns containing patient_id and RECIST labels
  { rev(.)[, c(2, 3)] } %>%
  { .[!is.na(patient.id)] } %>%
  setnames(., c('response', 'patient')) %>%
  ## Take out patients with dubious responses
  { .[!patient %in% c(61, 63, 64)] } 

response_data[, patient := as.character(paste0('pat_', patient))]

patient_labels[, response := NULL]
patient_labels[, patient := as.character(patient)]
patient_labels <- merge(patient_labels, response_data, by = 'patient')

exp_levels <- read.csv(file.path(p_root, 
                                 'data-raw/normalized_gene_expression.csv'),
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


danaher_scores <- read.csv(file.path(p_root,
                                     'data-raw/danaher_geneset_scores.csv'),
                    # verbose = T,
                    dec = ',', 
                    sep = ';', skip = 0L) %>% as.data.table
normalize_colnames(danaher_scores)

## Ensure everything's numeric
setdiff(colnames(danaher_scores), c('patient', 'arm', 'response')) %>% 
  { setNames(rep('numeric', length(.)), .) } %>%
  set_dt_types(danaher_scores, .)
danaher_scores.m <- melt(danaher_scores, 
                         id.vars = c('patient', 'arm', 'response'))
danaher_scores.m[, 'timepoint' := gsub('.*_(.*)', '\\1', variable)]
danaher_scores.m[, variable := gsub('(.*)_.*', '\\1', variable)]
danaher_scores.m[, patient := paste0('pat_', patient)]
danaher_scores.m[, timepoint := factor(timepoint, levels = names(timepoints))]
danaher_scores.m[, 'arm' := factor(arm, levels = levels(arm)[c(4,1,2,3,5)])]
set_dt_types(danaher_scores.m, c('patient' = 'factor', 'timepoint' = 'factor',
                                 'variable' = 'factor'))
# danaher_scores.m <- merge(danaher_scores.m, 
#                           patient_labels[, .(patient, response, arm)], 
#                           by = 'patient')


## Gene signatures
danaher_signature_genes <- read.csv(file.path(p_root, 
                                              'data-raw/danaher_signature_genes.csv'),
                           dec = ',', sep = ';') %>% as.data.table %>%
  normalize_colnames()
setnames(danaher_signature_genes, 
         gsub('\\.', '_', colnames(danaher_signature_genes)))
danaher_signature_genes[, 'signature_name' := gsub(' ', '_',
                                                   tolower(signature_name))] 
danaher_signature_genes[, 'signature_name' := gsub('\\.|-', '_',
                                                   tolower(signature_name))] 
danaher_signature_genes[, table(signature_name)]
danaher_signature_genes[, table(gene_name)]
