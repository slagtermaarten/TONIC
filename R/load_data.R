# library(edgeR)
# d <- DGEList(counts=yourcounts)
# d <- calcNormFactors(d)
# v <- voom(d, design)

timepoints <- c('baseline' = 'Baseline', 'post.induction' = 'Post-induction',
                'on.nivo' = 'On nivo')
timepoints_inv <- setNames(names(timepoints), timepoints)


clinical_annotation <- read.csv(file.path(p_root,
                         'data-raw/180222_clinical_data_all_patients.csv'),
                         dec = ',', sep = ';') %>% as.data.table %>%
  normalize_colnames()
setnames(clinical_annotation, gsub('\\.', '_', colnames(clinical_annotation)))
setnames(clinical_annotation, gsub('_+', '_', colnames(clinical_annotation)))
setnames(clinical_annotation, gsub('_$', '', colnames(clinical_annotation)))
clinical_annotation[, 'patient' := paste0('pat_', study_id)]
clinical_annotation[, cd8_mm2 := gsub(',', '.', cd8_mm2)]
clinical_annotation[, stil := gsub(',', '.', stil)]
setnames(clinical_annotation, 'stil', 's_til')
maartenutils::set_dt_types(clinical_annotation,
                           c('cd8_mm2' = 'numeric',
                             's_til' = 'numeric'))
clinical_annotation[, induction_therapy := factor(induction_therapy,
  levels = c("No induction", "Cisplatin", "Cyclofosfamide", "Doxorubicin",
             "Radiotherapy"))]
clinical_annotation[, 'response_bin' := as.integer(response %in% c('CR', 'PR'))]
clinical_annotation <- clinical_annotation[patient != 'pat_NA']

patient_labels <- read.csv(file.path(p_root, 'data-raw/patient_labels.csv'),
                           dec = ',', sep = ';') %>% as.data.table %>%
  normalize_colnames()
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


## Merge Adaptive TCR abundance measures
## first merge adaptive sample IDs...
pacman::p_load(readxl)
pacman::p_load(naturalsort)
adaptive_sample_annotation <-
  read_excel(file.path(data_dir, 'SampleManifest_adaptive.xlsx')) %>%
  as.data.table %>%
  maartenutils::normalize_colnames()

# adaptive_sample_annotation[, .('adaptive_sample_name' = sample_name, 'cf_number' = cf_nummer)]
patient_labels <- merge(patient_labels,
  adaptive_sample_annotation[, .('adaptive_sample_name' = sample_name,
                                 'cf_number' = cf_nummer,
                                 'patient' = paste0('pat_', study_id),
                                 'timepoint' =
                                   timepoints[as.integer(gsub('N15TON-(\\d)',
                                                              '\\1',
                                                              tijdspunt))])],
                        by = c('patient', 'timepoint'))
## ... then merge summary stats
tmp <- read_tsv(file = file.path(data_dir, 'samples.diversity.txt')) %>%
  as.data.table %>%
  maartenutils::normalize_colnames()
setnames(tmp, gsub('\'', '', colnames(tmp)))
# merge(patient_labels, tmp, by = c('patient', 'timepoint'))

setdiff(colnames(tmp), c('sample', 'efron_thisted_estimator',
                         'daley_smith_estimator', 'ichao1')) %>%
  { setNames(rep('numeric', length(.)), .) } %>%
  set_dt_types(tmp, .)
setnames(tmp, 'sample', 'adaptive_sample_name')
## Finally merge in adaptive info
patient_labels <- merge(patient_labels, tmp, by = 'adaptive_sample_name')
rm(tmp)
patient_labels[, timepoint := factor(timepoint, levels = timepoints)]

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
                    sep = ';', skip = 0L) %>% as.data.table %>%
  normalize_colnames()

## Ensure everything's numeric
setdiff(colnames(danaher_scores), c('patient', 'arm', 'response')) %>%
  { setNames(rep('numeric', length(.)), .) } %>%
  set_dt_types(danaher_scores, .)
danaher_scores.m <- melt(danaher_scores,
                         id.vars = c('patient', 'arm', 'response'))
danaher_scores.m[, 'timepoint' := gsub('.*_(.*)', '\\1', variable)]
danaher_scores.m[, variable := gsub('(.*)_.*', '\\1', variable)]
danaher_scores.m[, patient := paste0('pat_', patient)]
danaher_scores.m[, timepoint := factor(timepoints[as.character(timepoint)],
                                       levels = timepoints)]
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


## Genes to test in comparative analyses
test_genes <- exp_levels[, unique(gene_symbol)] %>%
  grep(pattern = 'POS_\\w', x = ., value = T, invert = T) %>%
  grep(pattern = 'NEG_\\w', x = ., value = T, invert = T) %>%
  grep(pattern = 'ERCC Controls', x = ., value = T, invert = T) %>%
  grep(pattern = 'Housekeeping', x = ., value = T, invert = T)

test_gene_sets <- danaher_scores.m[, levels(variable)]
