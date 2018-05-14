# library(edgeR)
# d <- DGEList(counts=yourcounts)
# d <- calcNormFactors(d)
# v <- voom(d, design)

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
stopifnot(all(patient_labels[!is.na(filename), filename] %in% colnames(exp_levels)))


danaher_scores <- read.csv(file.path(p_root, 'data-raw', 
                                     'danaher_geneset_scores.csv'),
                           # verbose = T,
                           dec = ',',
                           sep = ';', skip = 0L) %>% as.data.table
setnames(danaher_scores, tolower(colnames(danaher_scores)))

## Ensure everything's numeric
setdiff(colnames(danaher_scores), c('patient', 'arm', 'response')) %>%
  { setNames(rep('numeric', length(.)), .) } %>%
  set_dt_types(danaher_scores, .)
danaher_scores.m <- melt(danaher_scores,
                         id.vars = c('patient', 'arm', 'response'))
danaher_scores.m[, 'timepoint' := gsub('.*_(.*)', '\\1', variable)]
danaher_scores.m[, variable := gsub('(.*)_.*', '\\1', variable)]
danaher_scores.m[, patient := paste0('pat_', patient)]
danaher_scores.m[, timepoint := factor(plyr::mapvalues(timepoint, 
                                                       unique(timepoint), 
                                                       timepoints),
                                       levels = timepoints)]
danaher_scores.m[arm == 'Cyclofosfamide',  arm := 'Cyclophosphamide']
# danaher_scores.m[, table(arm)]
set_dt_types(danaher_scores.m, c('patient' = 'factor', 'timepoint' = 'factor',
                                 'variable' = 'factor'))
danaher_scores.m <- controlled_merge(danaher_scores.m,
                                     patient_labels[, .(patient, response, arm)])
danaher_scores.m[, arm := factor(arm, levels = treatment_arms)]
setnames(danaher_scores.m, 'response', 'clinical_response')
danaher_scores.m[, clinical_response := factor(clinical_response, 
                                               levels = c('NR', 'R'))]
                                     

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
# danaher_signature_genes[, table(signature_name)]
# danaher_signature_genes[, table(gene_name)]


## Genes to test in comparative analyses
test_genes <- exp_levels[, unique(gene_symbol)] %>%
  grep(pattern = 'POS_\\w', x = ., value = T, invert = T) %>%
  grep(pattern = 'NEG_\\w', x = ., value = T, invert = T) %>%
  grep(pattern = 'ERCC Controls', x = ., value = T, invert = T) %>%
  grep(pattern = 'Housekeeping', x = ., value = T, invert = T)

test_gene_sets <- danaher_scores.m[, levels(variable)]


## Read adaptive seq data if required
read_adaptive_seqs <- function(force_reload = F) {
  if (exists('arr', envir = globalenv()) && !force_reload) {
    return(invisible()) 
  }
  arr <- w_fread(file.path(p_root, 'data-raw', 'AdaptiveAllRearrangements.tsv'),
                 col_classes = c('productive_frequency' = 'numeric'))
  arr[, productive_frequency := as.numeric(productive_frequency)]
  arr[, reads := NULL]
  arr <- arr[!is.na(amino_acid) & amino_acid != 'na']
  setnames(arr, 'sample_name', 'adaptive_sample_name')
  arr1 <-
    controlled_merge(arr[grepl('T_', adaptive_sample_name)],
                     patient_labels[, .(patient, adaptive_sample_name,
                                        timepoint, arm, clinical_response,
                                        response)],
                     by_cols = 'adaptive_sample_name',
                     all.x = T, all.y = F)
  arr2 <-
    controlled_merge(arr[grepl('B_', adaptive_sample_name)],
                     blood_adaptive[, .(patient, adaptive_sample_name, arm,
                                        blood_timepoint,
                                        clinical_response, response)],
                     by_cols = 'adaptive_sample_name',
                     all.x = T, all.y = F)
  setnames(arr2, 'blood_timepoint', 'timepoint')
  arr <- rbind(arr1, arr2, fill = T)[patient != 'patient_64']
  
  arr[, timepoint := factor(as.character(timepoint), 
                            levels = c(timepoints, blood_timepoints))]

  ## Merge in 'Fraction T-cells of nucleated cells' 
  arr <- 
    controlled_merge(arr, 
                     patient_labels[, .(adaptive_sample_name, 
                                        adaptive_t_cells)],
                     by_cols = 'adaptive_sample_name')
  arr <- 
    controlled_merge(arr, 
                     blood_adaptive[, .(adaptive_sample_name, 
                                        adaptive_t_cells)],
                     by_cols = 'adaptive_sample_name')

  arr[, 'normalized_frequency' := productive_frequency * adaptive_t_cells]
  assign('arr', arr, envir = globalenv())
  return(invisible())
}

vcfs <- list.files(file.path(forge_mirror, 'calls'), pattern = '.vcf$') 
vcf_table <- data.table(vcf_fn = vcfs)
vcf_table[, 'tumor_cf' := gsub('.{2}_.{4}_.{1,2}_(CF\\d{5})_.*', 
                               '\\1', vcf_fn)]

# vcf_table[, 'cf' := gsub('.{2}_.{4}_.{1,2}_(CF\\d{5})_.*(CF\\d{5}).*', 
#                          '\\1-\\2', fn)]
# vcf_table[, 'tumor_cf' := gsub('(.*)-(.*)', '\\1', cf)]
# vcf_table[, 'normal_cf' := gsub('(.*)-(.*)', '\\2', cf)]
# vcf_table[, cf := NULL]

contra_fns <- list.files(file.path(forge_mirror, 'contra', 'CNATable'), 
                   pattern = '.txt$') 
contra_fns <- data.table(contra_fn = contra_fns)
contra_fns[, 'tumor_cf' := gsub('.{4}_.{1,2}_(CF\\d{5})_.*', 
                                '\\1', contra_fn)]
wes_table <- merge(contra_fns, vcf_table, all = T)
rm(vcf_table)
rm(contra_fns)


#' Filter out patients that have dubious annotation
#'
#'
filter_patients <- name <- function(p_dat, ...) { 
  comb_vars <- as.character(...)
  clinical_params <- c('response', 'clinical_response')
  if (is.null(comb_vars)) return(p_dat)
  if (any(comb_vars %in% clinical_params)) {
    # resp_var <- comb_vars[which(comb_vars %in% clinical_params)]
    p_dat <- p_dat[patient %nin% c('pat_63', 'pat_64')]
  }
  return(p_dat)
}
