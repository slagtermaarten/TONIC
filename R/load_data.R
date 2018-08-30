# library(edgeR)
# d <- DGEList(counts=yourcounts)
# d <- calcNormFactors(d)
# v <- voom(d, design)

if (T) {
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
}



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


vcfs <- list.files(file.path(forge_mirror, 'calls'), pattern = 'combined.vcf$') 
vcf_table <- data.table(vcf_fn = vcfs)
vcf_table[, 'tumor_cf' := gsub('.{2}_.{4}_.{1,2}_(CF\\d{5})_.*', 
                               '\\1', vcf_fn)]
vcf_table[, 'normal_cf' :=
          gsub('.{2}_.{4}_.{1,2}_CF\\d{5}_.{8}_vs_.{4}_.{1,2}_(CF\\d{5}).*', 
                               '\\1', vcf_fn)]
vcf_table[grepl('.vcf', normal_cf), 
          'normal_cf' :=
          gsub('.{2}_.{4}_.{1,2}_CF\\d{5}_.{8}_vs_.{4}_.{1,2}_(\\d{2}).*', 
                               '\\1', vcf_fn)]
wes_table <- vcf_table

# vcf_table[, 'cf' := gsub('.{2}_.{4}_.{1,2}_(CF\\d{5})_.*(CF\\d{5}).*', 
#                          '\\1-\\2', fn)]
# vcf_table[, 'tumor_cf' := gsub('(.*)-(.*)', '\\1', cf)]
# vcf_table[, 'normal_cf' := gsub('(.*)-(.*)', '\\2', cf)]
# vcf_table[, cf := NULL]

if (F) {
  contra_fns <- list.files(file.path(forge_mirror, 'contra', 'CNATable'),
                           pattern = '.txt$') 
  contra_fns <- data.table(contra_fn = contra_fns)
  contra_fns[, 'tumor_cf' := gsub('.{4}_.{1,2}_(CF\\d{5})_.*', 
                                  '\\1', contra_fn)]
  wes_table <- controlled_merge(contra_fns, vcf_table, by_cols = 'tumor_cf')
  rm(vcf_table)
  rm(contra_fns)
}

sequenza_fns <- list.files(file.path(forge_mirror, 'sequenza_plots', 'seqres'), 
                           pattern = 'segments.txt$') 
sequenza_fns <- data.table(sequenza_fn = sequenza_fns)
sequenza_fns[, 'tumor_cf' := gsub('.{1,2}_(CF\\d{5})_.*', '\\1', sequenza_fn)]
sequenza_fns[, sequenza_fn := file.path(forge_mirror, 'sequenza_plots', 'seqres', sequenza_fn)]
wes_table %<>% controlled_merge(sequenza_fns, dup_priority = 'a')
rm(sequenza_fns)

vcfs <- list.files(file.path(forge_mirror, 'haplotypecaller-q100'), 
                   pattern = 'q100.vcf$') 
vcf_table <- data.table(germline_vcf = vcfs)
vcf_table[, 'normal_cf' := gsub('.{2}_.{4}_.{1,2}_(CF\\d{5})_.*', 
                               '\\1', germline_vcf)]
vcf_table[grepl('\\.vcf$', normal_cf), 
          normal_cf := gsub('.{2}_.{4}_.{1,2}_(\\d{2})_.*', '\\1', 
                             germline_vcf)]
wes_table %<>% controlled_merge(vcf_table, by_cols = 'normal_cf', 
                                dup_priority = 'a')
rm(vcfs)
rm(vcf_table)
maartenutils::write_tsv(wes_table, 
                        file.path(p_root, 'ext', 'dnaseq_cf_numbers.tsv'))

# vcf_table[!grepl('.vcf', normal_cf)]

read_cibersort <- function(fn = file.path('data-raw', 'CIBERSORT.Output_Job2.csv')) {
  cibersort <- fread(fn) %>% normalize_colnames()
  setnames(cibersort, 'input_sample', 'cf_number')
  colnames(cibersort) <- gsub('_\\(tregs\\)', '', colnames(cibersort))
  cibersort[, cf_number := tolower(cf_number)]
  cibersort <-
    controlled_merge(cibersort, rna_sample_annotation[, .(cf_number, patient,
                                                        timepoint)],
                   by_cols = 'cf_number')
  cibersort <-
    controlled_merge(cibersort, 
                     rna_sample_annotation[, .(patient, arm, clinical_response)],
                     by_cols = 'patient')
  cibersort[, timepoint := factor(timepoint, levels = timepoints)]
  cibersort[, timepoint := droplevels(timepoint)]
  cibersort[, cf_number := NULL]
  # cibersort[is.na(clinical_response), clinical_response := 'NA']
  return(cibersort)
}

