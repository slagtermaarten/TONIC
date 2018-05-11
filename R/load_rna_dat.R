pacman::p_load(readxl)
pacman::p_load(naturalsort)

rna_sample_annotation <-
  read_excel(file.path(data_dir, 'Sample_annotations_RNAseq_TONICstageI.xlsx'))
setDT(rna_sample_annotation)
normalize_colnames(rna_sample_annotation)
## Get rid of last two rows
rna_sample_annotation <- rna_sample_annotation[-nrow(rna_sample_annotation)]
rna_sample_annotation <- rna_sample_annotation[-nrow(rna_sample_annotation)]
rna_sample_annotation[, timepoint_biopsy :=
                      gsub('.*-()', '\\1', timepoint_biopsy)]
rna_sample_annotation[, timepoint_biopsy :=
                      factor(timepoints[as.integer(timepoint_biopsy)],
                             levels = timepoints)]
rna_sample_annotation[, 'patient' := paste0('pat_', study_id)]
rna_sample_annotation[, patient_id := NULL]
rna_sample_annotation[, sequencing_number := NULL]
rna_sample_annotation[, study_id := NULL]
rna_sample_annotation[, `t-number` := NULL]
rna_sample_annotation <- rna_sample_annotation[patient != 'pat_NA']
setkey(rna_sample_annotation, patient)
setnames(rna_sample_annotation, 'timepoint_biopsy', 'timepoint')
rna_sample_annotation <- rna_sample_annotation[naturalsort::naturalorder(patient)]
rna_sample_annotation[, cf_number := tolower(cf_number)]
rna_sample_annotation <- merge(rna_sample_annotation, 
      unique(patient_labels[, .(patient, clinical_response, response)]), 
      all.x = T, all.y = F, by = 'patient')

## Read in RNA
rna_read_counts <- read_excel(file.path(data_dir, 'readcounts.xlsx'), sheet = 1)
setDT(rna_read_counts)
normalize_colnames(rna_read_counts)
## Select all annotation columns
rna_read_counts_ann <- rna_read_counts[, grep('4698', colnames(rna_read_counts),
                                          value = T, invert = T), with = F]
## TODO normalize read counts (edgeR)

rna_read_counts <- rna_read_counts[, lapply(.SD, sum), 
  # by = grep('4698', colnames(rna_read_counts), value = T, invert = T),
  # by = external_gene_id,
  by = ensembl_gene_id,
  .SDcols = grep('4698', colnames(rna_read_counts), value = T)]
# stopifnot(!rna_read_counts[, any(table(external_gene_id) > 1)])
# rna_read_counts <- merge(rna_read_counts, rna_read_counts_ann, 
#                          all.x = T, all.y = F)
setnames(rna_read_counts, 
         gsub('4698_\\d{1,2}_|_[tcga]{7}', '', colnames(rna_read_counts)))
rna_read_counts <- rna_read_counts %>% 
  column_to_rownames(var = 'ensembl_gene_id')

if (F) {
  which.duplicated <- function(vec) unique(vec[duplicated(vec)])
  rna_read_counts_ann %>%
    { .[external_gene_id %in% which.duplicated(external_gene_id)] } %>%
    { .[order(external_gene_id)] } %>%
    { DT::datatable(.) }
}

if (F) {
  ## TODO compute TPM instead
  library(edgeR)
  library(limma)
  d <- voom(rna_read_counts)
  # d <- DGEList(counts = rna_read_counts, genes = rownames(rna_read_counts))
  calcNormFactors(d)
  rna_read_counts_cpm <- cpm(d)
  rownames(rna_read_counts_cpm) <- rownames(rna_read_counts)
}

if (F) {
# patients = rna_sample_annotation[, unique(patient)]
# timepoints = 'Baseline'
# ann_subs <- rna_sample_annotation[patient %in% patients & 
#                                   timepoint %in% timepoints]
rna_read_counts <- rna_read_counts[, c('external_gene_id', 
                                rna_sample_annotation[, cf_number]), 
                            with = F] %>%
  melt(id.vars = 'external_gene_id', variable.name = 'cf_number') %>%
  merge(rna_sample_annotation)
rna_read_counts[, cf_number := NULL]
setkey(rna_read_counts, patient, external_gene_id, timepoint)
rna_read_counts[, t_number := NULL]
setcolorder(rna_read_counts, c('patient', 'induction_arm', 'timepoint',
                        'external_gene_id', 'value'))
rna_read_counts[, table(patient)]

# rna_FC <- rna_read_counts[, .('FC' = compute_FC(.SD)), 
#                    by = c('patient', 'external_gene_id')]
}
