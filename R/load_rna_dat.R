pacman::p_load(readxl)
pacman::p_load(naturalsort)

fn <- file.path(data_dir, 'Sample_annotations_RNAseq_TONICstageI.xlsx')
# sys_file_open(fn)
rna_sample_annotation <- read_excel(fn, sheet = 1, na = c("", "NA"))
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
rna_sample_annotation <-
  clean_columns(instance = 'load_rna_dat.R', rna_sample_annotation, 
                col_names = c('patient_id', 'induction_arm', 
                              'sequencing_number', 'study_id', 't-number'))
rna_sample_annotation <- rna_sample_annotation[patient != 'pat_NA']
setkey(rna_sample_annotation, patient)
setnames(rna_sample_annotation, 'timepoint_biopsy', 'timepoint')
rna_sample_annotation <-
  rna_sample_annotation[naturalsort::naturalorder(patient)]
rna_sample_annotation[, cf_number := tolower(cf_number)]
rna_sample_annotation <- controlled_merge(rna_sample_annotation,
      unique(patient_labels[, .(patient, arm, clinical_response, response)]),
      by = 'patient')

if (local_run) {
  fastqc_res <- parse_fastqc(fastq_dir = file.path(p_root, 
                                                   'fastq', 'rna', 'fastqc'))
  rna_sample_annotation <- 
    controlled_merge(rna_sample_annotation, fastqc_res, by = 'cf_number')

  gcf_rna_stats <- suppressWarnings(read_excel(file.path(data_dir,
                                                         'gcf_stats_4698.xlsx'), 
                                                 sheet = 1, na = c("", "NA"))) %>%
      as.data.table %>%
      normalize_colnames
  gcf_rna_stats <- gcf_rna_stats[, 1:18, with = F]
  gcf_rna_stats[, 'cf_number' := tolower(sample)]
  gcf_rna_stats[, cf_number := gsub('\\d{4}_\\d{1,2}_(cf.*)_\\w{7}', 
                                    '\\1', cf_number)]
  gcf_rna_stats <- clean_columns(instance = '', 
                                 fh = gcf_rna_stats, 
                                 col_names = c('x__1', 'sample'))
  setnames(gcf_rna_stats,
    c("bioanalyzer_nm_(mol/l)_(yield_after_libprep)",
    "nr_of_genes_with_normalized_expression_>_0",
    "nr_of_genes_with_normalized_expression_>_5",
    "nr_of_genes_with_normalized_expression_>_10",
    "nr_of_genes_with_normalized_expression_>_15"),
    c("bioanalyzer_nM",
    "nr_of_genes_with_normalized_expression_>_0",
    "nr_of_genes_with_normalized_expression_>_5",
    "nr_of_genes_with_normalized_expression_>_10",
    "nr_of_genes_with_normalized_expression_>_15"))

  rna_sample_annotation <- 
    controlled_merge(rna_sample_annotation, gcf_rna_stats, by = 'cf_number')

  rna_sample_annotation <- 
    controlled_merge(rna_sample_annotation, 
                     readRDS(file.path('rds', 'fastq_files.rds')))
}



if (F) {
  ## Read in RNA
  rna_read_counts <- suppressWarnings(read_excel(file.path(data_dir, 
                                                           'readcounts.xlsx'), 
                                                 sheet = 1, na = c("", "NA")))
  setDT(rna_read_counts)
  normalize_colnames(rna_read_counts)

  ## Select all annotation columns
  rna_read_counts_ann <-
    rna_read_counts[, grep('4698', colnames(rna_read_counts),
                           value = T, invert = T), with = F]
  rna_read_counts_ann[, 'gene_length' := end_position - start_position]

  rna_read_counts <- rna_read_counts[, lapply(.SD, sum),
    # by = grep('4698', colnames(rna_read_counts), value = T, invert = T),
    # by = external_gene_id,
    by = ensembl_gene_id,
    .SDcols = grep('4698', colnames(rna_read_counts), value = T)]

  ensgs <- rna_read_counts[, ensembl_gene_id]

  ## Select all sample columns, remove ensembl_id
  rna_read_counts <-
    rna_read_counts[, grep('4698', colnames(rna_read_counts),
                           value = T, invert = F), with = F]
  ## TODO normalize read counts (edgeR)
  # library(limma)
  # library(edgeR)
  # d <- DGEList(counts=t(dtf), group = group)
  # d <- calcNormFactors(d)
  # scaledMatrix <- limma::voom(d)

  ## Row names are not displayed in data.table objects, but there maintained
  ## nonetheless
  rownames(rna_read_counts) <- ensgs
  setnames(rna_read_counts,
           gsub('4698_\\d{1,2}_|_[tcga]{7}', '', colnames(rna_read_counts)))
  gene_symbols <- rna_read_counts_ann$external_gene_id %>%
    { .[match(ensgs, rna_read_counts_ann$ensembl_gene_id)] }
  rna_read_counts$gene_symbol <- gene_symbols
  rna_read_counts <- rna_read_counts[, lapply(.SD, sum), by = gene_symbol]
  rownames(rna_read_counts) <- rna_read_counts$gene_symbol
  rna_read_counts$gene_symbol <- NULL
}

if (F) {
  ## TODO 2018-07-05 12:55 Fix me when needed
  rna_read_counts$ensembl_gene_id <- rownames(rna_read_counts)
  tpms <- controlled_merge(rna_read_counts,
                           rna_read_counts_ann[, .(ensembl_gene_id,
                                                   gene_length)])
  rna_read_counts$ensembl_gene_id <- NULL
  ## Check first gene, second sample as representative for all genes
  # tpm_val <- tpms[1,cf10652] / tpms[1, gene_length]
  tpms <- tpms[, lapply(.SD, function(x) x / gene_length),
    .SDcols = colnames(rna_read_counts)]
  # stopifnot(tpms[1, cf10652] == tpm_val)

  if (F) {
    ## Checking the syntax I'm using
    dat <- matrix(1:9, nrow = 3, byrow = T)
    dat <- t(t(dat) / apply(dat, 2, sum))
    v <- c(3, 6, 9)
    v <- c(1, 4, 7)
    v / sum(v)
    dat <- t(t(dat) * 1e6 / apply(dat, 2, sum))
  }

  ## Normalise by library size
  tpms <- as.data.table(t(t(tpms) * 1e6 / apply(tpms, 2, sum)))
  tpms$hugo_symbol <- rna_read_counts_ann$external_gene_id[match(ensgs, rna_read_counts_ann$ensembl_gene_id)]
  tpms <- tpms[, lapply(.SD, sum), by = hugo_symbol]
  tpms <- column_to_rownames(tpms, 'hugo_symbol')
  # stopifnot(!rna_read_counts[, any(table(external_gene_id) > 1)])
  # rna_read_counts <- merge(rna_read_counts, rna_read_counts_ann,
  #                          all.x = T, all.y = F)
  # tpms <- as.data.frame(tpms)
}

if (T) {
  tpms_salmon <- fread(file.path(p_root, 'salmon_rna', 
                                 'salmon_tpm_mat.tsv')) %>%
    column_to_rownames('hugo_symbol') %>%
    normalize_colnames

  rna_read_counts_salmon <- fread(file.path(p_root, 'salmon_rna', 
                                  'salmon_count_mat.tsv')) %>%
    normalize_colnames
  rownames(rna_read_counts_salmon) <- rna_read_counts_salmon$hugo_symbol
  rna_read_counts_salmon$hugo_symbol <- NULL
}

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
