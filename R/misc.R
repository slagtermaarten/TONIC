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


#' Parse FastQC results
#'
#'
parse_fastqc <- function(fastq_dir = sprintf('%s/salmon_rna/fastqc', p_root)) {
  ## Merge FASTQC results
  fastqc_res <- 
    fread(sprintf('cat %s/CF*/**/summary.txt', fastq_dir), header = F)
  setnames(fastqc_res, c('fastqc_test_result', 'fastqc_test', 'filename'))
  fastqc_res[, 'cf_number' := 
             gsub('\\d{4}_\\d{1,2}_(CF.*)_\\w{7}_S\\d{1,2}_R\\d_\\d{3}\\.fastq\\.gz', '\\1', filename)]
  fastqc_res[, fastqc_test_result := fastqc_test_result == 'PASS']
  fastqc_res[, cf_number := tolower(cf_number)]
  fastqc_res[, filename := NULL]
  ## All samples must have same number of tests
  stopifnot(fastqc_res[, .N, by = cf_number][, uniqueN(N) == 1])
  ## Change to lowercase
  fastqc_res[, fastqc_test := variabilize_character(fastqc_test)]
  ## Widen data 
  fastqc_res <- dcast(fastqc_res, cf_number ~ fastqc_test, 
                      value.var = 'fastqc_test_result')
  return(fastqc_res)
}


#' Exclude some genes
#'
#'
filter_exp_mat <- function(dtf = tpms_salmon, 
                           search_term = 'ribosomal|mitochondrial',
                           search_var = 'hgnc_symbol',
                           invert = F) {
  entrez_table <- as.data.table(readRDS('rds/entrez_table.rds'))
  
  exclude_genes <- 
    entrez_table[grepl(search_term, get(search_var)), hgnc_symbol]
  if (invert) {
    exclude_genes <- setdiff(rownames(dtf), exclude_genes)
  }
  entrez_table[hgnc_symbol %in% exclude_genes, gene_biotype]
  messagef('Excluding %d/%d genomic features', length(exclude_genes),
           length(rownames(entrez_table)))

  row_names <- rownames(dtf) %>% { .[. %nin% exclude_genes] }
  dtf <- dtf[rownames(dtf) %nin% exclude_genes]
  rownames(dtf) <- row_names
  return(dtf)
}


#' Exclude some genes and subsequently renormalize TPMs
#'
#'
renormalize_tpms <- function(dtf = tpms_salmon, 
                             search_term = 'ribosomal|mitochondrial',
                             search_var = 'hgnc_symbol',
                             invert = F) {
  row_names <- rownames(dtf)
  dtf <- filter_exp_mat(dtf)
  ## Renormalize TPMs
  dtf <- dtf[, lapply(.SD, function(x) x * 1e6 / sum(x))]
  rownames(dtf) <- row_names
  stopifnot(all(eps(colSums(dtf), 1e6)))
  return(dtf)
}


#' Quantile normalisation of data.frame or data.table
#'
#' 
quantile_normalisation <- function(df) {
  df_rank <- apply(df, 2, rank, ties.method = 'min')
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
   
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
   
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}


#' Compute housekeeping gene quantiles
#'
#'
compute_hkg_quantiles <- function(dtf = tpms_salmon) {
  ## Read in housekeeping genes
  hkg <- fread(file.path(data_dir, 'HK_genes.txt'), skip = 1, header = F)
  setnames(hkg, c('hgnc_symbol', 'refseq'))

  hkg_dtf <- dtf[rownames(dtf) %in% hkg[, hgnc_symbol]]
  # per_sample_quants <- 
  #   hkg_dtf[, lapply(.SD, function(x) quantile(x, probs = c(.1, .5, .9)))]
  apply(hkg_dtf, 1, function(x) sd(x) / mean(x))
}
