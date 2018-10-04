BLOCK_SIZE <- 25

HALLMARK_pathways <- c(filter_gmt('h.all', 'HALLMARK'),
 filter_gmt('ifngamma', '*'),
# extra_pathways <- filter_gmt('custom.cancer.immunity', '*')
 suppressWarnings(filter_gmt('custom.savas', '*')),
 suppressWarnings(filter_gmt('custom.gutrantien', '*')))
# WNT_pathways <- grep('WNT', names(gmt), value = T)
# BRCA_pathways <- grep('BRCA', names(gmt), value = T)
# TP53_pathways <- grep('TP53', names(gmt), value = T)
# REACTOME_pathways <- grep('REACTOME', names(gmt), value = T)
# TARGET_pathways <- c('PAM50', grep('HALLMARK', names(gmt), value = T)) %>%
#   { grep('WNT|P53|PAM', ., value = T) }

HALLMARK_pathways <- filter_gmt('h.all', 'HALLMARK')

plot_bool <- interactive()

#' FC based ranking function
#'
#'
my_paired_FC_test <- function(x, y, abs = F) {
  ## Compute FCs per patient
  FCs <- rbindlist(lapply(unique(y$patient), function(l_patient) {
    pre_cf <- y[y$patient == l_patient &
                y$timepoint == 'Baseline', cf_number]
    post_cf <- y[y$patient == l_patient &
                 y$timepoint == 'Post-induction', cf_number]
    logFC <- unlist(x[match(post_cf, rownames(x)), ]) -
             unlist(x[match(pre_cf, rownames(x)), ])
    return(as.list(logFC))
  }), fill = T)

  ## Summarize genes by median FCs
  ret_val <- as.matrix(apply(FCs, 2, median))
  colnames(ret_val) <- 'FC'
  return(ret_val)
}


#' WC test ranking function
#'
#'
my_paired_WC_test <- function(x, y, abs = F) {
  cur_pats <- unique(y$patient)
  pre_idx <- vapply(cur_pats,
                    function(pat) which(y$patient == pat &
                                        y$timepoint == 'Baseline'),
                    integer(1))
  post_idx <- vapply(cur_pats,
                     function(pat) which(y$patient == pat &
                                         y$timepoint == 'Post-induction'),
                     integer(1))
  ret_val <- apply(x, 2, function(r) {
    w_test <- suppressWarnings(
              tryCatch(wilcox.test(x = r[pre_idx], y = r[post_idx],
                                   paired = T),
                       error = function(e) { print(e); browser() }))
    t_val <- (1 - w_test$p.value) *
      (log2(median(r[post_idx]) + 1) - log2(median(r[pre_idx]) + 1))
    return(t_val)
  })
  ret_val[is.na(ret_val)] <- 0
  # rownames(rna_read_counts_salmon_tmm_M)
  # ret_val[which(is.nan(ret_val))[1]]
  ret_val <- as.matrix(ret_val)
  colnames(ret_val) <- 't_val'
  if (abs) {
    ret_val <- abs(ret_val)
  }
  return(ret_val)
}


#' WC test ranking function
#'
#'
my_unpaired_WC_test <- function(x, y, abs = F) {
  R_idx <- y[, which(clinical_response == 'R')]
  NR_idx <- y[, which(clinical_response == 'NR')]

  ret_val <- apply(x, 2, function(r) {
    w_test <- tryCatch(wilcox.test(x = r[R_idx], y = r[NR_idx], paired = F, 
                                   exact = F),
                       error = function(e) { print(e); browser() })
    t_val <- sign(median(r[R_idx], na.rm = T) - median(r[NR_idx], na.rm = T)) *
              (1 - w_test$p.value)
    return(t_val)
  })
  # options(digits = 10)
  ret_val[is.na(ret_val)] <- 0
  ret_val <- as.matrix(ret_val)
  colnames(ret_val) <- 't_val'

  if (abs) {
    ret_val <- abs(ret_val)
  }
  return(ret_val)
}


#' WC test ranking function
#'
#'
my_unpaired_WC_test_ca15_3 <- function(x, y, abs = F) {
  HI_idx <- y[, which(baseline_ca15_3_bin == T)]
  LO_idx <- y[, which(baseline_ca15_3_bin == F)]

  ret_val <- apply(x, 2, function(r) {
    w_test <- tryCatch(wilcox.test(x = r[HI_idx], y = r[LO_idx], paired = F, 
                                   exact = F),
                       error = function(e) { print(e); browser() })
    t_val <- sign(median(r[HI_idx], na.rm = T) - median(r[LO_idx], na.rm = T)) *
              (1 - w_test$p.value)
    return(t_val)
  })
  ret_val[is.na(ret_val)] <- 0
  ret_val <- as.matrix(ret_val)
  colnames(ret_val) <- 't_val'

  if (abs) {
    ret_val <- abs(ret_val)
  }
  return(ret_val)
}


#' Run GSEA
#'
#'
gsea_wrapper <- function(gene_sets = HALLMARK_pathways,
                         patients = rna_sample_annotation[arm == 'Doxorubicin',
                                                          unique(patient)],
                         gene_score_fn = my_paired_WC_test,
                         paired_test = T,
                         preprocessing = 'none',
                         resp_exp = NULL,
                         exp_mat = rna_read_counts_salmon,
                         allowed_timepoints = NULL,
                         nperm = 25,
                         abs = F,
                         fn_extra = '') {
  starttime <- Sys.time()
  if (!exists('rna_sample_annotation')) source('R/load_rna_dat.R')

  ## Prepare X and Y matrices
  Y_mat <- rna_sample_annotation[patient %in% patients,
                                 .(patient, cf_number,
                                   timepoint, clinical_response)] %>%
    controlled_merge(patient_labels[timepoint == 'Baseline', 
                                    .(patient, 'baseline_ca15_3' = ca15_3)]) %>%
    mutate(baseline_ca15_3_bin = baseline_ca15_3 <= 
           patient_labels[!is.na(ca15_3) & clinical_response == 'R', 
                          max(ca15_3)])

  ## Subselect patients for which response var is not NA
  if (!is.null(resp_exp)) {
    Y_mat <- Y_mat[!is.na(get(resp_exp))]
    if (Y_mat[, class(get(resp_exp))] == 'factor') {
      Y_mat <- Y_mat[, (resp_exp) := droplevels(get(resp_exp))]
    }
  }

  allowed_timepoints <- allowed_timepoints %||% timepoints
  Y_mat <- Y_mat[timepoint %in% allowed_timepoints]
  if (null_dat(Y_mat)) return(NULL)

  if (paired_test && !null_dat(Y_mat)) {
    allowed_pats <- Y_mat[, .N, by = patient][N == 2, patient]
    if (length(allowed_pats) == 0) {
      return(NULL)
    }
    Y_mat <- Y_mat[patient %in% allowed_pats]
  }

  if (null_dat(Y_mat) || Y_mat[, any(table(get(resp_exp)) == 0)] ||
      nrow(Y_mat) < 6) {
    return(NULL)
  }

  idx <- match(Y_mat[, cf_number], colnames(exp_mat))
  if ('data.table' %in% class(exp_mat)) {
    RCs <- exp_mat[, idx, with = F]
    rownames(RCs) <- rownames(exp_mat)
  } else if ('matrix' %in% class(exp_mat)) {
    RCs <- exp_mat[, idx]
  } else {
    stop('Not implemented')
  }
  if (grepl('ENSG', rownames(exp_mat)[1])) {
    RCs$ensembl_gene_id <- rownames(exp_mat)
    RCs$gene_symbol <-
      rna_read_counts_ann[match(RCs$ensembl_gene_id, ensembl_gene_id),
                          external_gene_id]
    RCs$ensembl_gene_id <- NULL
    RCs <- RCs[, lapply(.SD, sum), by = gene_symbol]
    RCs <- as.data.frame(RCs) %>% column_to_rownames('gene_symbol')
  }

  ## Ensure we found data for all selected samples
  stopifnot(ncol(RCs) == nrow(Y_mat))
  if (preprocessing == 'cpm') {
    ## As we're not using the variance stab, this is redundant and the formula
    ## does not need to be consistent with the actual experiment
    design_mat <- stats::model.matrix(as.formula(sprintf('~ 0 + timepoint')),
                                      as.data.frame(Y_mat))
    RCs <- limma::voom(RCs, design_mat, plot = plot_bool, span = .1)
    RCs <- t(RCs$E)
    dimnames(RCs)
    colnames(RCs) <- rownames(exp_mat)
  } else if (preprocessing == 'none') {
    RCs <- t(DGEList(RCs)$counts)
    colnames(RCs) <- rownames(exp_mat)
  }
  browser()
  res <- ggsea(x = RCs,
               y = Y_mat,
               gene.sets = gene_sets,
               gene.score.fn = gene_score_fn,
               es.fn = ggsea_weighted_ks,
               sig.fun = ggsea_calc_sig,
               nperm = nperm,
               gs.size.min = 2,
               gs.size.max = 1e4,
               verbose = TRUE,
               block.size = BLOCK_SIZE,
               parallel = !local_run,
               return_values = c('leading_edge'),
               abs = abs)
  dtf <- as.data.frame(res[[1]])
  colnames(dtf) <- gsub('(.*)\\.', '', colnames(dtf))
  dtf$resp_exp <- resp_exp
  dtf$abs <- abs
  if (!is.null(res$leading_edge[[1]])) {
    dtf$leading_edge_genes <- sapply(res$leading_edge[[1]],
                                     function(x)
                                     paste(colnames(RCs)[x], collapse = ';'))
  }
  endtime <- Sys.time(); d <- endtime - starttime
  messagef(paste0('Finished. #samples: %d, response var: %s, ',
                  '#permutations: %d, gene ranks: %s, elapsed time: %.2f %s'),
           length(patients),
           resp_exp,
           nperm,
           ifelse(abs, 'absolute', 'relative'),
           as.list(d)[[1]],
           attr(d, 'units'))
  return(dtf)
}


gsea_all_arms <- function(gene_sets = HALLMARK_pathways,
                          patients = rna_sample_annotation[!is.na(clinical_response),
                                                           unique(patient)],
                          allowed_timepoints = c('Baseline', 'Post-induction'),
                          gene_score_fn = my_unpaired_WC_test,
                          exp_mat = rna_read_counts_salmon_tmm,
                          resp_exp = 'timepoint',
                          paired_test = F,
                          nperm = 1000,
                          abs = F,
                          fn_extra = '') {
  if (!exists('rna_sample_annotation')) source('R/load_rna_dat.R')
  arm_specific_gsea <-
    rbindlist(lapply(c(rna_sample_annotation[, levels(arm)], 'All arms'),
                     function(l_arm) {
      l_patients <- intersect(rna_sample_annotation[, patient], patients)
      if (l_arm != 'All arms') {
        l_patients <- intersect(l_patients,
                                rna_sample_annotation[arm %in% l_arm, patient])
      }
      res <- gsea_wrapper(gene_sets = gene_sets,
                          patients = l_patients,
                          allowed_timepoints = allowed_timepoints,
                          gene_score_fn = gene_score_fn,
                          exp_mat = exp_mat,
                          resp_exp = resp_exp,
                          paired_test = paired_test,
                          nperm = nperm,
                          abs = abs,
                          fn_extra = fn_extra)
      if (!is.null(res)) res$arm <- l_arm
      return(res)
    }), fill = T)

  if (!null_dat(arm_specific_gsea)) {
    arm_specific_gsea$arm <- factor(arm_specific_gsea$arm,
                                    levels = c('All arms', treatment_arms))
  }
  return(arm_specific_gsea)
}
