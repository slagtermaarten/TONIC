library(ggsea)
library(limma)
library(edgeR)
devtools::load_all('~/libs/GSEAgenesets')

BLOCK_SIZE <- 25

HALLMARK_pathways <- filter_gmt('h.all', 'HALLMARK')
# WNT_pathways <- grep('WNT', names(gmt), value = T)
# BRCA_pathways <- grep('BRCA', names(gmt), value = T)
# TP53_pathways <- grep('TP53', names(gmt), value = T)
# REACTOME_pathways <- grep('REACTOME', names(gmt), value = T)
# TARGET_pathways <- c('PAM50', grep('HALLMARK', names(gmt), value = T)) %>%
#   { grep('WNT|P53|PAM', ., value = T) }

plot_bool <- interactive()

#' FC based ranking function
#'
#'
my_paired_FC_test <- function(x, y, abs = F) {
  ## Compute FCs per patient
  FCs <- rbindlist(lapply(unique(y$patient), function(patient) {
    pre_cf <- y[y$patient == parent.frame(3)$patient &
                y$timepoint == 'Baseline', cf_number]
    post_cf <- y[y$patient == parent.frame(3)$patient &
                 y$timepoint == 'Post-induction', cf_number]
    logFC <- unlist(x[match(post_cf, rownames(x)), ]) - unlist(x[match(pre_cf, rownames(x)), ])
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
  ## TODO compute an exact WC here
  ret_val <- apply(x, 2, function(r) {
    w_test <- tryCatch(wilcox.test(x = r[pre_idx], y = r[post_idx], paired = T),
                       error = function(e) { print(e); browser() })
    t_val <- sign(median(r[post_idx]) - median(r[pre_idx])) *
      (1 - w_test$p.value)
    return(t_val)
  })
  ret_val <- as.matrix(ret_val)
  colnames(ret_val) <- 'p_val'

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

  # v <- x[, which(colnames(x) == 'FGR')]
  # median(v[R_idx])
  # median(v[NR_idx])

  ret_val <- apply(x, 2, function(r) {
    w_test <- tryCatch(wilcox.test(x = r[R_idx], y = r[NR_idx], paired = F),
                       error = function(e) { print(e); browser() })
    t_val <- sign(median(r[R_idx]) - median(r[NR_idx])) * (1 - w_test$p.value)
    return(t_val)
  })
  ret_val <- as.matrix(ret_val)
  colnames(ret_val) <- 'p_val'

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
                         resp_exp = NULL,
                         allowed_timepoints = NULL,
                         nperm = 25,
                         abs = F,
                         fn_extra = '') {
  starttime <- Sys.time()

  ## Prepare X and Y matrices
  resp_var <- rna_sample_annotation[patient %in% patients,
                                    .(patient, cf_number,
                                      timepoint, clinical_response)]
  ## Subselect patients for which response var is not NA
  if (!is.null(resp_exp)) {
    resp_var <- resp_var[!is.na(get(resp_exp))]
    if (resp_var[, class(get(resp_exp))] == 'factor') {
      resp_var <- resp_var[, (resp_exp) := droplevels(get(resp_exp))]
    }
  }

  if (is.null(allowed_timepoints)) {
    allowed_timepoints <- timepoints
  }
  resp_var <- resp_var[timepoint %in% allowed_timepoints]
  if (null_dat(resp_var)) {
    return(NULL)
  }

  if (paired_test && !null_dat(resp_var)) {
    allowed_pats <- resp_var[, .N, by = patient][N == 2, patient]
    if (length(allowed_pats) == 0) {
      return(NULL)
    }
    resp_var <- resp_var[patient %in% allowed_pats]
  }
  if (null_dat(resp_var)) {
    return(NULL)
  }

  ## As we're not using the variance stab, this is redundant and the formula
  ## does not need to be consistent with the actual experiment
  design_mat <- stats::model.matrix(as.formula(sprintf('~ 0 + timepoint')),
                                    as.data.frame(resp_var))
  idx <- match(resp_var[, cf_number], colnames(rna_read_counts))
  RCs <- rna_read_counts[, idx, with = F]
  RCs$ensembl_gene_id <- rownames(rna_read_counts)
  RCs$gene_symbol <-
    rna_read_counts_ann[match(RCs$ensembl_gene_id, ensembl_gene_id),
                        external_gene_id]
  RCs$ensembl_gene_id <- NULL
  RCs <- RCs[, lapply(.SD, sum), by = gene_symbol]
  RCs <- as.data.frame(RCs) %>% column_to_rownames('gene_symbol')

  ## Ensure we found data for all selected samples
  stopifnot(ncol(RCs) == nrow(resp_var))
  cpms <- limma::voom(RCs, design_mat, plot = plot_bool, span = .1)
  cpms <- t(cpms$E)

  res <- ggsea(x = cpms,
               y = resp_var,
               gene.sets = gene_sets,
               gene.score.fn = gene_score_fn,
               es.fn = ggsea_weighted_ks,
               sig.fun = ggsea_calc_sig,
               nperm = nperm,
               gs.size.min = 3,
               gs.size.max = 1e4,
               verbose = TRUE,
               block.size = BLOCK_SIZE,
               parallel = !local_run,
               # return_values = c('leading_edge'),
               abs = abs)
  dtf <- as.data.frame(res[[1]])
  colnames(dtf) <- gsub('(.*)\\.', '', colnames(dtf))
  dtf$resp_exp <- resp_exp
  dtf$abs <- abs
  if (!is.null(res$leading_edge[[1]])) {
    dtf$leading_edge_genes <- sapply(res$leading_edge[[1]],
                                     function(x)
                                     paste(colnames(cpms)[x], collapse = ';'))
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
                          resp_exp = 'timepoint',
                          paired_test = F,
                          nperm = 1000,
                          abs = F,
                          fn_extra = '') {
  arm_specific_gsea <-
    rbindlist(lapply(c(rna_sample_annotation[, levels(arm)], 'All arms'),
                     function(arm) {
      l_patients <- rna_sample_annotation[patient %in% parent.frame(3)$patients,
                                          patient]
      if (arm != 'All arms') {
        l_patients <- intersect(l_patients,
                                rna_sample_annotation[arm == parent.frame(3)$arm,
                                                      patient])
      }
      res <- gsea_wrapper(gene_sets = gene_sets,
                          patients = l_patients,
                          allowed_timepoints = allowed_timepoints,
                          gene_score_fn = gene_score_fn,
                          resp_exp = resp_exp,
                          paired_test = paired_test,
                          nperm = nperm,
                          abs = abs,
                          fn_extra = fn_extra)
      if (!is.null(res))
        res$arm <- arm
      return(res)
    }), fill = T)

  arm_specific_gsea$arm <- factor(res$arm, 
                                  levels = c('All arms', treatment_arms))
  return(arm_specific_gsea)
}
