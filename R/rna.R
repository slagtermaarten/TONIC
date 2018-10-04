prep_gene_set <- function(gs = NULL,
                          gene_symbols = c('TP53', 'ERBB2'),
                          gs_name = paste(gene_symbols, collapse = '_'),
                          cf_numbers = c(),
                          cf_ordered = F,
                          patients = c(),
                          induction_arms = c(),
                          timepoints = c('Baseline'),
                          cov_thresh = 0,
                          log_transform = 'log2',
                          exp_mat = tpms_salmon,
                          gene_normalisation = T,
                          sample_normalisation = NULL,
                          distfun = 'pearson',
                          hclustfun = 'average',
                          shuffle = F,
                          fn_extra = '') {
  if (class(gene_symbols) == 'list') {
    if (is.null(gs_name)) {
      gs_name <- names(gene_symbols)
    }
    gene_symbols <- gene_symbols[[1]]
  }

  ## Select CF numbers to plot
  if (!cf_ordered) {
    cf_numbers <- Reduce(union,
      c(cf_numbers,
        rna_sample_annotation$cf_number[match(patients,
                          rna_sample_annotation$patient)],
        rna_sample_annotation$cf_number[match(induction_arms,
                          rna_sample_annotation$induction_arm)],
        rna_sample_annotation$cf_number[as.character(rna_sample_annotation$timepoint) %in% timepoints]))

    if (length(cf_numbers) < 2) {
      warningf('Too little samples selected: %d. Timepoint: %s. Arms: %s', 
               length(cf_numbers), 
               ifelse(is.null(timepoints), 'unspecified', timepoints),
               ifelse(is.null(induction_arms), 'unspecified', induction_arms))
      return(invisible())
    }
  }

  # ensgs <- rna_read_counts_ann[which(external_gene_id %in% gene_symbols),
  #                              ensembl_gene_id]
  mean(gene_symbols %in% rownames(exp_mat))
  p_dat <-
    GSEAgenesets::preprocess_rna(exp_mat, 
                                 gs_genes = gene_symbols,
                                 log_transform = log_transform,
                                 gene_normalisation = gene_normalisation,
                                 sample_normalisation = sample_normalisation,
                                 cov_thresh = NULL, symbol_var = NULL,
                                 donor_ids = cf_numbers)
  browser(expr = is.null(p_dat) || nrow(p_dat) <= 1)
  if (F) {
    ## Convert ENSGs to gene symbols
    gene_idx <- rna_read_counts_ann[, which(ensembl_gene_id %in% rownames(p_dat))]
    rownames(p_dat) <- rna_read_counts_ann[gene_idx, external_gene_id]
    rownames(p_dat) <- NULL
    messagef('Gene set: %s. Omitting these genes due to conversion probs: %s',
             gs_name,
             paste(setdiff(gene_symbols, rownames(p_dat)), collapse = ', '))
    ## Get rid of duplicate gene symbols by averaging entries
    p_dat <- as.data.table(cbind(as.data.frame(p_dat), 
          'gene_symbol' = rna_read_counts_ann[gene_idx, external_gene_id]))
    p_dat <- p_dat[, lapply(.SD, mean), by = gene_symbol]
    row.names <- p_dat$gene_symbol
    p_dat <- as.data.frame(p_dat[, -1])
    rownames(p_dat) <- row.names
  }
  ## Ensure no duplicated symbols
  browser(expr = any(duplicated(rownames(p_dat))))
  return(p_dat)
}


plot_gene_set <- function(gs = NULL,
                          gene_symbols = c('TP53', 'ERBB2'),
                          gs_name = paste(gene_symbols, collapse = '_'),
                          cf_numbers = c(),
                          cf_ordered = F,
                          patients = c(),
                          induction_arms = c(),
                          gene_normalisation = T,
                          sample_normalisation = NULL,
                          timepoints = c('Baseline'),
                          col_vars = c('timepoint', 'arm',
                                       'clinical_response', 'gs_score'),
                          cov_thresh = NULL,
                          ann_row = NA,
                          exp_mat = rna_read_counts_salmon_tmm_M,
                          log_transform = 'log2',
                          distfun = 'euclidean',
                          distfun_horizontal = NULL,
                          distfun_vertical = NULL,
                          hclustfun = 'complete',
                          shuffle = F,
                          fn_extra = '') {
  if (class(gene_symbols) == 'list') {
    if (is.null(gs_name)) {
      gs_name <- names(gene_symbols)
    }
    gene_symbols <- gene_symbols[[1]]
  }
  exp_object_name <-deparse(substitute(exp_mat))

  distfun_horizontal <- distfun_horizontal %||% distfun
  distfun_vertical <- distfun_vertical %||% distfun

  if (!cf_ordered) {
    ## Select CF numbers to plot
    cf_numbers <- Reduce(union,
      c(cf_numbers,
        rna_sample_annotation$cf_number[match(patients,
                          rna_sample_annotation$patient)],
        rna_sample_annotation$cf_number[match(induction_arms,
                          rna_sample_annotation$induction_arm)],
        rna_sample_annotation$cf_number[as.character(rna_sample_annotation$timepoint) %in% timepoints]))

    if (length(cf_numbers) < 2) {
      warningf('Too little samples selected: %d', length(cf_numbers))
      return(invisible())
    }
  }

  p_dat <- prep_gene_set(gs = gs, 
                         gene_symbols = gene_symbols,
                         exp_mat = exp_mat,
                         gs_name = gs_name, 
                         cf_numbers = cf_numbers,
                         cf_ordered = cf_ordered,
                         patients = patients, 
                         induction_arms = induction_arms,
                         gene_normalisation = gene_normalisation,
                         sample_normalisation = sample_normalisation,
                         timepoints = timepoints, 
                         cov_thresh = cov_thresh,
                         log_transform = log_transform,
                         distfun = distfun, 
                         hclustfun = hclustfun,
                         shuffle = shuffle, 
                         fn_extra = fn_extra)

  if (null_dat(as.data.frame(p_dat))) return(NULL)

  # ann_col <-
  #   cbind(rna_sample_annotation[match(cf_numbers,
  #                                     rna_sample_annotation[, cf_number]),
  #                               .(cf_number, timepoint, arm, clinical_response)],
  #         patient_labels[match(cf_numbers, rna_sample_annotation[, cf_number]),
  #                        .(efron_thisted_estimator, sample_clonality,
  #                          adaptive_t_cells)])
  ann_col <- data.table(cf_number = cf_numbers) %>%
    controlled_merge(rna_sample_annotation[, .(cf_number, patient, arm,
                                               timepoint)]) %>%
    controlled_merge(patient_labels[, .(patient, timepoint,
                                        clinical_response,
                                        efron_thisted_estimator, 
                                        sample_clonality,
                                        adaptive_t_cells)])
  setkey(ann_col, cf_number)
  ann_col <- ann_col[cf_numbers]
  stopifnot(all(ann_col$cf_number == cf_numbers))

  if ('gs_score' %in% col_vars) {
    gs <- compute_gene_set_score(timepoint = timepoints, 
                                 gene_symbols = gene_symbols,
                                 log_transform = log_transform, 
                                 sum_func = median)
    ann_col <- controlled_merge(ann_col, gs[, .(cf_number, gs_score)])
  }
  ann_col <- ann_col[, col_vars, with = F]

  fn <- file.path(img_dir,
                  sprintf('heatmap_%s_%s%s%s%s%s%s%s.pdf',
                    sprintf('_%s', gs_name),
                    exp_object_name,
                    ifelse(is.null(timepoints) || is.na(timepoints),
                           '', sprintf('_%s',
                                       paste(timepoints, collapse = '-'))),
                    ifelse(is.null(log_transform) || is.na(log_transform),
                           '', sprintf('_%s', log_transform)),
                    sprintf('_%s', hclustfun),
                    sprintf('_%s', distfun),
                    ifelse(length(fn_extra) > 1, sprintf('_%s', fn_extra), ''),
                    ifelse(shuffle, '-shuf', '')
                  ))
  if (!all(is.na(ann_row))) {
    ## Ensure p_dat and ann_row are in same order
    i_names <- intersect(rownames(ann_row), rownames(p_dat))
    ann_row <- ann_row[match(i_names, rownames(p_dat)), ]
    if (ncol(ann_row) == 1) {
      ann_row <- ann_row[match(i_names, rownames(ann_row)), ]
      ann_row <- data.frame('subtype' = ann_row)
      rownames(ann_row) <- i_names
    } else {
      ann_row <- ann_row[match(i_names, rownames(ann_row)), ]
    }
  }

  if (cf_ordered) {
    stopifnot(all(cf_numbers == colnames(p_dat)))
    # cf_numbers == arrange(patient_labels_tmp, adaptive_t_cells)$rna_cf_number
  }
  # p_dat <- p_dat[rownames(p_dat) != 'NA', ]
  # p_dat <- p_dat[!apply(p_dat, 1, function(x) any(is.na(x))), ]

  if (F) {
    pheatmap::pheatmap(p_dat,
                       clustering_distance_rows = distfun_horizontal,
                       clustering_distance_cols = distfun_vertical,
                       # clustering_method = hclustfun,
                       # annotation_row = ann_row,
                       # annotation_col = ann_col,
                       # annotation_colors = tonic_color_palettes,
                       filename = fn,
                       main = tonic_cap(paste0(tolower(gs_name), fn_extra)))
  } else {
    NMF::aheatmap(p_dat,
                  labRow = NULL,
                  labCol = NULL,
                  Colv = ifelse(cf_ordered, NA, T),
                  distfun = distfun,
                  hclustfun = hclustfun,
                  annRow = ann_row,
                  annCol = ann_col,
                  # annColors = ann_colors,
                  filename = fn,
                  main = tonic_cap(paste0(tolower(gs_name), fn_extra)))
  }
  sys_file_open(fn)
}

plot_gene_sets <- function(sets = filter_gmt(gmt_pat = 'h.all'),
                           timepoints = c('Baseline', 'Post-induction',
                                          'On nivo')) {
  for (idx in seq_along(sets)) {
    for (timepoint in timepoints) {
      plot_gene_set(gs_name = NULL,
                    gene_symbols = sets[idx],
                    timepoints = timepoint)
    }
  }
}


compute_gene_set_score <- function(timepoint = 'Baseline',
                                   gene_symbols = filter_gmt(gmt_pat =
                                                             'h.all')[1],
                                   log_transform = NULL,
                                   exp_mat = rna_read_counts_salmon_tmm_M,
                                   sum_func = median) {
  t_dat <- prep_gene_set(gene_symbols = gene_symbols,
                         gene_normalisation = F,
                         exp_mat = exp_mat,
                         log_transform = log_transform,
                         timepoints = timepoint)
  if (null_dat(t_dat)) return(NULL)
  ret_val <- apply(t_dat, 2, sum_func)
  cf_numbers <- names(ret_val)
  ret_val <- data.table('gs_score' = ret_val, 'cf_number' = cf_numbers)
  ret_val <-
    controlled_merge(ret_val, 
                     rna_sample_annotation[, .(cf_number, patient, timepoint)])
  return(ret_val)
}


gen_gene_set_score_matrix <- function(sets = filter_gmt(gmt_pat = 'h.all'),
                                      sum_func = median,
                                      log_transform = NULL,
                                      exp_mat = rna_read_counts_salmon_tmm_M,
                                      timepoints = c('Baseline', 
                                                     'Post-induction',
                                                     'On nivo')) {
  dtf <- rbindlist(lapply(seq_along(sets), function(idx) {
    dtf <- rbindlist(lapply(timepoints, function(timepoint) {
      compute_gene_set_score(gene_symbols = sets[idx], 
                             exp_mat = exp_mat,
                             log_transform = log_transform,
                             timepoint = timepoint,
                             sum_func = sum_func)
    }), fill = T)
    if (null_dat(dtf)) return(NULL)
    dtf[, 'gene_set' := names(sets[idx])]
    return(dtf)
  }), fill = T)
  if (null_dat(dtf)) return(NULL)
  dtf_m <- melt(dtf, id.vars = c('timepoint', 'gene_set', 'cf_number', 
                                 'patient'),
                measure.vars = 'gs_score') 
  stopifnot(dtf_m[, uniqueN(variable)] == 1)
  dtf_m[, variable := NULL]
  dtf_m[, timepoint := factor(timepoint, levels = timepoints)]
  dtf_m <- dtf_m[naturalsort::naturalorder(patient)]
  dtf_m <- controlled_merge(dtf_m, 
                            patient_labels[, .(patient, 
                                               clinical_response, arm)])
  dtf_m[is.na(clinical_response), clinical_response := 'NA']
  dtf_m[, timepoint := factor(timepoint, levels = timepoints)]
  return(dtf_m)
}


#' Leading edge gene scores
#'
#' A bit of wishful thinking won't hurt (selecting leading edge genes, we expect
#' to see a signal but this signal is not unbiased anymore)
compute_leading_edge_gene_scores <- function(res, fdr_thresh = .25) {
  dtf <- rbindlist(lapply(res[, unique(arm)], function(l_arm) {
    leading_edge <- res[arm == l_arm & fdr <= fdr_thresh,
                        setNames(leading_edge_genes, GeneSet)]
    if (is.null(leading_edge) || length(leading_edge) == 0)
      return(NULL)
    gs <- strsplit(x = leading_edge, split = '\\;')
    if (is.null(gs) || length(gs) == 0)
      return(NULL)

    ## Iterate over gene sets
    dtf <- rbindlist(lapply(seq_along(gs), function(idx) {
      gdtf <- rbind(
        compute_gene_set_score(timepoint = 'Baseline',
                               gene_symbols = gs[[idx]],
                               log_transform = NULL, sum_func = median),
        compute_gene_set_score(timepoint = 'Post-induction',
                               gene_symbols = gs[[idx]],
                               log_transform = NULL, sum_func = median))
      gdtf[, 'gene_set' := names(gs)[idx]]
      return(gdtf)
    }))
    dtf[, 'gsea_arm' := l_arm]
    return(dtf)
  }), fill = T) %>%
  controlled_merge(patient_labels[, .(patient, arm, clinical_response)],
                   by = 'patient')
  dtf[, 'value' := gs_score]
  return(dtf)
}


#' Plot parallel coordinates of leading edge gene sets
#'
#' @param allow_nonsensical Only plot scores for gene set-arm combinations that
#' passed the FDR threshold
plot_leading_edge_gene_scores <- function(dtf, res_name = 'paired',
                                          allow_nonsensical = F) {

  if (!allow_nonsensical) {
    dtf <- dtf[as.character(gsea_arm) == 'All arms' |
               as.character(gsea_arm) == as.character(arm)]
  }
  plots <- dtf %>% 
    unique(by = c('gsea_arm', 'gene_set')) %>%
    dplyr::mutate(gsea_arm = as.character(gsea_arm)) %>%
    pmap(function(gene_set, gsea_arm, ...) {
      gs <- gene_set
      l_arm <- gsea_arm
      p_dat <- dtf[as.character(gene_set) == gs &
                   as.character(gsea_arm) == l_arm]
      # if (!allow_nonsensical) {
        if (l_arm == 'All arms') {
          facet_var <- NULL
        } else {
          facet_var <- 'arm'
        }
      # } else {
      #   facet_var <- 'arm'
      # }
      if (grepl('(?<!un)paired', res_name, perl = T)) {
        allowed_pats <- p_dat[, .N == 2, by = patient][V1 == T, patient]
        p_dat <- p_dat[patient %in% allowed_pats]
      }
      emph_labels <- rep('plain', 5)
      emph_labels[which(p_dat[, levels(arm)] == l_arm)] <- 'bold'
      plot_parallel_coords(p_dat = p_dat,
                           colour_var = 'clinical_response',
                           facet_var = facet_var) +
        ggtitle(sprintf('%s - %s', tonic_cap(gene_set), gsea_arm)) +
        ## Emphasize the arm in which the leading edge genes are defined
        ## Not working ATM
        theme(strip.text = element_text(face = emph_labels))
    })

  plot_panel_layout(plots, 
                    filename = sprintf('plots2/leading_edge_gene_scores_%s.pdf',
                                       res_name),
                    ncol = 1, nrow = 3, labels = NULL)
}


#' Wrapper function to limma functionality
#'
#'
my_tt <- function(e_fit, coef = ncol(model_mat), ...) {
  tt <- limma::topTable(e_fit, coef = coef, ...)
  if (F && any(grepl('ENSG', rownames(tt)[1:5]))) {
    tt$gene_symbol <-
      rna_read_counts_ann[match(rownames(tt),
                                rna_read_counts_ann[, ensembl_gene_id]),
                          external_gene_id]
  } else {
    tt$gene_symbol <- rownames(tt)
  }
  return(tt)
}


#' Venn diagram wrapper
#'
#'
my_venn <- function(results, ...) {
  vennDiagram(results, cex = .9, ...)
}
# my_tt(e_fit)


compute_evenness <- function(x) {
  stopifnot(class(x) == 'numeric' || class(x) == 'integer')
  x_norm <- x / sum(x, na.rm = T)
  ## Divide observed entropy by maximally obtainable entropy
  -sum(log((x_norm)^(x_norm))) / log(length(x))
}


#' Wrapper around prcomp
#'
#'
run_pca <- function(dtf = rna_read_counts_salmon_tmm_M_ensg, 
                    scale = T, center = T, filter_covariates = T) {
  t_dtf <- as.matrix(dtf) %>% t
  colnames(t_dtf) <- rownames(dtf)
  if (filter_covariates) {
    t_dtf %<>%  { .[, !eps(apply(., 2, var), 0, 1e-1)] }
    t_dtf %<>%  { .[, apply(., 2, mean) > 1] }
    t_dtf %<>%  { .[, apply(., 2, function(x) sd(x) / mean(x)) > 1] }
  }
  # dimnames(t_rna_read_counts_salmon)
  pca <- prcomp(t_dtf, center = center, scale = scale)
  return(pca)
}
