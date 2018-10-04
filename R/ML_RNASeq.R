prepare_ml_dat_rnaseq <- function(tp = 'Baseline', type = 'genes',
                                  cpm_normalize = F,
                                  exp_mat = rna_read_counts_salmon_tmm_M) {
  type <- match.arg(type, choices = c('genes', 'gene_sets'), several.ok = F)

  t_dat_info <- rna_sample_annotation[timepoint %in% tp,
                      .(cf_number, patient, arm, clinical_response)] %>%
    { .[!is.na(clinical_response)] }

  ## As we're not using the variance stab, this is redundant and the formula
  ## does not need to be consistent with the actual experiment
  design_mat <- stats::model.matrix(as.formula(sprintf('~ 0 + clinical_response')), 
                                    as.data.frame(t_dat_info))
  idx <- match(t_dat_info[, cf_number], colnames(exp_mat))
  RCs <- exp_mat[, idx, with = F]
  rownames(RCs) <- rownames(exp_mat)
  # RCs$ensembl_gene_id <- rownames(exp_mat)
  # RCs$gene_symbol <-
  #   rna_read_counts_ann[match(RCs$ensembl_gene_id, ensembl_gene_id),
  #                       external_gene_id]
  # RCs$ensembl_gene_id <- NULL
  # RCs <- RCs[, lapply(.SD, sum), by = gene_symbol]
  # RCs <- as.data.frame(RCs) %>% column_to_rownames('gene_symbol')
  stopifnot(ncol(RCs) == nrow(t_dat_info))
  if (cpm_normalize) {
    cpms <- limma::voom(RCs, design_mat, plot = F, span = .1)
    cpms <- t(cpms$E)
    colnames(cpms) <- rownames(exp_mat)
    str(dimnames(cpms))
  } else {
    cpms <- as.matrix(t(RCs))
    colnames(cpms) <- rownames(exp_mat)
  }

  if (type == 'genes') {
  } else if (type == 'gene_sets') {
    browser()
  }

  ## Get rid of genes with NA values for some patients
  cpms <- cpms[, apply(cpms, 2, function(x) all(!is.na(x)))]
  ## Get rid of genes with small variance
  cpms <- cpms[, apply(cpms, 2, function(x) !eps(var(x), 0, 5e-1))]

  ## Look up responses
  setkey(t_dat_info, cf_number)
  responses <- t_dat_info[match(rownames(cpms), cf_number), clinical_response]
  return(list('t_dat' = cpms, 'responses' = responses))
}


glmnet_class_test_rnaseq <- function(alphas = seq(0.0, 1, by = .1),
                                     ncol = 4,
                                     nrow = 3,
                                     type = 'genes',
                                     tp = 'Baseline',
                                     nfolds = nrow(t_dat)) {
  prep <- prepare_ml_dat_rnaseq(tp = tp, type = type)
  t_dat <- prep[['t_dat']]
  responses <- prep[['responses']]

  print(table(responses))

  if (any(table(responses) == 0)) {
    stop('no responders among patients with this timepoint')
  }

  if (length(alphas) > (ncol * nrow)) {
    stop('adjust amount of plots to number of alphas')
  }

  indices <- 1:(ceiling(length(alphas) / ncol / nrow) * ncol * nrow)
  if (length(alphas) > 1) {
    layout(matrix(indices, ncol = ncol, nrow = nrow, byrow = T))
  }
  null_model <- table(responses) %>% { .[which.min(.)] / sum(.) }

  for (a in alphas) {
    res <- cv.glmnet(x = t_dat,
                     y = responses,
                     alpha = a,
                     family = 'binomial',
                     type.measure = 'class',
                     nfolds = nfolds)
    plot(res)
    abline(a = null_model, b = 0, type = 'dashed')
    title(parse(text=sprintf('alpha==%s', a)))
  }
  return(res)
}
