ggsea_wrapper <- function(gene.sets = find.gene_set('c7'),
                          resp_var = 'PDL1_FC', gene_var = 'entrez',
                          nperm = 1000, abs = F) {
  ggsea_sm <- function(x, y, abs = F) {
      if (!is.matrix(y)) {
          y = matrix(y, dimnames = list(names(y), 'Response'))
      }
      y_r <- rank(y)
      n.response <- ncol(y)
      n.genes <- ncol(regress_vars[[gene_var]][['X_mat_r']])
      coef <- t(lm(regress_vars[[gene_var]][['X_mat_r']] ~ y)$coefficients[-1, , drop = F])
      browser(expr = any(is.na(coef)))
      colnames(coef) <- colnames(y)
      rownames(coef) <- colnames(regress_vars[[gene_var]][['X_mat_r']])
      coef <- apply(coef, 2, '/', regress_vars[[gene_var]][['X_mat_r_sd']])
      if (abs) {
          coef <- base::abs(coef)
      }
      coef
  }

  ggsea_res <- ggsea(
      x = regress_vars[[gene_var]][['X_mat_r']], 
      y = fh[, get(resp_var)], 
      gene.sets = gene.sets,
      gene.score.fn = ggsea_sm, 
      es.fn = ggsea_weighted_ks,
      sig.fun = ggsea_calc_sig_simple, 
      # gene.names = NULL, 
      nperm = nperm,
      gs.size.min = 3, 
      # gs.size.max = 300, 
      verbose = TRUE, 
      block.size = 100,
      parallel = TRUE, 
      abs = abs
      # return_values = character()
  )

  ggsea_res[[1]][[1]] %>% 
    arrange(p) %>% 
    `[`(, c('es', 'fdr', 'p', 'GeneSet'))
}
