ggsea_wrapper <- function(gene.sets = find.gene_set('c7'),
                          y_var = resp_var) {
  ggsea_sm <- function (x, y, abs = F) {
      if (!is.matrix(y)) {
          y = matrix(y, dimnames = list(names(y), 'Response'))
      }
      y_r <- rank(y)
      n.response <- ncol(y)
      n.genes <- ncol(X_mat_r)
      coef <- t(lm(X_mat_r ~ y)$coefficients[-1, , drop = F])
      browser(expr = any(is.na(coef)))
      colnames(coef) <- colnames(y)
      rownames(coef) <- colnames(X_mat_r)
      coef <- apply(coef, 2, '/', X_mat_r_sd)
      if (abs) {
          coef <- base::abs(coef)
      }
      coef
  }

  ggsea_res <- ggsea(
      x = X_mat_r, 
      y = fh[, get(y_var)], 
      gene.sets = gene.sets,
      gene.score.fn = ggsea_sm, 
      es.fn = ggsea_weighted_ks,
      sig.fun = ggsea_calc_sig_simple, 
      # gene.names = NULL, 
      nperm = 1000,
      gs.size.min = 3, 
      # gs.size.max = 300, 
      verbose = TRUE, 
      block.size = 100,
      parallel = TRUE, 
      abs = TRUE
      # return_values = character()
  )

  ggsea_res[[1]][[1]] %>% 
    arrange(p) %>% 
    `[`(, c('es', 'fdr', 'p', 'GeneSet'))
}
