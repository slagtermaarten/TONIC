X_mat <- t(RC_f_n)
X_mat <- X_mat[, which(apply(X_mat, 2, sd) != 0)]
# inspect_mat(X_mat)
## Rank matrix in order to compute Spearman correlation
X_mat_r <- apply(X_mat, 2, rank)
# inspect_mat(X_mat_r)
X_mat_r_sd <- apply(X_mat_r, 2, sd)
stopifnot(all(X_mat_r_sd != 0))
stopifnot(!any(is.na(X_mat_r_sd)))
stopifnot(!any(is.na(X_mat)))
stopifnot(!any(is.na(fh[, get(resp_var)])))

if (F) {
  X_mat_e <- t(RC_f_n_e)
  X_mat_e <- X_mat_e[, which(apply(X_mat_e, 2, sd) != 0)]
# inspect_mat(X_mat)
## Rank matrix in order to compute Spearman correlation
  X_mat_r_e <- apply(X_mat_e, 2, rank)
# inspect_mat(X_mat_r)
  X_mat_r_sd_e <- apply(X_mat_r_e, 2, sd)
  stopifnot(all(X_mat_r_sd_e != 0))
  stopifnot(!any(is.na(X_mat_r_sd_e)))
  stopifnot(!any(is.na(X_mat_e)))
}

regress_vars <- list(
  'symbol' = list(X_mat_r = X_mat_r,
                  X_mat = X_mat,
                  X_mat_r_sd = X_mat_r_sd),
  'external_gene_id' = list(X_mat_r = X_mat_r,
                  X_mat = X_mat,
                  X_mat_r_sd = X_mat_r_sd)
  # , 'entrezgene' = list(X_mat_r = X_mat_r_e,
  #                     X_mat = X_mat_e,
  #                     X_mat_r_sd = X_mat_r_sd_e)
)
rm(list = grep("X_mat", ls(), value = T))
rm(list = grep("RC_f", ls(), value = T))
