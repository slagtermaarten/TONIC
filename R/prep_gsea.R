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

