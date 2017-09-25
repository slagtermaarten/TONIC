inspect_mat <- function(mat, 
                        maxR=min(5, nrow(mat)), 
                        maxC=min(5, ncol(mat))) 
{
  print(class(mat)); print(dim(mat)); print(mat[1:maxR, 1:maxC])
}

harmonic_mean <- function(v1 = c(1, 2, 3), v2 = c(3, 3, 4)) {
  stopifnot(length(v1) == length(v2))
  sapply(1:length(v1), function(i) 1 / mean(1 / c(v1[i], v2[i])))
}
# harmonic_mean()

