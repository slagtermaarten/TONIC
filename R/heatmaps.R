#' Plot a heatmap from a geneset
#'
#' TODO show pval from GSEA in heatmap
gen_heatmap <- function(gene_set, pname='test',
                        dir_name='plots', pval=NA) {
  stopifnot(length(gene_set) > 1 & is.character(gene_set))
  setkey(fh, cell_line)
  annotation <- fh[rownames(X_mat), .(PDL1_0, PDL1_10, PDL1_FC)]
  aheatmap(log2(t(X_mat[, colnames(X_mat) %in% gene_set]) + 1),
           filename = file.path(dir_name, sprintf('heatmap_%s.pdf', pname)),
           main = gsub('_', ' ', pname),
           annCol = log2(annotation))
}


#' Plot a heatmap from a geneset
#'
#' @return nothing, side effect is creating files with plots in them
gen_multiple_heatmaps <- function(gene_sets, dir_name=NA,
                                  pvals, del_previous=T) {
  if (!is.na(dir_name)) {
    dir_name <- file.path('plots', dir_name)
    if (del_previous && dir.exists(dir_name)) {
      unlink(dir_name, recursive = T)
    }
    dir.create(dir_name, showWarnings = F)
  } else {
    dir_name <- file.path('plots')
  }
  pvals <- setNames(pvals, names(gene_sets))
  sapply(names(gene_sets), function(x) {
    gen_heatmap(gene_sets[[x]], pname = x, dir_name = dir_name, 
                pval = pvals[x])
  })
  invisible()
}
