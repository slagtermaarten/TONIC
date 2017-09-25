## Step 1: identify which ones are missing
analyse_missing_genes <- function(gene_sets = find_gene_set('c7')) {
  messagef('Parsing %s', gene_sets)
  query_genes <- ggsea::read_gmt(file.path(gene_set_dir, gene_sets)) %>%
                                 unlist %>%
                                 unique()
  ## Genes that were lost due to missing annotation
  mg_ann <- setdiff(query_genes, RC_f[, 'external_gene_id'])
  ## Genes that are missing in annotation or lost due to zero std
  mg_std <- setdiff(query_genes, colnames(X_mat_r))
  
  messagef(paste0('Missing %d/%d (%.2f) genes in gene set(s), ',
                  '(%.2f) no variance, ',
                  '(%.2f) missing annotation'),
           length(mg_std), length(query_genes),
           length(mg_std) / length(query_genes),
           # length(mg_std) - length(mg_ann), length(query_genes),
           (length(mg_std) - length(mg_ann)) / length(query_genes),
           # length(mg_ann), length(query_genes),
           (length(mg_ann)) / length(query_genes))
  return(mg_ann)
}
if (F) {
  missing_by_geneset <- sapply(list_gene_sets(), analyse_missing_genes)
  comb_missing <- sort(unique(unlist(missing_by_geneset)))
  cat(paste(comb_missing, collapse = ' '), file = 'dat/missing_genes_comb.txt')
  # queryMany(mg, scopes='symbol', fields=c('symbol', 'reporter'), 
  #           returnall = T, species='human')
}

## Step 2: Query [http://www.genenames.org/] and read in results
missing_symb_fh <- fread('dat/symbol_checker_missing_symb.tsv')
setnames(missing_symb_fh, colnames(missing_symb_fh), 
         tolower(gsub(" ", "_", colnames(missing_symb_fh))))
missing_symb_fh <- missing_symb_fh[, .(input, match_type, approved_symbol)]
recon_symb <- missing_symb_fh[approved_symbol != input & approved_symbol != "", 
                setNames(approved_symbol, input)] %>%
                { grep(pattern = 'withdrawn', x = ., value = T, invert = T) }
overlap_symb <- sort(setdiff(recon_symb, RC_f[, 'external_gene_id']))


#' Update gene sets with more up to date gene symbols
#'
#'
update_genesets <- function(querystr = 'c7') {
  gene_sets <- read_gmt(find_gene_set(querystr))
  return(lapply(gene_sets, function(gs) {
    saveable_genes <- intersect(gs, names(recon_symb))
    if (length(saveable_genes) > 0) {
      return(c(setdiff(gs, saveable_genes), recon_symb[saveable_genes]))
    } else {
      return(gs)
    }
  }))
}
