gene_set_dir <- file.path('~/libs', 'GSEA_genesets')

list_gene_sets <- function(gene_symb = 'entrez') {
  list.files(gene_set_dir, pattern = gene_symb)
}

find_gene_set <- function(pat = 'c7', gene_symb = 'symbols') {
  gene_symb <- match.arg(gene_symb, choices = c('entrez', 'symbols'), 
                         several.ok = F)
  pat_m <- sprintf('%s.*%s', pat, gene_symb)
  gene.sets <- grep(pat_m, list.files(gene_set_dir, full.names = T), value = T)
  if (length(gene.sets) == 0) {
    stop('No gene set found')
  }
  if (length(gene.sets) > 1) {
    warning('Found multiple gene sets: ', paste(gene.sets, collapse = ', '),
            '; returning first')
    gene.sets <- gene.sets[1]
  }
  return(gene.sets)
}
find_gene_set()
