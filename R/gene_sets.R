gene_set_dir <- file.path('/DATA', 'resources', 'GSEA_genesets')

list_gene_sets <- function() list.files(gene_set_dir)

find_gene_set <- function(pat = 'c7') {
  gene.sets <- grep(pat, list.files(gene_set_dir, full.names = T), value = T)
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
