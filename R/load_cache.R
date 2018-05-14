# p_vals <- maartenutils::cond_readRDS(file.path('rds', 'p_vals.rds'))
# if (!is.null(p_vals)) {
#   sig_gene_sets <- lapply(p_vals[['gene_set']], 
#                           function(x) x[p_val <= .05, gene_set]) %>% 
#     unlist %>% unique %>% rev 
# }

