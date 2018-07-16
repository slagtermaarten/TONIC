sr(seq_stat_overview <- 
  rbindlist(lapply(patient_labels[, unique(patient)], compute_dnaseq_stats),
            fill = T))
# saveRDS(seq_stat_overview, 'rds/seq_stat_overview.rds')

rank_percentile <- function(dtf, rank_var, continuous = F) {
  ranked_var <- sprintf('%s_ranked', rank_var)
  if (continuous) {
    dtf[order(get(var)), (ranked_var) := (seq(1, .N) - 1) / (.N - 1)]
  } else {
    dtf[, (ranked_var) := (frank(get(rank_var), 
                                 ties.method = 'max') - 1) / (.N - 1)]
  }
  return(dtf)
}

for (vn in grep('ranked', setdiff(colnames(seq_stat_overview), 'patient'),
                value = T, invert = T)) {
  rank_percentile(seq_stat_overview, vn)
}
invisible(seq_stat_overview[, 'DD_rank' := mut_load_ranked + weighted_gen_SCNA_score])
# plot(seq_stat_overview[, .(mut_load_ranked, weighted_gen_SCNA_score)])
