ton_percentile_comparison_of_fc <- function(N_reps = 1e4, 
                                            gene_sets = c('pd1', 'tis'),
                                            tp1 = 'Baseline',
                                            tp2 = 'On nivo') {
  p_dat <- prepare_test_gene_set_difference(gene_set = gene_sets, 
                                            tp1 = tp1,
                                            tp2 = tp2)
  p_dat <- p_dat[, .('logFC' = log2(get(tp2)) - log2(get(tp1)), arm),
                 by = patient]
  obs <- p_dat[, .('FC' = median(logFC)), by = arm]

  null_dist <- vapply(1:N_reps, function(x) {
    p_dat_perm <- p_dat[, .(patient, logFC, 'arm' = sample(arm))]
    obs_perm <- p_dat_perm[, median(logFC), by = arm]
    setkey(obs_perm, arm)
    obs_perm <- obs_perm[obs[, as.character(arm)]]
    return(obs_perm[, V1])
  }, numeric(5))

  compare_obs_to_null <- function(obs_val, null_dist) {
    # if (obs_val > 0) {
    #   return(1 - mean(null_dist >= obs_val) * 2)
    # } else {
    #   return(1 - mean(null_dist <= obs_val) * 2)
    # }
    if (obs_val > 0) {
      return(mean(null_dist >= obs_val))
    } else {
      return(mean(null_dist <= obs_val))
    }
    # return(mean(null_dist >= obs_val))
  }

  percentiles <- vapply(seq_along(obs[, FC]), function(idx) {
    obs_val <- obs[, FC][idx]
    l_null_dist <- null_dist[idx, ]
    return(compare_obs_to_null(obs_val = obs_val, null_dist = l_null_dist))
  }, numeric(1))
  obs$percentile <- percentiles

  hists <- lapply(seq_along(obs[, FC]), function(idx) {
    arm_l <- obs[idx, arm]
    obs_val <- obs[, FC[idx]]
    dtf <- data.frame(null_fc = null_dist[idx, ])
    # if (obs_val > 0) {
    #   dtf$col_val <- dtf$null_fc >= obs_val
    # } else {
    #   dtf$col_val <- dtf$null_fc <= obs_val
    # }
    if (obs_val > 0) {
      dtf$col_val <- factor(ifelse(dtf$null_fc <= obs_val, 'S', 'G'))
    } else {
      dtf$col_val <- factor(ifelse(dtf$null_fc <= obs_val, 'G', 'S'))
    }
    cols <- setNames(gen_color_vector(n = 2, name = 'FantasticFox1'), c('S', 'G'))
    # plot_palette(cols)

    p_val <- compare_obs_to_null(obs_val, null_dist = dtf$null_fc)
    
    ggplot(data = dtf, aes_string(x = 'null_fc', 
                                  y = '(..count..)/sum(..count..)', 
                                  fill = 'col_val')) +
      geom_histogram() + 
      ggtitle(as.formula(sprintf("%s~(italic(n)==%d~italic(p)==%s)", 
                                 gsub(' |-', '~', obs[, arm[idx]]), 
                                 p_dat[arm == arm_l, .N],
                                 p_val))) +
      # geom_vline(color = 'red', xintercept = obs_val) +
      scale_x_continuous(name = 'logFC null distribution') +
      ylab('Relative frequency') +
      scale_fill_manual(name = '', values = cols, guide = F)
  })

  
  plot_panel_layout(hists[order(obs$FC)], 
                    plot_direct = F,
                    filename = file.path(plot_dir, 
                      sprintf('ton_fc_permutation_comp_%s_%s.pdf', tp1, tp2)),
                    w = 17, h = 20)
}
