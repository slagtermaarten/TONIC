## Labels for timepoints 1, 2 and 3
timepoint_labels <- c('baseline', 'post_induction', 'on_nivo')

## More human readable labels for plotting purposes
axis_subs <- c('tp1' = 'Baseline', 'tp2' = 'Post-induction',
               'tp3' = 'On-nivo',
               'baseline' = 'Baseline', 'post.induction' = 'Post-induction',
               'on.nivo' = 'On-nivo',
               'post_induction' = 'Post-induction',
               'on_nivo' = 'On-nivo')

single_gene_timepoint_comp <- function(gene_symbol = 'CD274',
                                       timepoints = c(1, 2)) {
  allowed_patients <- patient_labels[timepoint %in% timepoints,
                                     .N >= 2, by = patient][V1 == TRUE, patient]
  gene_idx <- exp_levels[, which(gene_symbol == parent.frame(3)$gene_symbol)]
  setkey(patient_labels, patient)
  res <- lapply(timepoints, function(tp) {
    exp_levels[gene_idx,
               patient_labels[allowed_patients][timepoint == tp, 
                                                as.character(filename)],
               with = F] %>% unlist %>% setNames(., NULL)
  }) %>% { data.table('patient' = allowed_patients,
                      'tp_a' = .[[1]], 'tp_b' = .[[2]]) }
  res <- merge(res, unique(patient_labels[, .(patient, response, arm)]),
               by = 'patient')
  res[, 'fc' := log2(tp_b) - log2(tp_a)]
  setnames(res, c('tp_a', 'tp_b'), timepoint_labels[timepoints])
  res[, 'variable' := gene_symbol]
  return(res)
}


plot_timepoint_comparison <- function(res,
                                      x_var = NA, y_var = NA,
                                      colour_var = 'arm',
                                      N_lines = 17,
                                      linesize = 1,
                                      title_trans = toupper,
                                      x_lab = x_var,
                                      y_lab = y_var) {

  if (is.na(x_var) || is.na(y_var)) {
    inferred <- intersect(colnames(res), timepoint_labels)
    x_var <- inferred[1]
    y_var <- inferred[2]
  }

  perc_data = res[, .('median_x' = median(get(x_var), na.rm = T),
                      'median_y' = median(get(y_var), na.rm = T),
                      'q0_x' = quantile(get(x_var), probs = 0, na.rm = T),
                      'q0_y' = quantile(get(y_var), probs = 0, na.rm = T),
                      'q100_x' = quantile(get(x_var), probs = 1, na.rm = T),
                      'q100_y' = quantile(get(y_var), probs = 1, na.rm = T),
                      'q25_x' = quantile(get(x_var), probs = .25, na.rm = T),
                      'q75_x' = quantile(get(x_var), probs = .75, na.rm = T),
                      'q25_y' = quantile(get(y_var), probs = .25, na.rm = T),
                      'q75_y' = quantile(get(y_var), probs = .75, na.rm = T)),
                  by = colour_var]

  fcs <- res[, quantile(fc, probs = c(0, 1), na.rm = T)]
  fcs <- perc_data[, quantile(log2(median_y + 1) - log2(median_x + 1),
                              probs = c(0, 1), na.rm = T)]
  ## Max FC either up- or downward
  max_fc <- ceiling(max(abs(fcs)))
  max_fc <- 5

  # if (!is.na(N_lines) && is.na(slope_d)) {
  #   slope_d <- range_size / (N_lines - 1)
  # } else if (is.na(N_lines) && !is.na(slope_d)) {
  #   N_lines <-  range_size / slope_d + 1
  # } else if (is.na(N_lines) && is.na(slope_d)) {
  #   stop('Define either N_lines or slope_d')
  # }
  # slopes <- slope_d * seq(floor(fcs[1]), ceiling(fcs[2]), length.out = N_lines)
  slopes <- 2^seq(-max_fc, max_fc, length.out = N_lines)
  # message(slopes)
  x_lab <- ifelse(x_lab %in% names(axis_subs), axis_subs[x_lab], x_lab)
  y_lab <- ifelse(y_lab %in% names(axis_subs), axis_subs[y_lab], y_lab)

  p <- ggplot(res, aes_string(x = x_var, y = y_var, colour = colour_var)) +
    # geom_vline(linetype = 'dashed',
    #            xintercept = median(res[, get(x_var)], na.rm = T), alpha = .5) +
    # geom_hline(linetype = 'dashed',
    #            yintercept = median(res[, get(y_var)], na.rm = T), alpha = .5) +
    # geom_point(alpha = nrow(res)^(-1/3)) +
    # geom_density_2d(alpha = .1) +
    # theme_bw() +
    # maartenutils::theme_ms(base_size = 8, legend.position = 'right') +
    fasanalysis::theme_fas(base_size = 8, legend.position = 'right') +
    ggplot2::theme(aspect.ratio = 1) +
    # scale_x_continuous(name = x_lab,
    #                    limits = c(perc_data[, min(q0_x, q0_y)],
    #                               perc_data[, max(q100_x, q100_y)])) +
    # scale_y_continuous(name = y_lab,
    #                    limits = c(perc_data[, min(q0_x, q0_y)],
    #                               perc_data[, max(q100_x, q100_y)])) +
    scale_x_continuous(name = x_lab) +
    scale_y_continuous(name = y_lab) +
    coord_cartesian(xlim = c(perc_data[, min(q25_x, q25_y)],
                             perc_data[, max(q75_x, q75_y)]),
                    ylim = c(perc_data[, min(q25_y, q25_x)],
                             perc_data[, max(q75_y, q75_x)])) +
    scale_colour_discrete()

  if (res[, uniqueN(variable) == 1]) {
    p <- p + ggtitle(title_trans(res[, unique(variable)]))
  }

  dist_vals <- list(mean = 0, sd = 10)
  for (sl in slopes) {
    ## Fade away lines according to normal distribution
    ## Fix identity line at alpha of .3
    alpha_val <- .3 / do.call('dnorm', c(list(x = 1), dist_vals)) *
                  do.call('dnorm', c(list(x = sl), dist_vals))
    print(sl)
    print(alpha_val)
    ## Fix me
    p <- p + geom_abline(slope = sl, intercept = 0,
                         colour = ifelse(eps(log2(sl) %% 1, 0),
                                         'red', 'black'),
                         linetype = ifelse(eps(sl, 1), 'dashed', 'dotted'),
                         alpha = alpha_val)
  }

  p <- p + geom_segment(aes(x = q25_x, y = median_y,
                            xend = q75_x, yend = median_y),
                        lineend = 'butt', size = linesize, data = perc_data)
  p <- p + geom_segment(aes(x = median_x, yend = q25_y,
                            xend = median_x, y = q75_y),
                        lineend = 'butt', size = linesize, data = perc_data)
  return(p)
}


plot_bivariates <- function(x_var = 'baseline', y_var = 'post.induction',
                            colour_var = 'arm') {
  p_dat <- danaher_scores.m %>%
                 .[!is.na(value) & (timepoint == x_var | timepoint == y_var)]
  allowed_comps <- p_dat[, .(.N == 2), by = c('patient', 'variable')] %>%
                     .[V1 == T, .(patient, variable)] %>% unique
  setkeyv(p_dat, c('patient', 'variable'))
  p_dat <- dcast(p_dat[allowed_comps],
                 formula = patient + arm + response + variable ~ timepoint)
  p_dat[ , 'logfc' := (get(y_var) + 1) - (get(x_var) + 1)]
  p_dat[, uniqueN(patient)]
  p_dat <- p_dat[apply(p_dat, 1, function(x) !any(is.na(x)))]

  library(grid)
  library(maartenutils)

  ## Make plot for each signature
  plots <- lapply(p_dat[, unique(variable)], function(gs) {
    p <- p_dat %>% .[variable == gs] %>%
      plot_timepoint_comparison(x_var = x_var, y_var = y_var,
                                colour_var = colour_var)
    if (gs != p_dat[, unique(variable)[1]]) {
      p <- p + ggplot2::theme(legend.position = 'none')
    }
    return(p)
  })
  plots_per_page <- 15
  pages <- ceiling(length(plots) / plots_per_page)

  for (p in 1:pages) {
    plot_idx <- ((p - 1) * plots_per_page + 1) : ((p) * plots_per_page)
    filename <- sprintf('%s/bivariate_%s-%s_%s-%s.pdf', 
                        img_dir,
                        x_var, y_var,
                        colour_var, p)
    maartenutils::plot_panel_layout(plots[plot_idx], ncol = 3, label_size = 8,
                                    w = 17.4, h = 25, filename = filename)
    maartenutils::sys_file_open(filename)
  }
}


plot_PCA_axes <- function(dat, tp, color_var = '', labellings) {
  library(ggbiplot)
  pca_dat <- prcomp(dat, scale = T)

  ## All PC pairs: 1&2, 3&4, 5&6, etc.
  N_plots <- ceiling(ncol(pca_dat$rotation)/2)
  PC_pairs <- 
    base::split(1:ncol(pca_dat$rotation), rep(1:N_plots, each = 2))[1:12]

  plots <- lapply(PC_pairs, function(li) {
    if (length(li) != 2 || any(is.na(li))) return(NULL)
    g <- ggbiplot::ggbiplot(pca_dat, obs.scale = 1, var.scale = 1, 
                  choices = li,
                  groups = labellings, 
                  ellipse = TRUE, 
                  var.axes = F, 
                  circle = TRUE)
    g <- g + scale_color_discrete(name = '')
    g <- g + theme(legend.direction = 'horizontal', legend.position = 'top')
    return(g)
  })

  fn <- sprintf('plots/PCA_decomp_%s_%s.pdf', tp, color_var)
  plot_panel_layout(plots, ncol = 3, filename = fn, 
                    w = 40, h = 40, label_size = 8)
  plot_panel_layout(plots, ncol = 3, filename = NULL, 
                    w = 40, h = 40, label_size = 8)
  # maartenutils::sys_file_open(fn)
}


prep_geneset_parallel_coords <- function(geneset = 'apm', 
                                         colour_var = 'response') {
  ## Subselect 'complete' patients
  # allowed_comps <- danaher_scores.m[!is.na(value) & !is.na(response)] %>% 
  #   { .[, (.N == 3), by = c('patient', 'variable')] } %>%
  #   { .[V1 == T, .(patient, variable)] } %>% unique
  # setkeyv(danaher_scores.m, c('patient', 'variable'))
  # p_dat <- danaher_scores.m[allowed_comps]
  # p_dat[, lapply(.SD, class)]
  # p_dat[, table(arm, variable, timepoint)]
  # p_dat[, table(timepoint)]
  # setnames(p_dat, 'variable', 'geneset')
  sum_dat <- danaher_scores.m[variable == geneset] %>%
    { .[, .('value' = median(value, na.rm = T), 'patient' = 'median'), 
        by = c('timepoint', colour_var)] } %>%
    { .[!is.na(get(colour_var))] }
  return(sum_dat)
}


plot_parallel_coords <- function(p_dat, sum_dat, colour_var = 'response',
                                 title = '') {
  ggplot(p_dat, 
         aes_string(x = 'timepoint', y = 'value', 
                    fill = colour_var,
                    labels = NULL, group = 'patient', colour = colour_var)) +
    geom_point(alpha = .4) + geom_line(alpha = .2) +
    geom_point(aes_string(group = colour_var), shape = 21, 
               size = 4, alpha = .8, colour = 'black', 
               data = sum_dat, size = 2) + 
    geom_line(aes(group = NULL), alpha = .6, data = sum_dat) +
    scale_x_discrete(name = '', labels = timepoints) +
    scale_y_continuous(name = 'Gene expression score') +
    scale_fill_brewer(palette = "Set2", direction = 1) +
    scale_colour_brewer(palette = "Set2", direction = 1) +
    rotate_x_labels(45) +
    ggtitle(title)
}


plot_parallel_coords_geneset <- function(geneset = 'apm',
                                         colour_var = 'response') {
  sum_dat <- prep_geneset_parallel_coords(geneset = geneset, 
                                          colour_var = colour_var)
  p_dat <- danaher_scores.m[variable == geneset]
  plot_parallel_coords(p_dat = p_dat, sum_dat = sum_dat, 
                       colour_var = colour_var, title = geneset)
}


prep_gene_parallel_coords <- function(gene = 'CD274', 
                                      colour_var = 'response') {
  gene_idx <- exp_levels[, which(gene_symbol == gene)]
  p_dat <- exp_levels[gene_idx, rev(colnames(exp_levels)[-1]), with = F] %>% 
    { log2(. + 1) } %>%
    t %>%
    { data.frame(value = ., filename = rownames(.)) } %>%
    { merge(., patient_labels[, .(filename, patient, 
                                  timepoint, arm, response)]) } %>%
    as.data.table

  p_dat[, response := factor(response, 
                             levels = levels(response)[c(2, 4, 5, 3, 1)])]
  p_dat[, timepoint := factor(timepoint, levels = timepoints)]
  p_dat[, arm := factor(arm, levels = levels(arm)[c(4, 5, 2, 1, 3)])]
  ## Ensure response factor in correct order
  stopifnot(p_dat[, all(levels(response) == c('CR', 'PR', 'SD', 'PD', ''))])

  sum_dat <- p_dat %>% 
    { .[, .('value' = median(value, na.rm = T), 'patient' = 'median'), 
        by = c('timepoint', colour_var)] } %>%
    { .[!is.na(get(colour_var))] }
  return(list('p_dat' = p_dat, 'sum_dat' = sum_dat))
}


plot_parallel_coords_single_gene <- function(gene = 'CD274',
                                             colour_var = 'response') {

  prep <- prep_gene_parallel_coords(gene = gene, 
                                    colour_var = colour_var)
  plot_parallel_coords(p_dat = prep[['p_dat']], sum_dat = prep[['sum_dat']], 
                       colour_var = colour_var, title = gene)
}
