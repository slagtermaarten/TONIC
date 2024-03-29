single_gene_timepoint_comp <- function(gene_symbol = 'CD274',
                                       timepoints = c(1, 2)) {
  allowed_patients <- patient_labels[timepoint %in% timepoints,
                                     .N >= 2, by = patient][V1 == TRUE, patient]
  l_gene_symbol <- gene_symbol
  gene_idx <- exp_levels[, which(gene_symbol == l_gene_symbol)]
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


#' Golgotha plots to compare FC distributions between induction arms
#'
#'
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


prepare_nano_bivariates <- function(x_var = 'Baseline',
                                    y_var = 'Post-induction',
                                    colour_var = 'arm') {
  p_dat <- danaher_scores.m %>%
    { .[!is.na(value) & (timepoint == x_var | timepoint == y_var)] }
  p_dat <- filter_patients(p_dat, colour_var, 'clinical_response')
  allowed_comps <- p_dat[, .(.N == 2), by = c('patient', 'variable')] %>%
    { .[V1 == T, .(patient, variable)] } %>%
    unique
  setkeyv(p_dat, c('patient', 'variable'))
  p_dat <- dcast(p_dat[allowed_comps],
                 formula = patient + arm + clinical_response + variable ~ timepoint)
  p_dat[ , 'logfc' := (get(y_var) + 1) - (get(x_var) + 1)]
  p_dat[, uniqueN(patient)]
  p_dat <- p_dat[apply(p_dat, 1, function(x) !any(is.na(x)))]
  return(p_dat)
}


plot_bivariates <- function(x_var = 'baseline', y_var = 'post.induction',
                            colour_var = 'arm') {
  p_dat <- prepare_nano_bivariates(x_var = x_var, y_var = y_var,
                          colour_var = colour_var)
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
  dat <- filter_patients(dat, color_var, 'response')
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


prep_geneset_parallel_coords <- function(gene_set = 'apm',
                                         colour_var = 'response',
                                         facet_var = NULL) {
  ## Only compute medians for factor combinations with sufficient patients
  by_vars <- setdiff(c('timepoint', colour_var, facet_var),
                     c('clinical_response', 'response'))
  by_vars <- c('timepoint', colour_var, facet_var)
             
  # danaher_scores.m[, table(variable)]
  allowed <- danaher_scores.m[!is.na(value) & variable == gene_set] %>%
    { .[, .('N' = uniqueN(patient)), by = by_vars] } %>%
    { .[N > 1] }
  setkeyv(danaher_scores.m, by_vars)

  sum_dat <- filter_patients(
    danaher_scores.m[allowed][variable %in% gene_set],
    colour_var, facet_var)

  # danaher_scores.m[allowed][!is.na(value)][variable == gene_set][arm == 'Cisplatin']
  sum_dat <- sum_dat[variable == gene_set] %>%
    { .[, .('value' = median(value, na.rm = T), 'patient' = 'median'),
        by = by_vars] }

  if (!is.null(colour_var) && colour_var %in% colnames(sum_dat)) {
    sum_dat <- sum_dat[!is.na(get(colour_var))]
  }

  return(sum_dat)
}


plot_parallel_coords <- function(p_dat,
                                 sum_dat = NULL,
                                 colour_var = 'response',
                                 group_var = 'patient',
                                 facet_var = NULL,
                                 size_var = NULL,
                                 timepoint_v = 'timepoint',
                                 man_colors = tonic_color_palettes[[colour_var]],
                                 swarm_width = .2,
                                 sum_colors = T,
                                 unswarm_zero = F,
                                 overlay_boxplot = F,
                                 point_alpha = .7,
                                 line_alpha = point_alpha / 2 * 1,
                                 median_line_alpha = .8,
                                 filter_vals = T,
                                 title = '') {
  if (!require('vipor')) {
    devtools::install_github('sherrillmix/vipor')
  }
  p_dat <- p_dat[!is.na(value)]
  if (filter_vals) {
    p_dat <- filter_patients(p_dat, colour_var, group_var, facet_var, size_var)
  }
  p_dat[, vipor::offsetX(value, width = swarm_width)]
  p_dat[, as.integer(get(timepoint_v))]
  stopifnot(p_dat[, class(get(timepoint_v))] == 'factor')
  p_dat[, 'x_coord_base' := as.integer(get(timepoint_v)), by = c(timepoint_v, facet_var)]
  p_dat[, 'x_coord' := vipor::offsetX(value, width = swarm_width) +
                     as.integer(get(timepoint_v)),
        by = c(timepoint_v, facet_var)]
  stopifnot(p_dat[, mean(is.na(x_coord))] < 1)

  if (unswarm_zero) {
    p_dat[value == 0, x_coord := as.integer(get(timepoint_v)),
          by = c(timepoint_v, facet_var)]
  }

  if (!is.null(facet_var)) {
    fl <- length(facet_var)
  } else {
    fl <- 1
  }

  if ('clinical_response' %in% colnames(p_dat)) {
    p_dat[is.na(clinical_response), clinical_response := 'NA']
  }

  p <- ggplot(p_dat, aes_string(x = 'x_coord', y = 'value',
                                labels = NULL, group = group_var,
                                size = size_var)) +
    # ggbeeswarm::geom_quasirandom(alpha = point_alpha) +
    geom_point(alpha = point_alpha) +
    geom_line(alpha = line_alpha) +
    # scale_x_discrete(name = '', labels = timepoints) +
    ylab('Gene expression score') +
    ggtitle(title)

  if (overlay_boxplot && F) {
    ## FIX ME
    if (!is.null(colour_var)) {
      group_var <- colour_var
    } else {
      group_var <- 'x_coord_base'
    }
    p <- p + geom_boxplot(aes_string(x = 'x_coord_base', y = 'value', 
                                     group = group_var, 
                                     fill = 'paste(x_coord_base, colour_var)'),
                          alpha = .5, width = swarm_width, outlier.shape = NA)
  }

  p <- p + scale_x_continuous(name = '', minor_breaks = c(),
             breaks = p_dat[, seq_along(levels(get(timepoint_v)))],
             labels = p_dat[, levels(get(timepoint_v))])
  if (timepoint_v == 'timepoint') {
    p <- p + rotate_x_labels(45)
  } else if (timepoint_v == 'blood_timepoint') {
  }

  if (!is.null(sum_dat)) {
    if (sum_colors == F) {
      sum_dat[, 'x_coord' := as.integer(get(timepoint_v)),
              by = c(timepoint_v, facet_var)]
      if (colour_var %nin% colnames(sum_dat)) {
        sum_dat[, (colour_var) := NA]
      }

      p <- p + geom_point(aes_string(group = facet_var), shape = 21,
                          size = 3, alpha = median_line_alpha, colour = 'black',
                          fill = 'black',
                          data = sum_dat, size = 2)
      p <- p + geom_line(aes_string(group = facet_var), shape = 21,
                          size = 2, alpha = median_line_alpha,
                          lineend = 'round', colour = 'black',
                          data = sum_dat, size = 2)
    } else {
      sum_dat[, 'x_coord' := as.integer(get(timepoint_v)),
              by = c(timepoint_v, colour_var, facet_var)]
      if (colour_var %nin% colnames(sum_dat)) {
        sum_dat[, (colour_var) := NA]
      }
      
      if (all(c(facet_var, colour_var) %in% colnames(p_dat))) {
        sum_dat[, 'group_var' := interaction(get(facet_var), get(colour_var))]
      } else {
        sum_dat[, 'group_var' := 'patient']
      }

      p <- p + geom_point(aes_string(group = 'group_var', colour = colour_var), 
                          shape = 21,
                          size = 3, alpha = median_line_alpha, 
                          data = sum_dat, 
                          size = 2)
      p <- p + geom_line(aes_string(group = 'group_var', colour = colour_var), 
                         shape = 21,
                         size = 2, alpha = median_line_alpha,
                         lineend = 'round',
                         data = sum_dat, 
                         size = 2)
    }
  }

  if (!is.null(colour_var)) {
    p <- p + aes_string(fill = colour_var)
    p <- p + scale_fill_manual(name = var_to_label(colour_var),
                               values = man_colors)
    p <- p + aes_string(colour = colour_var)
    p <- p + scale_colour_manual(name = var_to_label(colour_var),
                                 values = man_colors)
  }

  if (!is.null(facet_var)) {
    p <- p + facet_grid(as.formula(sprintf('~ %s', facet_var)))
  }

  return(p)
}


plot_parallel_coords_geneset <- function(gene_set = 'apm',
                                         colour_var = 'clinical_response',
                                         facet_var = NULL, ...) {
  sum_dat <- prep_geneset_parallel_coords(gene_set = gene_set,
                                          colour_var = colour_var,
                                          facet_var = facet_var)
  p_dat <- filter_patients(danaher_scores.m[variable == gene_set],
                           colour_var, facet_var) %>%
    { .[naturalsort::naturalorder(timepoint)] } %>%
    { .[naturalsort::naturalorder(patient)] }
  p <- plot_parallel_coords(p_dat = p_dat,
                            sum_dat = sum_dat,
                            swarm_width = .1,
                            facet_var = facet_var,
                            colour_var = colour_var,
                            title = gene_set,
                            ...)
  return(p)
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

  p_dat <- filter_patients(p_dat, colour_var)
  p_dat[, response := factor(response,
                             levels = levels(response)[c(2, 4, 5, 3, 1)])]
  p_dat[, timepoint := factor(timepoint, levels = timepoints)]
  p_dat[, arm := factor(arm, levels = levels(arm)[c(4, 5, 2, 1, 3)])]
  ## Ensure response factor in correct order
  # p_dat[, levels(response)]
  # stopifnot(p_dat[, all(levels(response) == c('CR', 'PR', 'SD', 'PD'))])

  sum_dat <- p_dat %>%
    { .[, .('value' = median(value, na.rm = T), 'patient' = 'median'),
        by = c('timepoint', colour_var)] } %>%
    { .[!is.na(get(colour_var))] }
  return(list('p_dat' = p_dat, 'sum_dat' = sum_dat))
}


plot_parallel_coords_single_gene <- function(gene = 'CD274',
                                             colour_var = 'response', ...) {

  prep <- prep_gene_parallel_coords(gene = gene,
                                    colour_var = colour_var)
  plot_parallel_coords(p_dat = prep[['p_dat']],
                       sum_dat = prep[['sum_dat']],
                       colour_var = colour_var, title = gene, ...)
}


plot_p_values_induction <- function(m, size_var = '1/p_val',
                                    single_arms_only = T,
                                    tp1 = timepoints[1],
                                    tp2 = timeoints[3]) {
  logfc_range <- m[, max(abs(range(logFC)))] %>% { . * c(-1, 1) }
  if (single_arms_only) m <- m[arm != 'All arms']
  p <- ggplot(m, aes_string(x = 'arm', y = 'gene_set', fill = 'logFC',
                            size = size_var)) +
    geom_tile(size = 0) +
    geom_point() +
    rotate_x_labels(rotate_labels = 90) +
    scale_fill_gradient2(name = 'Log2 FC',
                         low = scales::muted('blue'),
                         high = scales::muted('red'),
                         limits = logfc_range,
                         na.value = 'white') +
    scale_size_continuous(name = 'Unadjusted p-value',
                          range = c(0.1, 2),
                          # trans = 'probability',
                          trans = 'identity',
                          # limits = c(1e-5, fdr_thresh),
                          # labels = function(x) x) +
                          labels = function(x) round(1/x, 3)) +
    scale_y_discrete(name = '', expand = c(0, 0), labels = tonic_cap) +
    scale_x_discrete(name = '', expand = c(0, 0)) +
    ggtitle(m[1, sprintf('%s vs. %s', tp1, tp2)]) +
    theme(legend.position = 'right', legend.direction = 'vertical') +
    rotate_x_labels(45)
  return(p)
}
