plot_cor_mat <- function(tp = timepoints,
                         cormethod = 'spearman',
                         sig_level = c(.05, .0005, 000005)) {
  library(corrplot)
  p_dat <-
    patient_labels[, .(timepoint,
                       tis_score,
                       productive_clonality,
                       adaptive_t_cells,
                       sample_amount_ng,
                       cd8_mm2,
                       pd_l1_tumor,
                       pd_l1_immunoinfiltrate,
                       ca15_3,
                       ldh,
                       crp,
                       s_til,
                       ichao1,
                       observed_richness,
                       efron_thisted_estimator,
                       daley_smith_estimator)] %>%
    as.data.frame() %>%
    dplyr::filter(timepoint %in% tp) %>%
    dplyr::select(-timepoint) %>%
    { .[apply(., 1, function(x) all(!is.na(x))), ] } %>%
    as.matrix

  ## Format column names a bit
  colnames(p_dat) <- gsub('^s_', 'stromal ', colnames(p_dat))
  colnames(p_dat) <- gsub('pd_l1', 'PD-L1', colnames(p_dat))
  colnames(p_dat) <- gsub('ca15_3', 'ca15.3', colnames(p_dat))
  colnames(p_dat) <- simple_cap(gsub('_', ' ', colnames(p_dat)),
    caplist = c('mm2', 'Ichao1', 'cells', 'TIS', 'score', 's', 'ng', 'LDH',
                'CRP', 'TIL'))

  corrplot(cor(p_dat, method = cormethod),
           p.mat = cor.mtest(p_dat, method = cormethod)$p,
           addrect = 6,
           tl.col = 'gray40',
           pch.col = 'white',
           pch.cex = .9,
           # order = 'AOE',
           order = "hclust",
           # tl.pos = 'ld',
           # type = 'lower',
           tl.cex = .9,
           sig.level = sig_level,
           insig = "label_sig",
           method = 'color',
           mar = c(0,0,6,0), # http://stackoverflow.com/a/14754408/54964
           title = sprintf('%s - %s',
                           paste(tp, collapse = ', '),
                           simple_cap(cormethod)))
}


plot_parallel_adaptive <- function(p_var = 'efron_thisted_estimator',  
                                   facet_var = 'arm', compartment = 'tumor',
                                   ...) {
  compartment <- match.arg(compartment, choices = c('blood', 'tumor'))
  if (compartment == 'tumor') {
    timepoint_v <- 'timepoint'
    p_dat <- patient_labels[, c('patient', timepoint_v, 'clinical_response',
                                'arm', p_var), with = F]
  } else if (compartment == 'blood') {
    timepoint_v <- 'blood_timepoint'
    p_dat <- blood_adaptive[, c('patient', timepoint_v, 'clinical_response',
                                'arm', p_var), with = F]
  }
  setnames(p_dat, p_var, 'value')
  p_dat <- p_dat[!is.na(value)]

  plot_parallel_coords(p_dat, facet_var = facet_var,
                       timepoint_v = timepoint_v,
                       colour_var = 'clinical_response') +
    scale_y_continuous(name = var_to_label(p_var))
}


prepare_adaptive_FC <- function(
  # facet_var = NULL,
  # facet_var = rev(c('arm', 'clinical_response')),
  # facet_var = 'clinical_response',
  dtf = patient_labels,
  facet_var = 'arm',
  alpha_lev = .8,
  # var1 = 'pielou_evenness',
  var1 = 'sample_clonality',
  var2 = 'efron_thisted_estimator') {
  aggregate_timepoints <- function(dtf = patient_labels, 
                                   lvar = 'efron_thisted_estimator') {
    p_dat <- Reduce(function(x, y) merge(x, y, all = TRUE, by = 'patient'),
     list(dtf[timepoint == timepoints[1],
                         c('patient', lvar), with = F][!is.na(get(lvar))],
          dtf[timepoint == timepoints[2],
                         c('patient', lvar), with = F][!is.na(get(lvar))],
          dtf[timepoint == timepoints[3],
                         c('patient', lvar), with = F][!is.na(get(lvar))]))
    p_dat <- setnames(p_dat, c('patient',
                               sprintf('%s_%s', lvar, c('tp1', 'tp2', 'tp3'))))
    p_dat[, (sprintf('%s_12', lvar)) := (log2(get(sprintf('%s_tp2', lvar))) -
                                           log2(get(sprintf('%s_tp1', lvar))))]
    p_dat[, (sprintf('%s_13', lvar)) := (log2(get(sprintf('%s_tp3', lvar))) -
                                           log2(get(sprintf('%s_tp1', lvar))))]
    p_dat <- 
      rbind(cbind(p_dat[, .(patient, 'FC' = 0)], 
                  'timepoint' = timepoints[1]),
            cbind(p_dat[, .(patient, 'FC' = get(sprintf('%s_12', lvar)))], 
                  'timepoint' = timepoints[2]),
            cbind(p_dat[, .(patient, 'FC' = get(sprintf('%s_13', lvar)))], 
                  'timepoint' = timepoints[3]))
    setnames(p_dat, 'FC', lvar)
    return(p_dat)
  }

  ## Merge the FCs of the two variables
  agg_tp <-
    merge(aggregate_timepoints(dtf = dtf, lvar = var1), 
          aggregate_timepoints(dtf = dtf, lvar = var2),
          by = c('patient', 'timepoint'), all = T) %>%
    { .[naturalsort::naturalorder(patient)] }

  ## Merge in patient labels
  agg_tp <- patient_labels %>%
    dplyr::select_('-timepoint') %>%
    dplyr::select_(sprintf('-%s', var1)) %>%
    dplyr::select_(sprintf('-%s', var2)) %>%
    { unique(., by = 'patient') } %>%
    { merge(agg_tp, ., by = c('patient'), all.y = T) }

  ## Restore this variable as this may have been incorrectly merged
  agg_tp[, comb_time_resp := sprintf('%s-%s', timepoint, clinical_response)]
  agg_tp <- agg_tp[!comb_time_resp == 'NA-NR']

  ## Patients with On nivo but no Post-induction, set post-induction to 0
  if ('timepoint' %in% colnames(dtf)) {
    agg_tp[, comb_time_resp := 
         factor(comb_time_resp, 
                levels = patient_labels[, levels(comb_time_resp)])]
    # agg_tp[, is.na(get(var1)), by = patient]
    missing_pats <- 
      agg_tp[!is.na(get(var1)), 
             (timepoints[3] %in% .SD[, timepoint]) & 
             (timepoints[2] %nin% .SD[, timepoint]), by = patient] %>%
             { .[V1 == T, patient] }
    agg_tp[patient %in% missing_pats & timepoint == timepoints[2], (var1) := 0]

    missing_pats <- 
    agg_tp[!is.na(get(var2)), 
           (timepoints[3] %in% .SD[, timepoint]) & 
           (timepoints[2] %nin% .SD[, timepoint]), by = patient] %>%
           { .[V1 == T, patient] }
    agg_tp[patient %in% missing_pats & timepoint == timepoints[2], (var2) := 0]
    agg_tp[, timepoint := factor(timepoint, levels = timepoints)]
  } else if ('blood_timepoint' %in% colnames(dtf)) {
    # agg_tp[, comb_time_resp := 
    #      factor(comb_time_resp, 
    #             levels = patient_labels[, levels(comb_time_resp)])]
    browser()
  }

  ## Patients with Post-induction but no On nivo, remove on-nivo entry
  # missing_pats <- 
  # agg_tp[!is.na(get(var1)), 
  #        (timepoints[2] %in% .SD[, timepoint]) & 
  #        (timepoints[3] %nin% .SD[, timepoint]), by = patient][V1 == T, patient]
  # agg_tp <- agg_tp[(patient %in% missing_pats) %nand% (timepoint == timepoints[3])] 
  # rm(missing_pats)
  # agg_tp <- agg_tp[timepoint != timepoints[1]]
  agg_tp <- agg_tp[order(patient, timepoint)]
  return(agg_tp)
}



plot_adaptive_FC <- function(dtf = patient_labels,
                             facet_var = 'arm',
                             alpha_lev = .8,
                             var1 = 'sample_clonality',
                             var2 = 'efron_thisted_estimator') {
  agg_tp <- prepare_adaptive_FC(dtf = dtf, facet_var = facet_var, 
                                alpha_lev = alpha_lev, var1 = var1, 
                                var2 = var2)
  s_plot <- ggplot(agg_tp, aes_string(x = var1, y = var2,
                                      group = 'patient',
                                      colour = 'comb_time_resp')) +
    geom_hline(yintercept = 0, color = 'gray20', linetype = 'dashed') + 
    geom_vline(xintercept = 0, color = 'gray20', linetype = 'dashed') + 
    geom_path(alpha = alpha_lev, show.legend = F) +
    geom_point(data = agg_tp[(eps(get(var1), 0) %nand% eps(get(var2), 0))]) +
    # geom_point(size = 2, alpha = alpha_lev) +
    # ggrepel::geom_text_repel(segment.color = 'grey90', color = 'black') +
    # geom_text(color = 'white') +
    # geom_smooth(se = F) +
    scale_fill_manual(values = comb_time_resp_palette) +
    scale_colour_manual(name = '', values = comb_time_resp_palette) +
    scale_shape_manual(name = '',
                       values = setNames(c(1, 2, 3, 1, 2, 3) + 14,
                                         names(comb_time_resp_palette))) +

    scale_x_continuous(var_to_label(var1), trans = 'identity', expand = c(0.1, 0.1)) +
    scale_y_continuous(var_to_label(var2), trans = 'identity',
                       expand = c(0, 0)) +
    theme(legend.position = 'top', legend.direction = 'vertical') +
    gg_legend_alpha_cancel

  if (!is.null(facet_var)) {
    s_plot <- s_plot + facet_wrap(as.formula(sprintf('~ %s',
                                                     paste(facet_var,
                                                           collapse = ' + '))),
                                  nrow = length(facet_var))
  }
  return(s_plot)
}
