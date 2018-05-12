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
    scale_y_continuous(name = var_to_label(p_var)) +
    ggtitle(sprintf('%s', compartment))
}


prepare_adaptive_FC <- function(# facet_var = NULL,
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


plot_FC <- function(agg_tp,
                    facet_var = 'arm',
                    alpha_lev = .8,
                    colour_var = 'comb_time_resp',
                    plot_points = 'nonzero',
                    var1 = 'sample_clonality',
                    var2 = 'efron_thisted_estimator') {
  if (is.null(agg_tp)) return(NULL)
  s_plot <- ggplot(agg_tp, aes_string(x = var1, y = var2,
                                      group = 'patient')) +
    geom_hline(yintercept = 0, color = 'gray20', linetype = 'dashed') +
    geom_vline(xintercept = 0, color = 'gray20', linetype = 'dashed')
    # geom_point(size = 2, alpha = alpha_lev) +
    # ggrepel::geom_text_repel(segment.color = 'grey90', color = 'black') +
    # geom_text(color = 'white') +
    # geom_smooth(se = F) +

  if (plot_points == 'nonzero') {
    s_plot <- s_plot +
      geom_point(data = agg_tp[(eps(get(var1)) %nand% eps(get(var2)))]) +
      geom_path(alpha = alpha_lev, show.legend = F)
  } else if (plot_points == 'all') {
    s_plot <- s_plot + geom_point()
  }

  if (!is.null(colour_var) && colour_var == 'comb_time_resp') {
    s_plot <- s_plot + aes_string(shape = colour_var)
    pal_mask <- c(2,3,5,6)
    pal_mask <- c(1, 2, 3, 4, 5, 6)
    s_plot <- s_plot +
      aes_string(colour = colour_var) +
      scale_fill_manual(values = comb_time_resp_palette[pal_mask]) +
      scale_colour_manual(name = '', values = comb_time_resp_palette[pal_mask]) +
      scale_shape_manual(name = '',
               values = setNames(c(1, 2, 3, 1, 2, 3) + 14,
                                 names(comb_time_resp_palette))[pal_mask])
  }

  s_plot <- s_plot +
    scale_x_continuous(var_to_label(var1), trans = 'identity',
                       expand = c(0.1, 0.1)) +
    scale_y_continuous(var_to_label(var2), trans = 'identity',
                       expand = c(0.1, 0.1)) +
    theme(legend.position = 'top', legend.direction = 'vertical') +
    gg_legend_alpha_cancel

  if (!is.null(facet_var)) {
    s_plot <- s_plot +
      facet_wrap(as.formula(sprintf('~ %s',
                                    paste(facet_var, collapse = ' + '))),
                 nrow = length(facet_var))
  }
  return(s_plot)
}

compute_TCR_FCs <- function(patient = 'pat_11',
                            tp1 = 'Baseline',
                            tp2 = 'Post-induction',
                            tp3 = '-2',
                            tp4 = '0',
                            y_var = 'productive_frequency') {
  if (!exists('arr')) read_adaptive_seqs()

  suf_data <- all(c(tp1, tp2, tp3, tp4) %in%
                  arr[patient == parent.frame(3)$patient,
                      as.character(timepoint)])

  # arr[patient == parent.frame(3)$patient, unique(timepoint)]

  if (!suf_data) {
    mymessage('compute_TCR_FCs',
              sprintf('Skipping %s for tps %s %s %s %s %s',
                      patient, tp1, tp2, tp3, tp4))
    return(NULL)
  }

  tp1 <- as.character(tp1)
  tp2 <- as.character(tp2)
  tp3 <- as.character(tp3)
  tp4 <- as.character(tp4)

  FCs_l <- arr[patient == parent.frame(3)$patient] %>%
    filter(as.character(timepoint) %in% c(tp3, tp4, tp1, tp2)) %>%
    { .[, {
          if (.SD[timepoint %in% c(tp1, tp2), .N] > 0) {
            fc1 = .SD[timepoint == tp2,
                      log2(replace_NA(get(y_var), 0) + 1)] -
                  .SD[timepoint == tp1,
                      log2(replace_NA(get(y_var), 0) + 1)]
          } else {
            fc1 = as.numeric(NA)
          }
          if (.SD[timepoint %in% c(tp3, tp4), .N] > 0) {
            fc2 = .SD[timepoint == tp4,
                      log2(replace_NA(get(y_var), 0) + 1)] -
                  .SD[timepoint == tp3,
                      log2(replace_NA(get(y_var), 0) + 1)]
          } else {
            fc2 = as.numeric(NA)
          }

          ad1 = .SD[timepoint == tp2, replace_NA(get(y_var), 0)] -
                .SD[timepoint == tp1, replace_NA(get(y_var), 0)]
          ad2 = .SD[timepoint == tp4, replace_NA(get(y_var), 0)] -
                .SD[timepoint == tp3, replace_NA(get(y_var), 0)]

          list('arm' = arm,
               'clinical_response' = clinical_response,
               'response' = response,
               'fc1' = fc1, 'fc2' = fc2,
               'exp1' = ad1, 'exp2' = ad2)
        }, by = .(patient, amino_acid)] }
  return(FCs_l)
}


compute_tp_comp_FCs <- function(
  tp1 = 'Baseline', tp2 = 'Post-induction', tp3 = '-2', tp4 = '0',
  overwrite = T, y_var = 'normalized_frequency') {
  if (!exists('arr')) read_adaptive_seqs()

  plyr::llply(arr[, naturalsort::naturalsort(auto_name(unique(patient)))],
              function(patient) {
    fn <- get_FC_fn(patient, tp1, tp2, tp3, tp4, y_var)
    if (file.exists(fn) && !overwrite) {
      return(readRDS(fn))
    }
    FCs <- compute_TCR_FCs(patient = patient,
                           tp1 = tp1, tp2 = tp2, tp3 = tp3, tp4 = tp4, 
                           y_var = y_var)
    if (is.null(FCs)) return(NULL)
    saveRDS(FCs, file = fn)
    # FC_p <- plot_FC(FCs, facet_var = NULL, alpha_lev = .2,
    #                 colour_var = NULL,
    #                 plot_points = 'all',
    #                 var1 = 'exp1', var2 = 'exp2')
  }, .parallel = !local_run)
}


timepoint_s_dat <- function(timepoint) {
  stopifnot(length(timepoint) == 1)
  if (timepoint %in% blood_timepoints) {
    s_dat <- blood_adaptive
    s_dat <- s_dat[as.character(blood_timepoint) == parent.frame(3)$timepoint]
  } else if (timepoint %in% timepoints) {
    s_dat <- patient_labels
    s_dat <- s_dat[as.character(timepoint) == parent.frame(3)$timepoint]
  } else {
    stopf('Unknown timepoint %s', timepoint)
  }
  return(s_dat)
}


#' Wilcoxon rank sum to test for association with response at single timepoint
#'
#'
test_adaptive_association <- function(measures = c('sample_clonality',
                                                   'efron_thisted_estimator',
                                                   'adaptive_t_cells'),
                                      timepoint = 'On nivo',
                                      y_var = 'clinical_response',
                                      facet_var = NULL,
                                      patient_ids = patient_labels[, unique(patient)]) {
  stopifnot(timepoint %in% c(timepoints, blood_timepoints))
  s_dat <- timepoint_s_dat(timepoint)
  p_dat <- s_dat[as.character(patient) %in% patient_ids]
  if (null_dat(p_dat)) return(NULL)
  comp_levels <- p_dat[, levels(get(y_var))]

  res <- rbindlist(lapply(measures, function(measure) {
    if (measure %nin% colnames(p_dat)) {
      warningf('%s not found', measure)
      return(NULL)
    }

    dtf <- p_dat[, .('measure' = measure,
                     'p_val' = tryCatch(wilcox.test(as.formula(sprintf('%s ~ %s',
                                                              measure, y_var)),
                                           data = .SD)$p.val,
                                  error = function(e) { print(e); return(NA) }),
                     'log2FC' = median(.SD[get(y_var) == comp_levels[2],
                                       log2(get(measure))],
                                       na.rm = T) -
                                median(.SD[get(y_var) == comp_levels[1],
                                       log2(get(measure))],
                                       na.rm = T)),
                by = facet_var]
    return(dtf)
  }), fill = T)
  # res[, 'p_val.bh' := p.adjust(p_val, method = 'BH'), by = facet_var]
  res[, 'p_val.fdr' := p.adjust(p_val, method = 'fdr'), by = facet_var]
  # res[, 'p_val.bonferroni' := p.adjust(p_val, method = 'bonferroni'),
  #     by = facet_var]
  res[, 'timepoint' := timepoint]
  first_cols <- c('timepoint', 'measure')
  setcolorder(res, c(first_cols, setdiff(colnames(res), first_cols)))
  return(res)
}


#' Plot single timepoint summary stat comparisons
#'
#'
plot_single_timepoint_r_assocation <- function(
  measures = c('sample_clonality', 'efron_thisted_estimator',
               'adaptive_t_cells'),
  l_timepoint = 'Baseline', facet_var = 'arm') {
  s_dat <- timepoint_s_dat(l_timepoint)
  lapply(measures, function(measure) {
    p <- ggplot(s_dat,
           aes_string(x = 'clinical_response', y = measure,
                      fill = 'clinical_response')) +
      geom_boxplot() +
      ggpubr::stat_compare_means(label = 'p.signif') +
      scale_fill_manual(name = '', values = resp_colors) +
      scale_y_continuous(name = sprintf('%s at %s', var_to_label(measure),
                                        l_timepoint)) +
      xlab('Clinical response')
    if (!is.null(facet_var)) {
      p <- p + facet_grid(as.formula(sprintf('~ %s ', facet_var)))
    }
    return(p)
  })
}


#' Plot TCR abundance comparison between blood and tumor compartments
#'
#'
plot_tp_comp_direct <- function(tp1 = 'Baseline', tp2 = '-2',
                                y_var = 'normalized_frequency') {
  plyr::llply(arr[, naturalsort::naturalsort(auto_name(unique(patient)))],
              function(patient) {
    suf_data <- arr[patient == parent.frame(3)$patient,
                    all(c(tp1, tp2) %in% timepoint)]
    if (!suf_data) return(NULL)
    t_dat <- arr[patient == parent.frame(3)$patient &
                 timepoint %in% c(tp1, tp2)]

    p_dat <- merge(t_dat[timepoint == tp1,
                   c(y_var, 'amino_acid')],
                   t_dat[timepoint == tp2,
                   c(y_var, 'amino_acid')],
                   all = TRUE, by = 'amino_acid')

    check_timepoint <- function(tp) {
      if (grepl('-*\\d{1,2}', tp2)) return(sprintf('day %s', tp))
      else return(tp)
    }
    tp1 <- check_timepoint(tp1)
    tp2 <- check_timepoint(tp2)
    max_val <- c(p_dat[, range(productive_frequency.x, na.rm = T)],
      p_dat[, range(productive_frequency.y, na.rm = T)]) %>%
    { max(.) }

    ggplot(p_dat, aes(x = productive_frequency.x, y = productive_frequency.y)) +
      geom_count(alpha = .1) +
      scale_x_continuous(name = sprintf('%s %s', var_to_label(y_var), tp1),
                         limits = c(0, max_val)) +
      scale_y_continuous(name = sprintf('%s %s', var_to_label(y_var), tp2),
                         limits = c(0, max_val)) +
      geom_hline(yintercept = 0, color = 'gray20', linetype = 'dashed') +
      geom_vline(xintercept = 0, color = 'gray20', linetype = 'dashed') +
      ggtitle(sprintf('%s - %s - %s',
                      patient,
                      patient_labels[patient == parent.frame(3)$patient,
                                     unique(arm)],
                      patient_labels[patient == parent.frame(3)$patient,
                                     unique(clinical_response)])) +
      theme(legend.position = 'right', legend.direction = 'vertical',
            aspect.ratio = 1)
  })
}


format_patient_title <- function(patient) {
  # patient_labels[patient == parent.frame(3)$patient]
  sprintf('%s - %s - %s',
          patient,
          patient_labels[patient == parent.frame(3)$patient,
                         unique(arm)],
          patient_labels[patient == parent.frame(3)$patient,
                         unique(clinical_response)]) %>%
  unique()
}


#' Get filename of cached FC computation
#'
#'
get_FC_fn <- function(patient, tp1, tp2, tp3, tp4, y_var = '') {
  file.path(rds_dir, 'FCs', sprintf('%s_FCs_%s_%s_%s%s.rds',
                                    patient, tp1, tp2, tp3, tp4,
                                    ifelse(y_var == '', '', 
                                           sprintf('_%s', y_var))))
}


#' Compare TCR abundance FCs in two compartments
#'
#'
plot_tp_comp_FCs <- function(tp1 = 'Baseline', tp2 = 'Post-induction',
                             tp3 = '-2', tp4 = '0', 
                             y_var = 'normalized_frequency') {
  read_adaptive_seqs()

  plyr::llply(arr[, naturalsort::naturalsort(auto_name(unique(patient)))],
              function(patient) {
    fn <- get_FC_fn(patient, tp1, tp2, tp3, tp4, y_var)
    if (file.exists(fn)) {
      FCs <- readRDS(fn)
    } else {
      return(NULL)
    }

    ggplot(FCs, aes(x = exp1, y = exp2)) + geom_count(alpha = .1) +
      xlab(sprintf('Expansion in tumor (%s vs. %s)', tp1, tp2)) +
      ylab(sprintf('Expansion in blood (day %s vs. day %s)', tp3, tp4)) +
      geom_hline(yintercept = 0, color = 'gray20', linetype = 'dashed') +
      geom_vline(xintercept = 0, color = 'gray20', linetype = 'dashed') +
      # ggtitle(sprintf('%s %s-%s vs %s-%s', patient, tp1, tp2, tp3, tp4)) +
      ggtitle(format_patient_title(patient)) +
      theme(legend.position = 'right', legend.direction = 'vertical',
            aspect.ratio = 1)
  })
}


#' Compare Adaptive summary stats
#'
#'
compare_adaptive_summary_stats <- function(tp1 = 'Baseline',
                                           tp2 = 'On nivo',
                                           x_var = 'clinical_response',
                                           colour_var = 'clinical_response',
                                           facet_var = 'arm',
                                           comp_measure = 'sample_clonality') {

  t_dat <- patient_labels[, {
    .SD[timepoint == tp2, get(comp_measure)] -
    .SD[timepoint == tp1, get(comp_measure)]
  }, by = patient]
  setnames(t_dat, 'V1', 'value')

  t_dat <- merge(t_dat[!is.na(value)], patient_labels, all.x = T, all.y = F) %>%
    unique(by = c('patient'))

  label_proc <- function(labels) {
    sapply(labels, function(lab) parse(text = lab))
  }

  t_dat[!is.na(value), 'label' := sprintf('%s~(italic(n)==%d)',
                                          gsub(' ', '~', get(x_var)), .N),
        by = c(facet_var, x_var)]
  if (x_var == 'arm') {
    labels <- t_dat[, unique(label)]
    labels_stripped <- t_dat[, gsub('~', ' ',
                                    gsub('(.*)~\\(.*', '\\1', unique(label)))]
    right_order <- labels[match(t_dat[, levels(arm)], labels_stripped)]
    t_dat[, label := factor(label, levels = right_order)]
  }

  p1 <- ggplot(t_dat[!is.na(label)],
               aes_string(x = 'label', y = 'value', fill = colour_var))
  p1 <- p1 + geom_boxplot()
  p1 <- p1 + ggbeeswarm::geom_quasirandom()
  p1 <- p1 + ggpubr::stat_compare_means()
  p1 <- p1 + scale_y_continuous(name = sprintf('%s %s vs. %s',
                                               var_to_label(comp_measure),
                                               tp2, tp1)) +
    scale_x_discrete(name = '', labels = label_proc) +
    rotate_x_labels(45) +
    theme(legend.position = 'none')

  if (!is.null(facet_var)) {
    p1 <- p1 + facet_grid(as.formula(sprintf('~ %s', facet_var)),
                          scales = 'free_x')
  }
  return(p1)
}


replace_NA <- function(vec, replace_val = 0) {
  if (length(vec) > 0) {
    vec[is.na(vec)] <- replace_val
  } else {
    vec <- 0
  }
  return(vec)
}


plot_adaptive_FC <- function(dtf = patient_labels,
                             facet_var = 'arm',
                             alpha_lev = .8,
                             colour_var = 'comb_time_resp',
                             p_timepoints = timepoints,
                             var1 = 'sample_clonality',
                             var2 = 'efron_thisted_estimator') {
  agg_tp <- prepare_adaptive_FC(dtf = dtf, facet_var = facet_var,
                                alpha_lev = alpha_lev, var1 = var1,
                                var2 = var2)
  plot_FC(agg_tp[timepoint %in% p_timepoints],
          facet_var = facet_var, alpha_lev = alpha_lev,
          colour_var = colour_var, var1 = var1, var2 = var2)
}


#' Plot all relevant comps
#'
#'
plot_all_comps <- function(f = compare_adaptive_summary_stats,
                           colour_var = 'arm',
                           x_var = 'arm',
                           facet_var = NULL) {
  list(
  f(tp1 = 'Baseline', tp2 = 'Post-induction',
    comp_measure = 'sample_clonality',
    colour_var = colour_var,
    x_var = x_var,
    facet_var = facet_var),
  f(tp1 = 'Baseline', tp2 = 'On nivo',
    comp_measure = 'sample_clonality',
    colour_var = colour_var,
    x_var = x_var,
    facet_var = facet_var),

  f(tp1 = 'Baseline', tp2 = 'Post-induction',
    comp_measure = 'adaptive_t_cells',
    colour_var = colour_var,
    x_var = x_var,
    facet_var = facet_var),
  f(tp1 = 'Baseline', tp2 = 'On nivo',
    comp_measure = 'adaptive_t_cells',
    colour_var = colour_var,
    x_var = x_var,
    facet_var = facet_var),

  f(tp1 = 'Baseline', tp2 = 'Post-induction',
    comp_measure = 'efron_thisted_estimator',
    colour_var = colour_var,
    x_var = x_var,
    facet_var = facet_var),
  f(tp1 = 'Baseline', tp2 = 'On nivo',
    comp_measure = 'efron_thisted_estimator',
    colour_var = colour_var,
    x_var = x_var,
    facet_var = facet_var))
}


plot_TCR_chronological <- function(patient = 'pat_11',
                                   timepoint_v = 'timepoint',
                                   facet_var = NULL,
                                   colour_var = 'shared_timepoints',
                                   allowed_timepoints = timepoints,
                                   # p_var = 'productive_frequency',
                                   p_var = 'normalized_frequency',
                                   compartment = 'tumor') {
  read_adaptive_seqs(force_reload = F)
  pat_arr <- unique(arr[patient %in% parent.frame(3)$patient &
                        get(timepoint_v) %in% allowed_timepoints,
                        .(adaptive_sample_name, amino_acid, productive_frequency,
                          normalized_frequency, timepoint)])
  if (pat_arr[, uniqueN(timepoint)] == 1) return(NULL)
  if (null_dat(pat_arr)) return(NULL)
  pat_arr[, productive_frequency := sum(productive_frequency, na.rm = T),
          by = .(amino_acid, timepoint)]
  pat_arr[, normalized_frequency := sum(normalized_frequency, na.rm = T),
          by = .(amino_acid, timepoint)]
  pat_arr <- unique(pat_arr)
  # pat_arr[amino_acid %in% 'CSVQGAGTEAFF']
  # pat_arr[amino_acid %in% 'CSVPDPLGNTEAFF']

  setkey(pat_arr, amino_acid, timepoint)
  subs <- expand.grid('amino_acid' = pat_arr[, unique(amino_acid)],
                      'timepoint' = pat_arr[, unique(timepoint)]) %>%
    as.data.table
  # subs[, .N, by = amino_acid][N != length(timepoints)]
  pat_arr <- pat_arr[subs, ]

  # pat_arr[amino_acid %in% pat_arr[, .N, by = amino_acid][N != 3, amino_acid]] %>%
  #   { .[order(amino_acid, timepoint)] } %>%
  #   { .[, .N == 0] } %>%
  #   stopifnot()

  pat_arr[, timepoint := factor(timepoint, levels = timepoints)]
  pat_arr[is.na(normalized_frequency), normalized_frequency := 0]
  pat_arr[is.na(productive_frequency), productive_frequency := 0]
  pat_arr[, 'value' := get(p_var)]
  pat_arr[, 'shared_timepoints' := sum(value > 0), by = amino_acid]
  pat_arr[, shared_timepoints := factor(shared_timepoints)]
  wide_dat <- dcast(pat_arr, amino_acid ~ timepoint, value.var = 'value')
  timepoints_present <- intersect(timepoints, pat_arr[, unique(timepoint)])
  wide_dat[, 'clone_N' := .N, by = timepoints_present]
  cluster_assigns <- unique(wide_dat, by = timepoints_present) %>%
    { .[, 'cluster_ID' := 1:.N] }
  wide_dat <-
    controlled_merge(wide_dat,
                     cluster_assigns[, c('cluster_ID', timepoints_present),
                                     with = F])
  pat_arr <- controlled_merge(pat_arr,
                              wide_dat[, .(amino_acid, cluster_ID, clone_N)])
  pat_arr <- pat_arr[shared_timepoints != 0]
  size_breaks <- setNames(unlist(pat_arr[, .(min(clone_N), median(clone_N),
                                    max(clone_N))]), NULL)
  pat_arr <- unique(pat_arr, by = c('cluster_ID', 'timepoint'))
  cols <- gen_color_vector(n = 3, name = 'Zissou1') %>%
    { .[1:pat_arr[, uniqueN(shared_timepoints)]] }

  p <- plot_parallel_coords(pat_arr,
    facet_var = facet_var,
    size_var = 'clone_N',
    timepoint_v = timepoint_v,
    point_alpha = .5,
    line_alpha = .1,
    swarm_width = .3,
    group_var = 'amino_acid',
    man_colors = cols,
    colour_var = colour_var) +
    scale_y_continuous(name = var_to_label(p_var), trans = 'log10') +
    scale_size_continuous(name = '# Clones', breaks = size_breaks) +
    ggtitle(sprintf('%s', compartment))
  return(p)
}
