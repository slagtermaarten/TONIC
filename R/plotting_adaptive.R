#' Cross data type correlation matrix
#'
#'
plot_datatype_cor_mat <- function(tp = 'Baseline',
                         cormethod = 'spearman',
                         sig_level = c(.05, .0005, 000005)) {
  library(corrplot)

  patient_labels <- controlled_merge(patient_labels, seq_stat_overview)
  cibersort <- read_cibersort()
  cibersort[, rmse := NULL]
  cibersort[, p_value := NULL]
  cibersort_celltypes <- setdiff(colnames(cibersort), c('patient', 'timepoint'))
  t_cell_types <- setdiff(grep('^t_cells', cibersort_celltypes, value = T),
                          't_cells_gamma_delta')
  setnames(cibersort, cibersort_celltypes,
           sprintf('cibersort_%s', cibersort_celltypes))
  cibersort$cibersort_total_t_cells <-
    apply(cibersort[, sprintf('cibersort_%s', t_cell_types), with = F], 1, sum)
  patient_labels <- controlled_merge(patient_labels, cibersort,
                                     by_cols = c('patient', 'timepoint'))

  p_dat <-
    patient_labels[, .(timepoint,
                       tis_score,
                       productive_clonality,
                       adaptive_t_cells,
                       sample_amount_ng,
                       cd8_mm2,
                       cibersort_t_cells_cd8,
                       cibersort_total_t_cells,
                       # cibersort_t_cells_cd4_naive,
                       # cibersort_t_cells_cd4_memory_resting,
                       # cibersort_t_cells_cd4_memory_activated,
                       # cibersort_t_cells_follicular_helper,
                       # cibersort_t_cells_regulatory,
                       # cibersort_neutrophils,
                       pd_l1_tumor,
                       pd_l1_immunoinfiltrate,
                       ca15_3,
                       ldh,
                       crp,
                       s_til,
                       mut_load,
                       DD_rank,
                       gen_SCNA_score,
                       weighted_gen_SCNA_score,
                       perc_LOH,
                       weighted_perc_LOH,
                       ichao1,
                       observed_richness,
                       efron_thisted_estimator,
                       daley_smith_estimator)] %>%
    as.data.frame

  p_dat <- p_dat %>%
    dplyr::filter(timepoint %in% tp) %>%
    dplyr::select(-timepoint) %>%
    # { .[, apply(., 2, function(x) all(!is.na(x))), ] } %>%
    { .[apply(., 1, function(x) all(!is.na(x))), ] } %>%
    as.matrix

  print(dim(p_dat))

  ## Format column names a bit
  colnames(p_dat) <- gsub('^s_', 'stromal ', colnames(p_dat))
  colnames(p_dat) <- gsub('pd_l1', 'PD-L1', colnames(p_dat))
  colnames(p_dat) <- gsub('ca15_3', 'ca15.3', colnames(p_dat))
  colnames(p_dat) <- gsub('mut', 'mutational', colnames(p_dat))
  colnames(p_dat) <- gsub('t_cell', 'T cell', colnames(p_dat))
  colnames(p_dat) <- gsub('gen', 'gen.', colnames(p_dat))
  colnames(p_dat) <- tonic_cap(gsub('_', ' ', colnames(p_dat)),
                               cap_first_word_only = T)

  p_mat <- cor.mtest(p_dat, method = cormethod)$p
  p_mat <- matrix(p.adjust(p_mat, method = 'fdr'),
                  nrow = nrow(p_mat), byrow = T)

  corrplot(cor(p_dat, method = cormethod),
           p.mat = p_mat,
           addrect = 6,
           tl.col = 'gray40',
           pch.col = 'white',
           pch.cex = .6,
           # order = 'AOE',
           order = 'hclust',
           # tl.pos = 'ld',
           # type = 'lower',
           tl.cex = .9,
           sig.level = sig_level,
           insig = 'label_sig',
           method = 'color',
           mar = c(0, 0, 6, 0), # http://stackoverflow.com/a/14754408/54964
           title = sprintf('%s - %s',
                           paste(tp, collapse = ', '),
                           tonic_cap(cormethod, cap_first_word_only = T)))
}


plot_parallel_adaptive <- function(p_var = 'efron_thisted_estimator',
                                   facet_var = 'arm', compartment = 'tumor',
                                   colour_var = 'clinical_response',
                                   x_axis_tps = c(timepoints, blood_timepoints),
                                   ...) {
  compartment <- match.arg(compartment, choices = c('blood', 'tumor'))
  if (compartment == 'tumor') {
    timepoint_v <- 'timepoint'
    p_dat <- copy(patient_labels[, c('patient', timepoint_v, 'clinical_response',
                                 'arm', p_var), with = F])
  } else if (compartment == 'blood') {
    timepoint_v <- 'blood_timepoint'
    p_dat <- copy(blood_adaptive[, c('patient', timepoint_v, 'clinical_response',
                                 'arm', p_var), with = F])
  }
  setnames(p_dat, p_var, 'value')
  p_dat <- p_dat[!is.na(value)] %>%
    { .[get(timepoint_v) %in% x_axis_tps] }
  p_dat <- filter_patients(p_dat, colour_var, facet_var)
  p_dat[, (timepoint_v) := droplevels(get(timepoint_v))]

  sum_dat <- p_dat %>%
    { .[, .('value' = median(value, na.rm = T),
            'patient' = 'median',
            'clinical_response' = NA),
        by = c(timepoint_v, facet_var)] }

  plot_parallel_coords(p_dat, facet_var = facet_var,
                       group_var = 'patient',
                       timepoint_v = timepoint_v,
                       # sum_dat = sum_dat,
                       colour_var = colour_var) +
    scale_y_continuous(name = var_to_label(p_var, label_reps)) +
    ggtitle(sprintf('%s', compartment))
}


prepare_adaptive_FC <- function(# facet_var = NULL,
                                # facet_var = 'clinical_response',
                                dtf = patient_labels,
                                facet_var = 'arm',
                                alpha_lev = .8,
                                epsilon = 1,
                                timepoints = timepoints,
                                # var1 = 'pielou_evenness',
                                var1 = 'sample_clonality',
                                var2 = 'efron_thisted_estimator') {
  aggregate_timepoints <- function(dtf, lvar = 'efron_thisted_estimator') {
    p_dat <- Reduce(function(x, y) merge(x, y, all = TRUE, by = 'patient'),
     lapply(seq_along(timepoints), function(tpi) {
       res <- dtf[timepoint == timepoints[tpi], c('patient', lvar), with = F]
       res <- res[!is.na(get(lvar))]
       setnames(res, lvar, sprintf('%s_tp%d', lvar, tpi))
       return(res)
     }))

    res <- p_dat[, .(patient, 'FC' = 0, 'timepoint' = timepoints[1])]
    for (i in 2:length(timepoints)) {
      ## Create vars with format (var)_1(tp)
      ## Comparing all later timepoints to first timepoint
      cur_fc <- sprintf('%s_1%d', lvar, i)
      p_dat[, (cur_fc) := 
            (log2(get(sprintf('%s_tp%d', lvar, i)) + epsilon) - 
             log2(get(sprintf('%s_tp%d', lvar, 1)) + epsilon))]
      res <- rbind(res, 
                   p_dat[, .(patient, 'FC' = get(cur_fc), 
                             'timepoint' = timepoints[i])])
    }
    res <- unique(res[!is.na(FC)])
    setnames(res, 'FC', lvar)
    return(res)
  }

  ## Merge the FCs of the two variables
  agg_tp <-
    merge(aggregate_timepoints(dtf = dtf, lvar = var1),
          aggregate_timepoints(dtf = dtf, lvar = var2),
          by = c('patient', 'timepoint'), all = T) %>%
    { .[naturalsort::naturalorder(patient)] }

  ## Merge in patient labels
  agg_tp <- patient_labels[, .(clinical_response, patient, arm)] %>%
    { unique(., by = 'patient') } %>%
    { merge(agg_tp, ., by = c('patient'), all.y = T, all.x = T) }

  agg_tp <- filter_patients(agg_tp, facet_var)

  ## Restore this variable as this may have been incorrectly merged
  agg_tp[, comb_time_resp := sprintf('%s-%s', timepoint, clinical_response)]
  missing_levs <- agg_tp[grepl('NA', comb_time_resp), comb_time_resp]
  messagef('Omitting records with %s', paste(missing_levs, collapse = ', '))
  agg_tp <- agg_tp[comb_time_resp %nin% missing_levs]

  ## Do this code block if we're dealing with tumor data only
  if ('timepoint' %in% colnames(dtf) &&
      length(timepoints) == 3 &&
      all(timepoints == patient_labels[, levels(timepoint)])) {
    agg_tp[, comb_time_resp := factor(comb_time_resp,
                  levels = patient_labels[, levels(comb_time_resp)])]
    ## Patients with On nivo but no Post-induction, set post-induction to 0
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
  }

  agg_tp <- agg_tp[order(patient, timepoint)]
  agg_tp <- agg_tp[order(timepoint)]
  return(agg_tp)
}


plot_FC <- function(agg_tp,
                    facet_var = 'arm',
                    alpha_lev = .8,
                    colour_var = 'comb_time_resp',
                    # plot_points = 'nonzero',
                    plot_points = 'all',
                    plot_counts = T,
                    var1 = 'sample_clonality',
                    var2 = 'efron_thisted_estimator') {
  if (is.null(agg_tp)) return(NULL)
  agg_tp <- filter_patients(agg_tp, facet_var, colour_var)
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
    s_plot <- s_plot + geom_point() +
      geom_path(alpha = alpha_lev, show.legend = F)
  }

  if (!is.null(colour_var) && colour_var == 'comb_time_resp') {
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

  if (!is.null(colour_var) && colour_var == 'timepoint') {
    s_plot <- s_plot + aes_string(colour = colour_var) +
      # scale_fill_manual(values = timepoint_colors) +
      scale_colour_manual(name = '',
                          values = setNames(as.character(timepoint_colors),
                                            names(timepoint_colors))) +
      scale_shape_manual(name = '',
               values = setNames(c(1, 2, 3) + 14,
                                 names(timepoint_colors)))
  }

  s_plot <- s_plot +
    scale_x_continuous(sprintf('%s [log2 FC from Baseline]',
                               var_to_label(var1, label_reps)),
                       trans = 'identity',
                       expand = c(0.1, 0.1)) +
    scale_y_continuous(sprintf('%s [log2 FC from Baseline]',
                               var_to_label(var2, label_reps)),
                       trans = 'identity',
                       expand = c(0.1, 0.1)) +
    theme(legend.position = 'top', legend.direction = 'vertical') +
    gg_legend_alpha_cancel

  if (!is.null(facet_var)) {
    s_plot <- s_plot +
      facet_wrap(as.formula(sprintf('~ %s',
                                    paste(facet_var, collapse = ' + '))),
                 nrow = length(facet_var))
  }

  if (plot_counts) {
    counts <- agg_tp[, .SD[timepoint == dplyr::last(timepoint)], by = patient]
    if (!is.null(facet_var)) {
      counts <- counts[, .N, by = .('v1' = sign(get(var1)), 
                                    'v2' = sign(get(var2)), 
                                    'fv' = get(facet_var))]
      counts <- counts[!(v1 == 0 & v2 == 0)]
      counts[, 'perc' := N / sum(N), by = fv]
      setnames(counts, 'fv', facet_var)
      fv_u <- counts[, uniqueN(get(facet_var))]
    } else {
      ## Untested
      counts <- counts[, .N, by = .('v1' = sign(get(var1)), 
                                    'v2' = sign(get(var2)))] 
      counts <- counts[!(v1 == 0 & v2 == 0)]
      counts[, 'perc' := N / sum(N), by = fv]
      fv_u <- 1
    }
    counts[, 'total' := sum(N), by = facet_var]
    degree <- .95
    counts[, v1 := ifelse(v1 == -1, 1-degree, degree)]
    counts[, v2 := ifelse(v2 == -1, 1-degree^(1/fv_u), degree^(1/fv_u))]
    counts[, (var1) := interpolate_in_gg_range(s_plot, 'x', degree = v1)]
    counts[, (var2) := interpolate_in_gg_range(s_plot, 'y', degree = v2)]
    # counts[, label := sprintf('%d/%d (%s)', N, total, scales::percent(perc))]
    # counts[, label := sprintf('%d/%d', N, total)]
    counts[, label := scales::percent(perc)]
    s_plot <- s_plot + geom_text(data = counts, 
                       aes_string(group = NULL, colour = NULL,
                                  x = var1, y = var2, label = 'label'),
                       size = 3, guide = F)
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

  l_patient <- patient
  suf_data <- all(c(tp1, tp2, tp3, tp4) %in%
                  arr[patient == l_patient, as.character(timepoint)])

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

  FCs_l <- arr[patient == l_patient] %>%
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
  l_timepoint <- timepoint
  if (timepoint %in% blood_timepoints) {
    s_dat <- blood_adaptive
    s_dat <- s_dat[as.character(blood_timepoint) == l_timepoint]
  } else if (timepoint %in% timepoints) {
    s_dat <- patient_labels
    s_dat <- s_dat[as.character(timepoint) == l_timepoint]
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
  p_dat <- filter_patients(p_dat, y_var, facet_var)
  if (null_dat(p_dat)) return(NULL)
  comp_levels <- p_dat[, levels(get(y_var))]
  if (is.null(comp_levels)) comp_levels <- p_dat[, unique(get(y_var))]

  res <- rbindlist(lapply(measures, function(measure) {
    if (measure %nin% colnames(p_dat)) {
      warningf('%s not found', measure)
      return(NULL)
    }

    dtf <- p_dat[, .('measure' = var_to_label(measure, label_reps),
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
  s_dat <- filter_patients(s_dat, 'clinical_response')
  single_plots <- lapply(measures, function(measure) {
    y_label <-
      sprintf('%s - %s', var_to_label(measure, label_reps), l_timepoint)
    p <- ggplot(s_dat,
                aes_string(x = 'clinical_response', y = measure,
                           fill = 'clinical_response')) +
      geom_boxplot() +
      ggbeeswarm::geom_quasirandom() +
      ggpubr::stat_compare_means(label = 'p.signif',
                                 label.y = s_dat[, 1.2 *
                                                 max(get(measure), na.rm= T)]) +
      scale_fill_manual(name = '', values = resp_colors, guide = F) +
      scale_y_continuous(name = y_label) +
      xlab('Clinical response')
    if (!is.null(facet_var)) {
      p <- p + facet_grid(as.formula(sprintf('~ %s ', facet_var)))
    }
    return(p)
  })

  # ACC_plots <- lapply(seq(1, length(measures) - 1),
  #                     function(i) lapply(seq(i + 1, length(measures)),
  #                                        function(j) c(measures[i], measures[j]))) %>%
  #   unlist(recursive = F) %>%
  #   { lapply(., function(x) {
  #       x_label <- sprintf('%s - %s', var_to_label(x[1], label_reps),
  #                          l_timepoint)
  #       y_label <- sprintf('%s - %s', var_to_label(x[2], label_reps),
  #                          l_timepoint)
  #       s_dat <- s_dat[!is.na(get(x[1])) & get(x[1]) != 'NA']
  #       s_dat <- s_dat[!is.na(get(x[2])) & get(x[2]) != 'NA']
  #       s_dat[, 'comb_rank' := frank(frank(get(x[1]), ties.method = 'first') + 
  #                                    frank(get(x[2]), ties.method = 'first'), 
  #                                    ties.method = 'first')]
  #       sscurves <- evalmod(scores = s_dat[, c(x, 'comb_rank'), with = F], 
  #                           labels = s_dat[, as.integer(clinical_response == 'R')])
  #       sscurves <- evalmod(scores = s_dat[, comb_rank], 
  #                           labels = s_dat[, as.integer(clinical_response == 'R')])
  #       evalmod()
  #       s_dat
  #       browser()
  #       N_resp <- s_dat[, sum(clinical_response == 'R')]
  #       roc(s_dat[, clinical_response == 'R'], s_dat[, get(x[1])], ci = T, plot = T, col = 'red', main = l_timepoint)
  #       auc(s_dat[, clinical_response == 'R'], s_dat[, get(x[1])], ci = T)
  #       roc(s_dat[, clinical_response == 'R'], s_dat[, get(x[2])], ci = T, plot = T, col = 'blue', add = T)
  #       roc(s_dat[, clinical_response == 'R'], s_dat[, comb_rank], ci = T, plot = T, add = T)
  #       ROC(s_dat[, ], s_dat[, clinical_response == 'R'])
  #       glm(as.formula(sprintf('clinical_response ~ %s + %s', x[1], x[2])), family = 'logistic', data = s_dat)
  #
  #       ROC(s_dat[, get(x[1])], s_dat[, clinical_response == 'R'], plot = 'sp')
  #       ROC(s_dat[, get(x[2])], s_dat[, clinical_response == 'R'])
  #       ROC(s_dat[, comb_rank], s_dat[, clinical_response == 'R'])
  #       setkeyv(s_dat, x)
  #       s_dat <- s_dat[!is.na(get(x[1])) & !is.na(get(x[2]))]
  #       s_dat[order(-comb_rank), .(get(x[1]), get(x[2]), comb_rank, clinical_response)]
  #       s_dat[order(-comb_rank), 
  #         .(get(x[1]), get(x[2]), comb_rank, clinical_response)] %>%
  #         .[, .(cumsum(clinical_response == 'R'), 
  #               'acc' = cumsum(clinical_response == 'R') / cumsum(!is.na(clinical_response)))]
  #
  #   }) }

  # pacman::p_load(pso)

  # dual_plots[[2]]
  dual_plots <- lapply(seq(1, length(measures) - 1),
                      function(i) lapply(seq(i + 1, length(measures)),
                                         function(j) c(measures[i], measures[j]))) %>%
    unlist(recursive = F) %>%
    { lapply(., function(x) {
        x_label <- sprintf('%s - %s', var_to_label(x[1], label_reps),
                           l_timepoint)
        y_label <- sprintf('%s - %s', var_to_label(x[2], label_reps),
                           l_timepoint)
        p <- ggplot(data = s_dat, aes_string(x = x[1], y = x[2],
                    colour = 'clinical_response')) +
          geom_point() +
          scale_colour_manual(name = '', values = resp_colors, guide = F) +
          labs(list(x = x_label, y = y_label))
        if (!is.null(facet_var)) {
          p <- p + facet_grid(as.formula(sprintf('~ %s ', facet_var)))
        } else {
          ## Initial thresholds
          init_ts <-
            s_dat[, .(quantile(get(x[1]), probs = .5, na.rm = T),
                      quantile(get(x[2]), probs = .5, na.rm = T))] %>%
            unlist
          to_optimize <- function(ts) {
            quad_N <- c(
              s_dat[get(x[1]) >= ts[1] & get(x[2]) >= ts[2], .N],
              s_dat[get(x[1]) >= ts[1] & get(x[2]) < ts[2], .N],
              s_dat[get(x[1]) < ts[1] & get(x[2]) >= ts[2], .N],
              s_dat[get(x[1]) < ts[1] & get(x[2]) < ts[2], .N]) + 1
            # if (any(quad_N == 0)) return(1)
            quad_success <- c(
              s_dat[get(x[1]) >= ts[1] & get(x[2]) >= ts[2],
                    sum(clinical_response == 'R')],
              s_dat[get(x[1]) >= ts[1] & get(x[2]) < ts[2],
                    sum(clinical_response == 'R')],
              s_dat[get(x[1]) < ts[1] & get(x[2]) >= ts[2],
                    sum(clinical_response == 'R')],
              s_dat[get(x[1]) < ts[1] & get(x[2]) < ts[2],
                    sum(clinical_response == 'R')])
            suppressWarnings(prop.test(quad_success, quad_N)$p.value)
          }
          optimal_thresholds <- pso::psoptim(par = init_ts, fn = to_optimize,
            lower = unlist(s_dat[, .(quantile(get(x[1]), probs = 0.2, na.rm = T),
                                     quantile(get(x[2]), probs = 0.2, na.rm = T))]),
            upper = unlist(s_dat[, .(quantile(get(x[1]), probs = .8, na.rm = T),
                                     quantile(get(x[2]), probs = .8, na.rm = T))]))
          print((optimal_thresholds$par - init_ts) / init_ts)
          print((to_optimize(optimal_thresholds$par) - to_optimize(init_ts)) / 
                 to_optimize(init_ts))
          p <- p + geom_vline(xintercept = optimal_thresholds$par[1])
          p <- p + geom_hline(yintercept = optimal_thresholds$par[2])
          p <- p + annotate(geom = 'text',
                            label = sprintf('italic(p)==%s',
                                            fancy_scientific(optimal_thresholds$value)),
                            hjust = 1, 
                            vjust = 1, 
                            parse = T,
                            x = interpolate_in_gg_range(p, axis = 'x', 
                                                        degree = .95),
                            y = interpolate_in_gg_range(p, axis = 'y', 
                                                        degree = .95))
        }
        return(p)
    }) }
  return(c(single_plots, dual_plots))
}


#' Plot TCR abundance comparison between blood and tumor compartments
#'
#'
plot_tp_comp_direct <- function(tp1 = 'Baseline', tp2 = '-2',
                                y_var = 'normalized_frequency') {
  read_adaptive_seqs()
  plyr::llply(arr[, naturalsort::naturalsort(auto_name(unique(patient)))],
              function(patient) {
    l_patient <- patient
    suf_data <- arr[patient == l_patient,
                    all(c(tp1, tp2) %in% timepoint)]
    if (!suf_data) return(NULL)
    t_dat <- arr[patient == l_patient &
                 timepoint %in% c(tp1, tp2)]

    p_dat <- merge(t_dat[timepoint == tp1, c(y_var, 'amino_acid'), with = F],
                   t_dat[timepoint == tp2, c(y_var, 'amino_acid'), with = F],
                   all = TRUE, by = 'amino_acid')

    check_timepoint <- function(tp) {
      if (grepl('-*\\d{1,2}', tp)) return(sprintf('day %s', tp))
      else return(tp)
    }
    tp1 <- check_timepoint(tp1)
    tp2 <- check_timepoint(tp2)
    max_val <- c(p_dat[, range(get(sprintf('%s.x', y_var)), na.rm = T)],
                 p_dat[, range(get(sprintf('%s.y', y_var)), na.rm = T)]) %>%
      { max(.) }

    ggplot(p_dat, aes_string(x = sprintf('%s.x', y_var),
                             y = sprintf('%s.y', y_var))) +
      geom_count(alpha = .1) +
      scale_x_continuous(name = sprintf('%s %s',
                                        var_to_label(y_var, label_reps), tp1),
                         limits = c(0, max_val)) +
      scale_y_continuous(name = sprintf('%s %s',
                                        var_to_label(y_var, label_reps), tp2),
                         limits = c(0, max_val)) +
      geom_hline(yintercept = 0, color = 'gray20', linetype = 'dashed') +
      geom_vline(xintercept = 0, color = 'gray20', linetype = 'dashed') +
      ggtitle(sprintf('%s - %s - %s',
                      patient,
                      patient_labels[patient == l_patient,
                                     unique(arm)],
                      patient_labels[patient == l_patient,
                                     unique(clinical_response)])) +
      theme(legend.position = 'right', legend.direction = 'vertical',
            aspect.ratio = 1)
  })
}


format_patient_title <- function(patient) {
  # patient_labels[patient == parent.frame(3)$patient]
  l_patient <- patient
  sprintf('%s - %s - %s',
          patient,
          patient_labels[patient == l_patient,
                         unique(arm)],
          patient_labels[patient == l_patient,
                         unique(clinical_response)]) %>%
  unique()
}


#' Get filename of cached FC computation
#'
#'
get_FC_fn <- function(patient, tp1, tp2, tp3, tp4, y_var = '') {
  # list.files(file.path(rds_dir, 'FCs'))
  file.path(rds_dir, 'FCs', sprintf('%s_FCs_%s_%s_%s_%s.rds',
                                    patient, tp1, tp2, tp3, tp4,
                                    ifelse(y_var == '', '',
                                           sprintf('_%s', y_var))))
}


#' Compare TCR abundance FCs in two compartments
#'
#'
plot_tp_comp_FCs <- function(tp1 = 'Baseline', tp2 = 'Post-induction',
                             tp3 = '-2', tp4 = '0',
                             y_var = 'normalized_frequency',
                             tcr_subset = NULL) {
  plyr::llply(arr[, naturalsort::naturalsort(auto_name(unique(patient)))],
              function(l_patient) {
    fn <- get_FC_fn(l_patient, tp1, tp2, tp3, tp4, y_var)
    if (file.exists(fn)) {
      FCs <- readRDS(fn)
    } else {
      return(NULL)
    }

    ## 2018-06-18 14:52 added it turns out this is not what I need as no TCRs
    ## are filtered out.
    if (!is.null(tcr_subset)) {
      if (tcr_subset == 'til_baseline') {
        FCs <-
          controlled_merge(FCs,
                           arr[timepoint == 'Baseline' &
                               patient == l_patient,
                               .(timepoint, patient, normalized_frequency)])
        FCs <- FCs[normalized_frequency > 0]
      }
    }

    ggplot(FCs, aes(x = exp1, y = exp2)) + geom_count(alpha = .1) +
      xlab(sprintf('Expansion in tumor (%s vs. %s)', tp1, tp2)) +
      ylab(sprintf('Expansion in blood (day %s vs. day %s)', tp3, tp4)) +
      geom_hline(yintercept = 0, color = 'gray20', linetype = 'dashed') +
      geom_vline(xintercept = 0, color = 'gray20', linetype = 'dashed') +
      # ggtitle(sprintf('%s %s-%s vs %s-%s', l_patient, tp1, tp2, tp3, tp4)) +
      ggtitle(format_patient_title(l_patient)) +
      theme(legend.position = 'right', legend.direction = 'vertical',
            aspect.ratio = 1)
  })
}


compute_x_label <- function(t_dat, x_var, facet_var) {
  t_dat <- filter_patients(t_dat, x_var, facet_var)
  t_dat[!is.na(value),
        'label' := sprintf('%s~(italic(n)==%d)',
                           gsub(' ', '~', get(x_var)), .N),
        by = c(facet_var, x_var)]
  if (x_var == 'arm') {
    labels <- t_dat[, unique(label)]
    labels_stripped <- t_dat[, gsub('~', ' ',
                                    gsub('(.*)~\\(.*', '\\1', unique(label)))]
    right_order <- labels[match(t_dat[, levels(arm)], labels_stripped)]
    t_dat[, label := factor(label, levels = right_order)]
  }
  return(t_dat)
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

  t_dat <- controlled_merge(t_dat[!is.na(value)],
                            patient_labels[, .(patient, arm,
                                               clinical_response)]) %>%
    unique(by = c('patient'))

  t_dat <- filter_patients(t_dat, x_var, colour_var, facet_var)
  t_dat <- compute_x_label(t_dat, x_var, facet_var)
  return(t_dat)
}

plot_adaptive_summary_stats <- function(tp1 = 'Baseline',
                                        tp2 = 'On nivo',
                                        x_var = 'clinical_response',
                                        colour_var = 'clinical_response',
                                        facet_var = 'arm',
                                        comp_measure = 'sample_clonality') {

  t_dat <- compare_adaptive_summary_stats(tp1 = tp1, tp2 = tp2,
                                          x_var = x_var,
                                          colour_var = colour_var,
                                          facet_var = facet_var,
                                          comp_measure = comp_measure)
  label_proc <- function(labels) {
    sapply(labels, function(lab) parse(text = lab))
  }

  p1 <- ggplot(t_dat[!is.na(label)],
               aes_string(x = 'label', y = 'value', fill = colour_var))
  p1 <- p1 + geom_boxplot()
  p1 <- p1 + ggbeeswarm::geom_quasirandom()
  if (x_var == 'arm') {
    if (F) {
      p1 <- p1 + ggpubr::stat_compare_means(label.y = t_dat[, 1.2 * max(value)])
      p1 <- p1 + ggpubr::stat_compare_means(aes(label = ..p.signif..),
                                    # method = "t.test",
                                    ref.group = t_dat[arm == 'No induction',
                                                      unique(label)])
    } else {
      p1 <- p1 + ggpubr::stat_compare_means(label.y = t_dat[, 1.5 * max(value)])
      pval <- ggpubr::compare_means(formula = value ~ arm, t_dat,
                                    method = "kruskal.test")$p
      if (pval <= 0.05) {
        my_comparisons <-
          lapply(t_dat[, setdiff(levels(arm), 'No induction')],
                 function(x) t_dat[arm %in% c(x, 'No induction'),
                                   as.character(unique(label))])
        p1 <- p1 + ggpubr::stat_compare_means(comparisons = my_comparisons,
                     label.y = t_dat[, (1 + seq_along(my_comparisons) * .1) *
                                     max(value)],
                                              label = 'p.signif')
      }
    }
  } else if (x_var == 'clinical_response') {
    p1 <- p1 + ggpubr::stat_compare_means(label.y = t_dat[, 1.5 * max(value)])
    p1 <- p1 + ggpubr::stat_compare_means(label.y = t_dat[, 1.2 * max(value)],
                                          comparisons = t_dat[, unique(label)])
  }
  p1 <- p1 + scale_fill_manual(name = '',
                               values = tonic_color_palettes[[colour_var]])
  p1 <- p1 + scale_colour_manual(name = '',
                                 values = tonic_color_palettes[[colour_var]])
  p1 <- p1 + scale_y_continuous(name = sprintf('%s %s vs. %s',
                                               var_to_label(comp_measure,
                                                            label_reps),
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
                             plot_counts = T,
                             var1 = 'sample_clonality',
                             var2 = 'efron_thisted_estimator') {
  agg_tp <- prepare_adaptive_FC(dtf = dtf, 
                                facet_var = facet_var,
                                timepoints = p_timepoints,
                                alpha_lev = alpha_lev, 
                                var1 = var1,
                                var2 = var2)
  sr(comb_time_resp_palette)
  agg_tp_f <- agg_tp[timepoint %in% p_timepoints & !is.na(clinical_response)]
  # agg_tp_f[, timepoint]
  # agg_tp[, levels(comb_time_resp)]
  # agg_tp_f[, levels(comb_time_resp)]
  plot_FC(agg_tp_f,
          facet_var = facet_var, alpha_lev = alpha_lev,
          plot_counts = plot_counts,
          colour_var = colour_var, var1 = var1, var2 = var2)
}


#' Plot all relevant comps for Adaptive summary statistics
#'
#'
plot_all_comps <- function(f = plot_adaptive_summary_stats,
                           colour_var = 'arm',
                           x_var = 'arm',
                           facet_var = NULL) {
  list(
  'sample_clonality_BL_vs_PI' =
    f(tp1 = 'Baseline', tp2 = 'Post-induction',
    comp_measure = 'sample_clonality',
    colour_var = colour_var,
    x_var = x_var,
    facet_var = facet_var),
  'sample_clonality_BL_vs_ON' =
  f(tp1 = 'Baseline', tp2 = 'On nivo',
    comp_measure = 'sample_clonality',
    colour_var = colour_var,
    x_var = x_var,
    facet_var = facet_var),
  'Adaptive_T_cells_BL_vs_PI' =
  f(tp1 = 'Baseline', tp2 = 'Post-induction',
    comp_measure = 'adaptive_t_cells',
    colour_var = colour_var,
    x_var = x_var,
    facet_var = facet_var),
  'Adaptive_T_cells_BL_vs_ON' =
  f(tp1 = 'Baseline', tp2 = 'On nivo',
    comp_measure = 'adaptive_t_cells',
    colour_var = colour_var,
    x_var = x_var,
    facet_var = facet_var),
  'Rep_size_BL_vs_PI' =
  f(tp1 = 'Baseline', tp2 = 'Post-induction',
    comp_measure = 'efron_thisted_estimator',
    colour_var = colour_var,
    x_var = x_var,
    facet_var = facet_var),
  'Rep_size_BL_vs_ON' =
  f(tp1 = 'Baseline', tp2 = 'On nivo',
    comp_measure = 'efron_thisted_estimator',
    colour_var = colour_var,
    x_var = x_var,
    facet_var = facet_var))
}


plot_TCR_chronological <- function(patient = 'pat_11',
                                   timepoint_v = 'timepoint',
                                   facet_var = NULL,
                                   p_var = 'normalized_frequency',
                                   colour_var = 'shared_timepoints',
                                   x_axis_tps = timepoints,
                                   grep_it_timepoints = NULL,
      highlight_function = function(arr) {
        arr[it_timepoint == 'On nivo' & freq_rank >= .8, hl_var := T]
      },
                                   compartment = 'tumor') {
  arr <- prepare_TCR_chrono(patient = patient, timepoint_v = timepoint_v,
                            facet_var = facet_var, colour_var = colour_var,
                            highlight_function = highlight_function,
                            x_axis_tps = x_axis_tps,
                            p_var = p_var, compartment = compartment)
  if (null_dat(arr)) return(NULL)

  size_breaks <- 
    arr[, .(min(cluster_N), ceiling(mean(range(cluster_N))), 
            max(cluster_N))] %>%
    unlist %>%
    setNames(NULL)

  if (colour_var == 'shared_timepoints') {
    cols <- gen_color_vector(n = 3, name = 'Zissou1') %>%
      darken(c(1, 1.15, 1)) %>%
      { .[1:arr[, uniqueN(shared_timepoints)]] } %>%
      setNames(NULL) %>%
      attr_pass('class', 'color_vector')
  } else if (colour_var == 'it_timepoints') {
    cols <- gen_color_vector(n = arr[, uniqueN(it_timepoints)],
                             name = 'Zissou1') %>%
      setNames(NULL) %>%
      attr_pass('class', 'color_vector')
  } else if (colour_var == 'hl_var') {
    cols <- gen_color_vector(n = 2, name = 'Zissou1') %>%
      setNames(NULL) %>%
      attr_pass('class', 'color_vector')
  }

  if (!is.null(grep_it_timepoints)) {
    arr <- arr[grepl(grep_it_timepoints, it_timepoints)]
  }

  p <- plot_parallel_coords(arr,
                            facet_var = facet_var,
                            filter_vals = F,
                            size_var = 'cluster_N',
                            timepoint_v = timepoint_v,
                            point_alpha = .5,
                            line_alpha = .1,
                            swarm_width = .3,
                            group_var = 'amino_acid',
                            man_colors = cols,
                            colour_var = colour_var) +
    scale_y_continuous(name = var_to_label(p_var, label_reps),
                       trans = 'log10') +
    scale_size_continuous(name = '# Clones', breaks = size_breaks,
                          trans = 'identity', range = c(.5, 8)) +
    ggtitle(sprintf('%s', compartment)) +
    guides(
     size = guide_legend(order = 1),
     colour = guide_legend(order = 2),
     fill = guide_legend(order = 2)
    )
  return(p)
}
