#' Wilcoxon rank sum to test for association with response at single timepoint
#'
#'
test_gene_association <- function(gene = 'CD274', timepoint = 'On nivo',
                                  patient_ids = patient_labels[, unique(patient)],
                                  y_var = 'response') {
  l_timepoint <- timepoint
  stopifnot(l_timepoint %in% timepoints)
  p_dat <- prep_gene_parallel_coords(gene = gene,
                                     colour_var = y_var)[['p_dat']]
  p_dat <- p_dat[timepoint == l_timepoint & 
                 patient %in% patient_ids]
  p_dat <- filter_patients(p_dat, y_var)
  if (p_dat[, all(is.na(value))]) return(NULL)

  test_var <- y_var
  if (y_var == 'response') {
    p_dat[response != 'SD', 'response_bin' := response %in% c('CR', 'PR')]
    test_two_levels <- p_dat[, length(table(response_bin)) != 2]
    browser(expr = test_two_levels)
    if (test_two_levels) return(NULL)
    test_var <- 'response_bin'
  }
  return(wilcox.test(as.formula(sprintf('value ~ %s', test_var)),
                     data = p_dat)$p.val)
}


#' Wilcoxon rank sum to test for association with response at single timepoint
#'
#'
test_gene_set_association <- function(gene_set = 'tis', timepoint = 'On nivo',
                                      y_var = 'clinical_response', 
                                      facet_var = NULL,
                                      patient_ids = patient_labels[, unique(patient)]) {
  stopifnot(timepoint %in% timepoints)
  p_dat <- danaher_scores.m[variable %in% gene_set]
  l_timepoint <- timepoint
  p_dat <- p_dat[as.character(timepoint) == l_timepoint]
  p_dat <- p_dat[as.character(patient) %in% patient_ids] 
  # p_dat <- filter_patients(p_dat, y_var, facet_var)
  if (null_dat(p_dat) || p_dat[, all(is.na(value))]) {
    return(NULL)
  }

  comp_levels <- p_dat[, levels(get(y_var))]
  res <- p_dat[, .('p_val' = wilcox.test(as.formula(sprintf('value ~ %s', y_var)),
                     data = .SD, exact = F)$p.val,
                   'log2FC' = median(.SD[get(y_var) == comp_levels[2], value],
                                     na.rm = T) - 
                              median(.SD[get(y_var) == comp_levels[1], value],
                                     na.rm = T)), 
              by = c('variable', facet_var)]
  res[, 'p_val.bh' := p.adjust(p_val, method = 'BH'), by = facet_var]
  res[, 'p_val.fdr' := p.adjust(p_val, method = 'fdr'), by = facet_var]
  res[, 'p_val.bonferroni' := p.adjust(p_val, method = 'bonferroni'), 
      by = facet_var]
  return(res)
}


prepare_test_gene_set_difference <- function(gene_set = sig_gene_sets[1], 
                                             facet_var = 'arm',
                                             tp1 = 'Baseline',
                                             tp2 = 'On nivo') {
  ## Select relevant data
  l_gene_set <- gene_set
  data_subs <- danaher_scores.m[variable %in% l_gene_set] %>%
    { .[timepoint %in% c(tp1, tp2)] } %>%
    { .[!is.na(value)] } %>%
    filter_patients(facet_var)
  
  ## Determine which pats are sufficiently covered for statistical testing
  allowed_pats <- data_subs %>%
    { .[, .N, by = c('patient')] } %>%
    { .[N == length(gene_set) * 2, patient] }

  stopifnot(length(allowed_pats) > 0)

  ## This is required for Ton's approach
  ## In case of multiple gene sets, reduce to median gene set score for each
  ## patient and timepoint combination
  if (length(gene_set) > 1) {
    data_subs <- 
      data_subs[, .(arm, 'variable' = 'median', 
                    'value' = median(value, na.rm = T)), 
                by = .(patient, timepoint)] %>% unique
  }
  if (null_dat(data_subs)) return(NULL)

  t_dat <- dcast(data = data_subs[patient %in% allowed_pats], 
                 formula = patient ~ timepoint,
                 value.var = 'value') %>%
    merge(unique(patient_labels[, .(patient, arm)]))

  t_dat[arm == 'Cyclofosfamide', arm := 'Cyclophosphamide']
  t_dat[, arm := factor(arm, levels = treatment_arms)]
  set_dt_types(danaher_scores.m, c('patient' = 'factor', 
                                   'timepoint' = 'factor',
                                   'variable' = 'factor'))

  return(t_dat[naturalsort::naturalorder(patient)])
}


#' Compute FC and p value of timepoint comparison
#'
#' This function was developed to be run inside of the j part of a data.table
#' object
#'
#' @param dtf .SD
compute_FC_p_val <- function(dtf, tp1 = 'Baseline', tp2 = 'On nivo') {
  if (nrow(dtf) < 2) return(NA)
  return(list('p_val' = wilcox.test(x = dtf[, get(tp1)],
                                    y = dtf[, get(tp2)],
                                    conf.int = F, paired = T)$p.value,
              'logFC' = dtf[, log2(median(get(tp2)))] - 
                        dtf[, log2(median(get(tp1)))]))
}


#' Test difference in scores between time points (signed rank test)
#'
#'
test_gene_set_difference <- function(gene_set = sig_gene_sets[1], 
                                     facet_var = 'arm',
                                     tp1 = 'Baseline',
                                     tp2 = 'On nivo') {
  t_dat <- prepare_test_gene_set_difference(gene_set = gene_set, 
                                            facet_var = facet_var, tp1 = tp1,
                                            tp2 = tp2)
  t_dat <- filter_patients(t_dat, facet_var)

  if (!is.null(facet_var)) {
    tst <- rbind(t_dat[, compute_FC_p_val(.SD, tp1 = tp1, tp2 = tp2), 
                       by = facet_var],
                 t_dat[, compute_FC_p_val(.SD, tp1 = tp1, tp2 = tp2)], 
                 fill = T)
    if (facet_var == 'arm') comb_var <- 'All arms'
    else stop('not implemented')
    tst[is.na(get(facet_var)), (facet_var) := comb_var]
  } else {
    tst <- t_dat[, compute_FC_p_val(.SD)]
  }
  tst[, 'gene_set' := gene_set]
  tst[, 'tp1' := tp1]
  tst[, 'tp2' := tp2]
  return(tst)
}


if (T) {
  t_dat <- test_gene_set_association(
    gene_set = danaher_scores.m[, levels(variable)], timepoint = 'On nivo')
  t_dat <- t_dat[log2FC > 0 & p_val <= 0.05]
  sig_gene_sets <- setdiff(t_dat[, as.character(variable)], 'mmr.loss')
}
