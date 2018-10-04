cond_rm('blood_adaptive')
cond_rm('patient_labels')


#' Filter out patients that have dubious annotation
#'
#'
filter_patients <- name <- function(p_dat, ...) { 
  comb_vars <- as.character(...)
  clinical_params <- c('response', 'clinical_response', 'comb_time_resp')
  if (is.null(comb_vars)) return(p_dat)
  if (any(comb_vars %in% clinical_params)) {
    # resp_var <- comb_vars[which(comb_vars %in% clinical_params)]
    p_dat <- p_dat[patient %nin% c('pat_63', 'pat_64')]
    for (varn in clinical_params) {
      if (varn %in% colnames(p_dat)) {
        N <- p_dat[is.na(get(varn)), .N]
        if (N > 0) {
          messagef('Removing %d donors due to absence of %s', N, varn)
          p_dat <- p_dat[!is.na(get(varn))]
        }
      }
    }
  }
  if (F) {
    pres_cols <- 
      intersect(colnames(p_dat), c('patient', 'arm', 'timepoint', 
                                   'blood_timepoint',
                                   'filename', 'adaptive_sample_name'))
    p_dat <- unique(p_dat, by = pres_cols)
  }
  return(p_dat)
}


merge_tests <- function(idx = 1) {
  if (patient_labels[, uniqueN(arm) > 1, by = patient][, any(V1)]) {
    print('prob 1')
    browser()
  }

  if (patient_labels[patient %in% c('pat_64', 64), .N < 3]) {
    print('prob 1.5')
    browser()
  }

  if (patient_labels[arm == 'Cyclophosphamide', .N] == 0) {
    print('prob 2')
    print(idx)
    browser()
  }

  if (exists('blood_adaptive') && 'arm' %in% colnames(blood_adaptive) &&
      blood_adaptive[is.na(arm), .N] > 0) {
    print('prob 3')
    blood_adaptive[is.na(arm)]
    browser()
  }

  if (patient_labels[patient == 'pat_17' & is.na(arm), .N] > 0) {
    patient_labels[patient == 'pat_17']
    patient_labels[patient == 'pat_17', .N]
    print('prob 4')
    browser()
  }

  if (patient_labels[patient == 'pat_NA', .N > 0]) {
    patient_labels[patient == 'pat_NA']
    print('prob 5')
    browser()
  }

  if ('adaptive_sample_name' %in% colnames(patient_labels) &&
      patient_labels[grepl('B', adaptive_sample_name), .N] > 0) {
    print('prob 6')
    patient_labels[grepl('B', adaptive_sample_name)]
    browser()
  }

  if (exists('blood_adaptive') && 
      blood_adaptive[grepl('T', adaptive_sample_name), .N] > 0) {
    print('prob 7')
    browser()
  }

  bool <- filter_patients(patient_labels, 'clinical_response') %>%
    dplyr::filter(as.character(timepoint) == 'On nivo' & 
                  is.na(clinical_response)) %>%
    as.data.table %>%
    { .[, .N > 0] }
  if (length(bool) > 0 && bool) {
    print('prob 8')
    filter_patients(patient_labels, 'clinical_response') %>%
      dplyr::filter(timepoint == 'On nivo' & is.na(clinical_response))
    browser()
  }

  s_dat <- patient_labels[!is.na(filename)]
  if ('adaptive_sample_name' %in% colnames(s_dat)) {
    s_dat <- s_dat[!is.na(filename) | !is.na(adaptive_sample_name)]
  }
  if (s_dat[, any(timepoint == 'NA' | is.na(timepoint))]) {
    print('prob 9')
    # patient_labels[timepoint == 'NA' | is.na(timepoint)]
    patient_labels[is.na(timepoint)]
    browser()
  }

  if (F && any(duplicated(patient_labels[, filename]))) {
    print('prob 10: duplicatd filenames')
    browser()
  }
  
  if (patient_labels[, uniqueN(arm) > 1, by = patient][, any(V1)]) {
    browser()
    patient_labels[, uniqueN(arm) > 1, by = patient]
  }

  if (patient_labels[, uniqueN(arm) > 1, by = patient][, any(V1)]) {
    browser()
    patient_labels[, uniqueN(arm) > 1, by = patient]
  }

  if (exists('blood_adaptive') && 
      blood_adaptive[, uniqueN(arm) > 1, by = patient][, any(V1)]) {
    browser()
    blood_adaptive[, uniqueN(arm) > 1, by = patient]
  }

  if (patient_labels[, any(grepl('Pat', patient))]) {
    print('prob 11')
    browser()
  }

  if (patient_labels[, .N, by = .(patient, timepoint)] %>%
              dplyr::filter(!is.na(patient)) %>%
              dplyr::filter(N > 1) %>% nrow != 0) {
    print('prob 12: multiple entries per patient/timepoint combination')
    patient_labels[, .N, by = .(patient, timepoint)] %>%
      .[N > 1, unique(patient)] %>%
      { patient_labels[patient == .] }
    browser()
  }
  
  if (idx > 1 && 'cf_number' %in% colnames(patient_labels)) {
    sc <- patient_labels[patient == 'pat_33' & cf_number == 'CF15398', .N]
    if (sc > 0) {
      print('prob 13: incorrect sample annotation')
      browser()
    }
  }
}


