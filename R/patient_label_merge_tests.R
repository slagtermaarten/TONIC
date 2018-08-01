cond_rm('blood_adaptive')
cond_rm('patient_labels')

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
    filter(timepoint == 'On nivo' & is.na(clinical_response)) %>%
    { .[, .N > 0] }
  if (bool) {
    print('prob 8')
    filter_patients(patient_labels, 'clinical_response') %>%
      filter(timepoint == 'On nivo' & is.na(clinical_response))
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
}


