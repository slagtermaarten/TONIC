pacman::p_load(readxl)
pacman::p_load(naturalsort)
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

source('R/patient_label_merge_tests.R')

tonic_cleanup <- function(dtf) {
  if (null_dat(dtf)) return(NULL)
  dtf[, response := zoo::na.locf(response), by = patient]
  dtf[, clinical_response := zoo::na.locf(clinical_response), 
      by = patient]
  return(dtf)
}

## A primary source of the most relevant clinical information
patient_labels <- read.csv(file.path(p_root, 'data-raw/patient_labels.csv'),
                           dec = ',', sep = ';') %>% as.data.table %>%
  normalize_colnames()
patient_labels[arm == 'Cyclofosfamide', arm := 'Cyclophosphamide']
# patient_labels[patient == 'pat_33']
# patient_labels[patient == 'pat_17']
setnames(patient_labels, 'x', 'filename')
maartenutils::set_dt_types(patient_labels,
                           c('mean_log2_hk' = 'numeric',
                             'tis_score' = 'numeric'))
patient_labels <- patient_labels[!is.na(patient)]
patient_labels[, patient := paste0('pat_', patient)]
merge_tests(idx = 0)
# stopifnot(patient_labels[is.na(arm), .N <= 3])

if (T) {
  ## A second, more elaborate overview of clinical parameters
  clinical_annotation <- read.csv(file.path(p_root,
                           'data-raw/180222_clinical_data_all_patients.csv'),
                           dec = ',', sep = ';') %>% as.data.table %>%
    normalize_colnames()
  setnames(clinical_annotation, gsub('\\.', '_', colnames(clinical_annotation)))
  setnames(clinical_annotation, gsub('_+', '_', colnames(clinical_annotation)))
  setnames(clinical_annotation, gsub('_$', '', colnames(clinical_annotation)))
  clinical_annotation[, 'patient' := paste0('pat_', study_id)]
  clinical_annotation[, cd8_mm2 := gsub(',', '.', cd8_mm2)]
  clinical_annotation[, stil := gsub(',', '.', stil)]
  setnames(clinical_annotation, 'stil', 's_til')
  maartenutils::set_dt_types(clinical_annotation,
                             c('cd8_mm2' = 'numeric', 's_til' = 'numeric'))
  clinical_annotation[induction_therapy == 'Cyclofosfamide',
                      induction_therapy := 'Cyclophosphamide']
  clinical_annotation[, arm := factor(induction_therapy, 
                                      levels = treatment_arms)]
  clinical_annotation[, induction_therapy := NULL]
  clinical_annotation <- clinical_annotation[patient != 'pat_NA']
  clinical_annotation[, study_id := NULL]
  clinical_annotation[, response_bin := NULL]
  clinical_annotation[, timepoint := 'Baseline']
  clinical_annotation[, 'clinical_response' := 
                      ifelse(is.na(response) | response == 'NA', 'NA', 
                             ifelse(response %in% c('CR', 'PR'), 'R', 'NR'))]
  # clinical_annotation[patient %in% c('pat_50', 'pat_64'), 
  #                     .(patient, clinical_response)]
  # clinical_annotation[, response := NULL]
  # clinical_annotation[patient == 'pat_50']
  # stopifnot(length(intersect(colnames(clinical_annotation),
  #                            colnames(patient_labels))) <= 1)
  merge_tests(idx = 0.5)
  ## Only try nd fill missing values, don't add new patients.
  ## As this is and should be the last merging step, no omics data is available
  ## for patients that would be added by this step (i.e. patient 17)
  patient_labels <- controlled_merge(patient_labels, 
                                     clinical_annotation,
                                     # clinical_annotation[, .(arm, patient,
                                     #                         ca15_3,
                                     #                         cd8_mm2,
                                     #                         s_til,
                                     #                         ldh,
                                     #                         pd_l1_tumor,
                                     #                         pd_l1_immunoinfiltrate,
                                     #                         who_performance_status)], 
                                     clean_up_f = tonic_cleanup,
                                     by_cols = c('patient', 'timepoint'), all.x = T,
                                     all.y = T,
                                     dup_priority = 'a')
  patient_labels[patient %in% c('pat_50', 'pat_64'), 
                 .(patient, clinical_response, response)]
  merge_tests(idx = 0.75)
}


if (T) {
  ## 2018-05-13 15:51
  ## A third addition to the patient_labels object, placed before other merging
  ## steps in order to maximize data coverage (i.e. filling up NA fields)
  ann_corrections <-
    read_excel(file.path(data_dir, 'sample_annotation_corrections.xlsx'),
               na = c('', 'NA')) %>%
    as.data.table %>%
    maartenutils::normalize_colnames()
  ann_corrections[!is.na(response), 'clinical_response' := 
                  ifelse(response %in% c('PR', 'CR'), 'R', 'NR')]
  ann_corrections[, clinical_response := factor(clinical_response, 
                                                levels = c('NR', 'R'))]
  patient_labels <- controlled_merge(patient_labels, ann_corrections, all = T,
                                     clean_up_f = tonic_cleanup,
                                     dup_priority = 'f')
  patient_labels <- patient_labels[patient != 'pat_NA']
  merge_tests(1)
}

if (T) {
  ## Merge DNA Seq CF numbers
	dna_seq_cfs <-
		read_excel(file.path(data_dir, 'Samples_WES_TONICstageI_baseline.xlsx'),
               na = c('', 'NA')) %>%
		as.data.table %>%
		maartenutils::normalize_colnames()
  dna_seq_cfs[, 'patient' := paste0('pat_', study_id)]
  dna_seq_cfs[, c('sequencing_number', 'induction_arm', 
                  't_number', 'study_id') := list(NULL, NULL, NULL, NULL)]
  setnames(dna_seq_cfs, 'tumor/blood', 'sample_type')
  setnames(dna_seq_cfs, 'cf_number_(name_seq)', 'cf_number')
  patient_labels <- 
    controlled_merge(patient_labels,
                     dna_seq_cfs[sample_type == 'Tumor', 
                        .(cf_number, patient, timepoint = 'Baseline')],
                     all = T,
                     clean_up_f = tonic_cleanup,
                     dup_priority = 'f',
                     by_cols = c('patient', 'timepoint'))
  merge_tests(idx = 2)
}

if (F) {
  ## More recent overview of clinical response data, merge into existing overview
  response_data <- read.csv(file.path(p_root, 'data-raw/response_data.csv'),
                            dec = ',', sep = ';') %>% 
    as.data.table %>%
    normalize_colnames %>%
    ## Select last columns containing patient_id and RECIST labels
    { rev(.)[, c(2, 3)] } %>%
    { .[!is.na(patient.id)] } %>%
    setnames(., c('response', 'patient'))
    # ## Take out patients with dubious responses
    # { .[!patient %in% c(61, 63, 64)] }
  response_data[, patient := as.character(paste0('pat_', patient))]
  patient_labels[, response := NULL]
  patient_labels[, patient := as.character(patient)]
  patient_labels <- merge(patient_labels, response_data, by = 'patient')
  patient_labels[patient == 'pat_63']
}
# stopifnot(patient_labels[is.na(arm), .N <= 3])

if (T) {
	## Merge Adaptive TCR abundance measures
	## First clean up data...
	adaptive_sample_annotation <-
		read_excel(file.path(data_dir, 'SampleManifest_adaptive.xlsx'),
               na = c('', 'NA')) %>%
		as.data.table %>%
		maartenutils::normalize_colnames()
	setnames(adaptive_sample_annotation,
					 gsub('__', '_', colnames(adaptive_sample_annotation)))

  homogenize_arms <- function(vec) {
    vec <- tolower(vec)
    # vec[vec == 'cyclofosfamide'] = 'cyclophosphamide'
    vec[vec == 'rt'] = 'radiotherapy'
    vec[vec == 'cis'] = 'cisplatin'
    vec[vec == 'doxorubicin'] = 'doxo'
    vec[vec == 'doxorubicine'] = 'doxo'
    reps <- setNames(c('No induction', 'Radiotherapy', 'Cyclophosphamide',
                       'Cisplatin', 'Doxorubicin'), 
                     c('no induction', 'radiotherapy', 'cyclofosfamide',
                       'cisplatin', 'doxo'))
    return(reps[vec])
  }

  adaptive_sample_annotation[, arm := homogenize_arms(induction_arm)]
  adaptive_sample_annotation[, induction_arm := NULL]

	tumor_adaptive <- adaptive_sample_annotation[,
		c(colnames(adaptive_sample_annotation)[1:8], 'arm'), with = F]
	# patient_labels[, uniqueN(adaptive_sample_name)]
	# tumor_adaptive[, uniqueN(adaptive_sample_name)]
	tumor_adaptive <-
		tumor_adaptive[, .('adaptive_sample_name' = as.character(sample_name),
											 arm,
											 'cf_number' = cf_nummer,
											 'patient' = paste0('pat_', study_id),
											 'timepoint' =
							timepoints[as.integer(gsub('N15TON-(\\d)', '\\1', tijdspunt))])]

  ## Not all patients have all three time points assayed
  # tumor_adaptive[, .N, by = patient]
	# patient_labels[, class(patient)]
	# patient_labels[, class(timepoint)]
	# tumor_adaptive[, class(timepoint)]
	# tumor_adaptive[, class(patient)]

	## Then merge adaptive sample IDs...
  merge_tests(idx = 3)
  # tumor_adaptive[patient == 'pat_66']
  # patient_labels[patient == 'pat_66']
  ## TODO pat_33 is not merged properly
  tumor_adaptive[patient == 'pat_33']
  patient_labels[patient == 'pat_33']
	patient_labels <- controlled_merge(patient_labels, tumor_adaptive,
                                     by_cols = c('patient', 'timepoint'),
                                     dup_priority = 'f', all = T, 
                                     clean_up_f = tonic_cleanup)
  merge_tests(idx = 4)
  # stopifnot(patient_labels[is.na(arm), .N <= 3])
	if (F) {
		## Analyze missing samples
		missing_samples <- setdiff(tumor_adaptive[, adaptive_sample_name],
						patient_labels[, adaptive_sample_name])
		tumor_adaptive[adaptive_sample_name %in% missing_samples]
		missing_cfs <- tumor_adaptive[adaptive_sample_name %in% missing_samples,
																			 naturalsort(unique(cf_number))]
		missing_patients <- tumor_adaptive[adaptive_sample_name %in% missing_samples,
																			 naturalsort(unique(patient))]
		patient_labels[patient %in% missing_patients]
		patient_labels[cf_number %in% missing_cfs]

    # patient_labels[, clinical_response := zoo::na.locf(clinical_response,
    #                                                    fromLast = F),
    #                by = patient]
    # patient_labels[, clinical_response := zoo::na.locf(clinical_response,
    #                                                    fromLast = T),
    #                by = patient]
    # patient_labels[, clinical_response := zoo::na.locf(clinical_response,
    #                                                    fromLast = F),
    #                by = patient]
    # patient_labels[is.na(clinical_response)]

		## Fill missing values
		patient_labels[is.na(clinical_response), unique(patient)]
		patient_labels[is.na(arm)]
		patient_labels[patient == 'pat_33']
		pacman::p_load('zoo')
		# patient_labels[, arm := na.locf(arm), by = patient]
		patient_labels[, clinical_response := zoo::na.locf(clinical_response), 
                   by = patient]
	}

  ## Another patient_labels-like object but for blood timepoints
	blood_adaptive <- adaptive_sample_annotation[, 10:16]
	blood_adaptive <- blood_adaptive[!is.na(sample_name),
								 .('adaptive_sample_name' = sample_name,
                   'arm' = induction_arm_1, 
                   'patient' = paste0('pat_', study_id_1),
									 'blood_timepoint' = tijdspunt_in_weken)]
  blood_adaptive[, arm := factor(homogenize_arms(arm), 
                                 levels = treatment_arms)]
	blood_adaptive[, blood_timepoint := factor(blood_timepoint,
                            levels = sort(unique(blood_timepoint)))]
  merge_tests(idx = 5)
  blood_adaptive <- controlled_merge(blood_adaptive,
    unique(patient_labels, by = 'patient')[, .(patient, arm, response,
                                               clinical_response, tis_score)],
                                     dup_priority = 'f',
                                     clean_up_f = tonic_cleanup,
                                     all.x = T, by_cols = 'patient')
  merge_tests(idx = 6)

  ## ... Finally merge summary stats
  tmp <- read_tsv(file = file.path(data_dir, 'AdaptiveAllDiversity.txt')) %>%
		as.data.table %>%
		maartenutils::normalize_colnames()
	setnames(tmp, gsub('\'', '', colnames(tmp)))
	# merge(patient_labels, tmp, by = c('patient', 'timepoint'))

	setdiff(colnames(tmp), c('sample', 'efron_thisted_estimator',
													 'daley_smith_estimator', 'ichao1')) %>%
		{ setNames(rep('numeric', length(.)), .) } %>%
		set_dt_types(tmp, .)
  tmp[, daley_smith_estimator := as.numeric(daley_smith_estimator)]
	setnames(tmp, 'sample', 'adaptive_sample_name')
  setdiff(patient_labels[!is.na(adaptive_sample_name), adaptive_sample_name],
          tmp[, adaptive_sample_name]) %>%
					{ length(.) } %>%
					{ messagef('Missing summary stats for %d samples', .) }
	## Finally merge in adaptive info
	patient_labels <- controlled_merge(patient_labels, 
                                     tmp[grepl('T', adaptive_sample_name)], 
                                     dup_priority = 'f',
                                     clean_up_f = tonic_cleanup,
                                     by_cols = 'adaptive_sample_name',
													           all = T)
  merge_tests(idx = 7)
  blood_adaptive <- controlled_merge(blood_adaptive, 
                                     tmp[grepl('B', adaptive_sample_name)], 
                                     dup_priority = 'f',
                                     clean_up_f = tonic_cleanup,
                                     by_cols = 'adaptive_sample_name',
													           all = T)
  merge_tests(idx = 8)
	rm(tmp)

  tmp <- read_tsv(file = file.path(data_dir, 'AdaptiveSampleOverview.tsv')) %>%
		as.data.table %>%
		maartenutils::normalize_colnames()
  tmp[, total_t_cells := NULL]

  ## This is Fraction T/B Cells of Nucleated Cells Estimate
  setnames(tmp, 'fraction_productive_of_cells_mass_estimate',
           'adaptive_t_cells')
	setnames(tmp, 'sample_name', 'adaptive_sample_name')
	setnames(tmp, gsub('\'', '', colnames(tmp)))
	patient_labels <- controlled_merge(patient_labels, 
                                     tmp[grepl('T', adaptive_sample_name)], 
                                     dup_priority = 'f',
                                     clean_up_f = tonic_cleanup,
                                     by_cols = 'adaptive_sample_name',
													           all = T)
  merge_tests(idx = 9)
  blood_adaptive <- controlled_merge(blood_adaptive, 
                                     tmp[grepl('B', adaptive_sample_name)], 
                                     clean_up_f = tonic_cleanup,
                                     dup_priority = 'f',
                                     by_cols = 'adaptive_sample_name',
													           all = T)
  merge_tests(idx = 10)
	setkey(patient_labels, patient, timepoint)
}

## Final formatting
## One laesion responds, the other one doesn't
## Leave out patient
manual_clinical_corrections <- function(dtf) {
  # dtf <- dtf[patient != 'pat_64']
  dtf <- dtf[patient == 'pat_33', clinical_response := 'R']
  dtf <- dtf[patient == 'pat_33', response := 'PR']
  dtf[, response := factor(response, levels = c('PD', 'SD', 'PR', 'CR'))]
  dtf[arm == 'Cyclofosfamide',  arm := 'Cyclophosphamide']
  dtf[, arm := factor(arm, levels = treatment_arms)]
  dtf[patient == 'pat_64', clinical_response := NA]
  dtf[patient == 'pat_64', response := NA]
  # patient_labels[patient == 'pat_64', .(patient, response, clinical_response)]
  return(dtf)
}
merge_tests(idx = 10.1)
patient_labels <- manual_clinical_corrections(patient_labels)
patient_labels[, 'timepoint_number' := match(timepoint, timepoints)]
patient_labels[, timepoint := factor(timepoint, levels = timepoints)]
patient_labels <- patient_labels[patient != 'pat_NA']
patient_labels[, clinical_response := factor(clinical_response, 
                                             levels = c('NR', 'R'))]
blood_adaptive <- manual_clinical_corrections(blood_adaptive)

if (T) {
  ## Combine timepoint and clinical response into one variable and create
  ## corresponding palette
  patient_labels[, 'comb_time_resp' :=
                 sprintf('%s-%s', timepoint, clinical_response)]
  levs <- patient_labels[, expand.grid('timepoint' = levels(timepoint),
                                       'clinical_response' = levels(clinical_response))]
  levs$timepoint <- factor(levs$timepoint, levels = timepoints)
  levs$comb <- apply(levs, 1, paste, collapse = '-')
  levs$cf <- (1.3^(as.integer(levs$timepoint) - 2))
  levs$color <- maartenutils::lighten(c(rep(resp_colors[2], 3),
                                       rep(resp_colors[1], 3)),
                                     (1.4^(as.integer(levs$timepoint) - 2)))
  comb_time_resp_palette <- setNames(levs$color, levs$comb)
  # plot_palette(comb_time_resp_palette)
  patient_labels[, comb_time_resp := factor(comb_time_resp, levels = levs$comb)]
  cond_rm(levs)
}

if (F) {
  ## Check Leonie's remarks
  setkey(patient_labels, patient)
  setkey(blood_adaptive, patient)
  patient_labels[patient == 'pat_14']
  patient_labels[patient == 'pat_19']
  patient_labels[patient == 'pat_25']
  patient_labels[patient == 'pat_3']
  patient_labels[patient == 'pat_31']
  patient_labels[patient == 'pat_36']
  patient_labels[patient == 'pat_41']
  patient_labels[patient == 'pat_42']
  patient_labels[patient == 'pat_48']
  patient_labels[patient == 'pat_61']
  patient_labels[patient == 'pat_71']

  patient_labels[patient == 'pat_33']
  patient_labels[patient == 'pat_64']
  patient_labels[patient == 'pat_69']
  patient_labels[patient == 'pat_73']
  patient_labels[patient %in% c('pat_69', 'pat_73') & is.na(arm), .N]
}

merge_tests(idx = 12)

if (T) {
  hla_types <-
    read_excel(file.path(data_dir, 'haplotyping_TONIC_rna_Steven.xlsx'),
               na = c('', 'NA')) %>%
    as.data.table %>%
    maartenutils::normalize_colnames() 
  setnames(hla_types, gsub('\\(|\\)', '', colnames(hla_types)))
  setnames(hla_types, 'study_id', 'patient')
  hla_types[, patient := sprintf('pat_%s', patient)]
  patient_labels <- 
    controlled_merge(patient_labels, hla_types[, .(patient, hla_haplotype_rna)])
  patient_labels[, c('A1', 'A2', 'B1', 'B2', 'C1', 'C2') := 
                 tstrsplit(hla_haplotype_rna, ', ')]
  unique_hlas <- patient_labels[, {
    vec <- strsplit(.SD[, hla_haplotype_rna[1]], ', ')[[1]]
    .('unique_HLAs' = ifelse(is.na(vec), as.integer(NA), uniqueN(vec)))
  }, by = patient] %>% unique
  patient_labels <- controlled_merge(patient_labels, unique_hlas)
  patient_labels[, hla_haplotype_rna := NULL]
  merge_tests(idx = 13)
}



if (F) {
  write_tsv(patient_labels, '~/Downloads/TONIC_pat_labels.tsv')
  write_tsv(blood_adaptive, '~/Downloads/TONIC_blood_pat_labels.tsv')
}
