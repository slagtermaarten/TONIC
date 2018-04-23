pacman::p_load(readxl)
pacman::p_load(naturalsort)
# library(edgeR)
# d <- DGEList(counts=yourcounts)
# d <- calcNormFactors(d)
# v <- voom(d, design)

timepoints <- c('baseline' = 'Baseline', 'post.induction' = 'Post-induction',
                'on.nivo' = 'On nivo')
timepoints_inv <- setNames(names(timepoints), timepoints)
blood_timepoints <- auto_name(c(-2, 0, 6, 10, 12))
treatment_arms <- c('No induction', 'Radiotherapy', 'Cyclophosphamide',
                    'Cisplatin', 'Doxorubicin')

resp_colors <- maartenutils::gen_color_vector('Zissou1', 2) %>%
  darken(factor = c(1.0, 1.2)) %>%
  setNames(c('NR', 'R'))
# plot_palette(resp_colors)

## A primary source of the most relevant clinical information
patient_labels <- read.csv(file.path(p_root, 'data-raw/patient_labels.csv'),
                           dec = ',', sep = ';') %>% as.data.table %>%
  normalize_colnames()
setnames(patient_labels, 'x', 'filename')
setnames(patient_labels, gsub('\\.', '_', colnames(patient_labels)))
maartenutils::set_dt_types(patient_labels,
                           c('mean_log2_hk' = 'numeric',
                             'tis_score' = 'numeric'))
patient_labels[, patient := paste0('pat_', patient)]
stopifnot(patient_labels[is.na(arm), .N <= 3])

if (F) {
  ## More recent overview of clinical response data, merge into existing overview
  response_data <- read.csv(file.path(p_root, 'data-raw/response_data.csv'),
                            dec = ',', sep = ';') %>% as.data.table %>%
    normalize_colnames %>%
    ## Select last columns containing patient_id and RECIST labels
    { rev(.)[, c(2, 3)] } %>%
    { .[!is.na(patient.id)] } %>%
    setnames(., c('response', 'patient')) %>%
    ## Take out patients with dubious responses
    { .[!patient %in% c(61, 63, 64)] }
  response_data[, patient := as.character(paste0('pat_', patient))]
  patient_labels[, response := NULL]
  patient_labels[, patient := as.character(patient)]
  patient_labels <- merge(patient_labels, response_data, by = 'patient')
  patient_labels[patient == 'pat_63']
}
stopifnot(patient_labels[is.na(arm), .N <= 3])


if (T) {
	## Merge Adaptive TCR abundance measures
	## First clean up data...
	adaptive_sample_annotation <-
		read_excel(file.path(data_dir, 'SampleManifest_adaptive.xlsx')) %>%
		as.data.table %>%
		maartenutils::normalize_colnames()
	setnames(adaptive_sample_annotation, 
					 gsub('__', '_', colnames(adaptive_sample_annotation)))

  homogenize_arms <- function(vec) {
    vec <- tolower(vec)
    vec[vec == 'rt'] = 'radiotherapy'
    vec[vec == 'cis'] = 'cisplatin'
    vec[vec == 'doxorubicin'] = 'doxo'
    vec[vec == 'doxorubicine'] = 'dox'
    reps <- setNames(treatment_arms, c('no induction', 'radiotherapy',
                                       'cyclofosfamide', 'cisplatin', 'doxo'))
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
	patient_labels[, timepoint := as.character(timepoint)]
	# patient_labels[, class(patient)]
	# patient_labels[, class(timepoint)]
	# tumor_adaptive[, class(timepoint)]
	# tumor_adaptive[, class(patient)]
	## Then merge adaptive sample IDs...
	patient_labels <- merge(patient_labels, tumor_adaptive,
													by = c('patient', 'timepoint', 'arm'), 
													all.x = T, all.y = T)
  stopifnot(patient_labels[is.na(arm), .N <= 3])
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
	}

  patient_labels[, clinical_response := zoo::na.locf(clinical_response,
                                                     fromLast = F), 
                 by = patient]
  # patient_labels[, clinical_response := zoo::na.locf(clinical_response,
  #                                                    fromLast = T), 
  #                by = patient]
  # patient_labels[, clinical_response := zoo::na.locf(clinical_response,
  #                                                    fromLast = F), 
  #                by = patient]
  # patient_labels[is.na(clinical_response)]

	if (F) {
		## Fill missing values
		patient_labels[is.na(clinical_response), unique(patient)]
		patient_labels[is.na(arm)]
		patient_labels[patient == 'pat_33']
		pacman::p_load('zoo')
		patient_labels[, arm := na.locf(arm), by = patient]
		patient_labels[, clinical_response := na.locf(clinical_response), by = patient]
	}

  ## Another patient_labels-like object but for blood timepoints
	blood_adaptive <- adaptive_sample_annotation[, 10:16]
	blood_adaptive <- blood_adaptive[!is.na(sample_name), 
								 .('adaptive_sample_name' = sample_name, 
                   'patient' = paste0('pat_', study_id_1), 
									 'blood_timepoint' = tijdspunt_in_weken)]
	# blood_adaptive[, arm := factor(arm, levels = treatment_arms)]
	blood_adaptive[, blood_timepoint := factor(blood_timepoint, 
                            levels = sort(unique(blood_timepoint)))]
  blood_adaptive <- merge(blood_adaptive,
    unique(patient_labels, by = 'patient')[, .(patient, arm, response,
                                               clinical_response, tis_score)],
    by = 'patient', all.x = T, all.y = F)
	# blood_adaptive
	# dcast(blood_adaptive, formula = patient + arm + blood_adaptive ~ .,
	#       value.var = 'sample_name')
	# dcast(blood_adaptive, formula = patient + arm + blood_adaptive ~ sample_name)
	# dcast(blood_adaptive, formula = . ~ patient + arm + blood_adaptive,
	#       drop = T, fill = NA)
	# patient_labels <- merge(patient_labels, blood_adaptive,
	#                         by = c('patient'))


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
	## This line messes up patient_labels object
	patient_labels <- merge(patient_labels, tmp, by = 'adaptive_sample_name', 
													all.x = T, all.y = F)
  blood_adaptive <- merge(blood_adaptive, tmp, by = 'adaptive_sample_name', 
													all.x = T, all.y = F)
	rm(tmp)

  tmp <- read_tsv(file = file.path(data_dir, 'AdaptiveSampleOverview.tsv')) %>%
		as.data.table %>%
		maartenutils::normalize_colnames()
  tmp[, total_t_cells := NULL]
  # tmp[order(-fraction_productive_of_cells_mass_estimate)]
  setnames(tmp, 'fraction_productive_of_cells_mass_estimate', 
           'adaptive_t_cells')
	setnames(tmp, 'sample_name', 'adaptive_sample_name')
	setnames(tmp, gsub('\'', '', colnames(tmp)))
	patient_labels <- merge(patient_labels, tmp, by = 'adaptive_sample_name', 
													all.x = T, all.y = F)
  blood_adaptive <- merge(blood_adaptive, tmp, by = 'adaptive_sample_name', 
													all.x = T, all.y = F)
	setkey(patient_labels, patient, timepoint)
}

## Final formatting
## One laesion responds, the other one doesn't
## Leave out patient
manual_clinical_corrections <- function(dtf) {
  dtf <- dtf[patient != 'pat_64']
  dtf <- dtf[patient == 'pat_33', clinical_response := 'R']
  dtf <- dtf[patient == 'pat_33', response := 'PR']
  dtf[, response := factor(response, levels = c('PD', 'SD', 'PR', 'CR'))]
  dtf[arm == 'Cyclofosfamide',  arm := 'Cyclophosphamide']
  dtf[, arm := factor(arm, levels = treatment_arms)]
  return(dtf)
}
patient_labels <- manual_clinical_corrections(patient_labels)
patient_labels[, 'timepoint_number' := match(timepoint, timepoints)]
patient_labels[, timepoint := factor(timepoint, levels = timepoints)]
patient_labels <- patient_labels[patient != 'pat_NA']
stopifnot(patient_labels[is.na(arm), .N <= 3])

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
  levs$color <- maartenutils::lighten(c(rep(resp_colors[1], 3), 
                                       rep(resp_colors[2], 3)),
                                     (1.4^(as.integer(levs$timepoint) - 2)))
  comb_time_resp_palette <- setNames(levs$color, levs$comb)
  # plot_palette(comb_time_resp_palette)
  patient_labels[, comb_time_resp := factor(comb_time_resp, levels = levs$comb)]
  cond_rm(levs)
}

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
                             c('cd8_mm2' = 'numeric',
                               's_til' = 'numeric'))
  clinical_annotation[induction_therapy == 'Cyclofosfamide',  
                      induction_therapy := 'Cyclophosphamide']
  # clinical_annotation[, table(induction_therapy)]
  clinical_annotation[, induction_therapy := factor(induction_therapy,
                                                    levels = treatment_arms)]
  clinical_annotation[, 'response_bin' := as.integer(response %in% c('CR', 'PR'))]
  clinical_annotation <- clinical_annotation[patient != 'pat_NA']
  clinical_annotation[, study_id := NULL]
  clinical_annotation[, response_bin := NULL]
  clinical_annotation[, induction_therapy := NULL]
  clinical_annotation[, response := NULL]
  stopifnot(length(intersect(colnames(clinical_annotation),
                             colnames(patient_labels))) <= 1)
  patient_labels <- controlled_merge(patient_labels, clinical_annotation, by_cols = 'patient')
}


exp_levels <- read.csv(file.path(p_root,
                                 'data-raw/normalized_gene_expression.csv'),
                       # verbose = T,
                       dec = ',',
                       sep = ';', skip = 1L) %>% as.data.table

setnames(exp_levels, gsub('.RCC$', '', colnames(exp_levels)))
setnames(exp_levels, gsub('^X', '', colnames(exp_levels)))
setnames(exp_levels, 'File.Name', 'gene_symbol')
setdiff(colnames(exp_levels), 'gene_symbol') %>%
  { setNames(rep('numeric', length(.)), .) } %>%
  set_dt_types(exp_levels, .)
maartenutils::set_dt_types(exp_levels,
                           c('mean_log2_hk' = 'numeric',
                             'tis_score' = 'numeric'))

# length(intersect(patient_labels[, filename], colnames(exp_levels)))
# length(colnames(exp_levels))
stopifnot(all(patient_labels[!is.na(filename), filename] %in% colnames(exp_levels)))


danaher_scores <- read.csv(file.path(p_root,
                                     'data-raw/danaher_geneset_scores.csv'),
                    # verbose = T,
                    dec = ',',
                    sep = ';', skip = 0L) %>% as.data.table %>%
  normalize_colnames()

## Ensure everything's numeric
setdiff(colnames(danaher_scores), c('patient', 'arm', 'response')) %>%
  { setNames(rep('numeric', length(.)), .) } %>%
  set_dt_types(danaher_scores, .)
danaher_scores.m <- melt(danaher_scores,
                         id.vars = c('patient', 'arm', 'response'))
danaher_scores.m[, 'timepoint' := gsub('.*_(.*)', '\\1', variable)]
danaher_scores.m[, variable := gsub('(.*)_.*', '\\1', variable)]
danaher_scores.m[, patient := paste0('pat_', patient)]
danaher_scores.m[, timepoint := factor(timepoints[as.character(timepoint)],
                                       levels = timepoints)]
danaher_scores.m[arm == 'Cyclofosfamide',  arm := 'Cyclophosphamide']
# danaher_scores.m[, table(arm)]
danaher_scores.m[, arm := factor(arm, levels = treatment_arms)]
set_dt_types(danaher_scores.m, c('patient' = 'factor', 'timepoint' = 'factor',
                                 'variable' = 'factor'))
# danaher_scores.m <- merge(danaher_scores.m,
#                           patient_labels[, .(patient, response, arm)],
#                           by = 'patient')

## Gene signatures
danaher_signature_genes <- read.csv(file.path(p_root,
                                              'data-raw/danaher_signature_genes.csv'),
                           dec = ',', sep = ';') %>% as.data.table %>%
  normalize_colnames()
setnames(danaher_signature_genes,
         gsub('\\.', '_', colnames(danaher_signature_genes)))
danaher_signature_genes[, 'signature_name' := gsub(' ', '_',
                                                   tolower(signature_name))]
danaher_signature_genes[, 'signature_name' := gsub('\\.|-', '_',
                                                   tolower(signature_name))]
# danaher_signature_genes[, table(signature_name)]
# danaher_signature_genes[, table(gene_name)]


## Genes to test in comparative analyses
test_genes <- exp_levels[, unique(gene_symbol)] %>%
  grep(pattern = 'POS_\\w', x = ., value = T, invert = T) %>%
  grep(pattern = 'NEG_\\w', x = ., value = T, invert = T) %>%
  grep(pattern = 'ERCC Controls', x = ., value = T, invert = T) %>%
  grep(pattern = 'Housekeeping', x = ., value = T, invert = T)

test_gene_sets <- danaher_scores.m[, levels(variable)]
