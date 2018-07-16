read_CONTRA <- function(patient = 'pat_1') {
  rel_subs <-
    patient_labels[patient == parent.frame(3)$patient &
                   timepoint == 'Baseline' &
                   !is.na(cf_number)]
  if (nrow(rel_subs) == 0) return(NULL)
  cf_number <- rel_subs[, cf_number]
  full_fn <- wes_table[tumor_cf == cf_number,
                       file.path(forge_mirror, 'contra', 'CNATable',
                                 contra_fn)]
  if (length(full_fn) == 0) {
    return(NULL)
  }
  contra <- suppressWarnings(tryCatch(fread(full_fn),
              error = function(e) { print(e); browser() }))
  contra <- normalize_colnames(contra)
}


lookup_DNA_cf <- function(patient = 'pat_69', timepoint = 'Baseline') {
  ## Find DNA CF
  l_patient = patient
  # patient_labels[timepoint == 'Baseline' & !is.na(cf_number), .N]
  # patient_labels[timepoint == 'Baseline' & !is.na(cf_number), patient]
  # patient_labels[patient == parent.frame(3)$patient]
  rel_subs <-
    patient_labels[tolower(patient) == tolower(parent.frame(3)$patient) &
                   timepoint == parent.frame(3)$timepoint &
                   !is.na(cf_number)]
  if (nrow(rel_subs) == 0) return(NULL)
  cf_number <- rel_subs[, unique(cf_number)]
  browser(expr = length(cf_number) > 1)
  return(cf_number)
}


lookup_WES_fn <- function(patient = 'pat_65', timepoint = 'Baseline') {
  cf_number <- lookup_DNA_cf(patient = patient, timepoint = timepoint)
  if (is.null(cf_number)) return(NULL)
  browser(expr = wes_table[tumor_cf == cf_number, .N == 2])
  full_fn <-
    wes_table[tumor_cf == cf_number, file.path(forge_mirror, 'calls', vcf_fn)]
  return(full_fn)
}


filter_VCF <- function(patient = 'pat_1', timepoint = 'Baseline',
                       read_annotated = T,
                       tlod_threshold = 40, nlod_threshold = 10) {
  full_fn <- lookup_WES_fn(patient = patient, timepoint = timepoint)
  if (length(full_fn) == 0) { return(NULL) }
  # vcf <- filter_VCF(full_fn, 
  #                   tlod_threshold = tlod_threshold, 
  #                   nlod_threshold = nlod_threshold)
  # if (null_dat(vcf)) return(NULL)
  if (read_annotated) {
    full_fn <- gsub('\\.vcf$', '.annotated.vcf', full_fn)
  }
  if (!file.exists(full_fn)) {
    warningf('%s; %s does not exist', patient, full_fn)
    return(NULL)
  }
  com <- sprintf("grep SOMATIC %s | grep -v '#'", full_fn)
  vcf <- suppressWarnings(
           tryCatch(fread(com, fill = T, sep = '\t'),
           error = function(e) { print(e); return(NULL) }))
  if (null_dat(vcf)) return(NULL)
  # vcf <- vcf[!grepl('^#|GERMLINE', 1)]
  # vcf <- vcf[!grepl('.*##.*', 1)]
  # vcf <- vcf[!grepl('.*GERMLINE.*', 1)]
  vcf_colnames <- system(sprintf('cat %s | grep CHROM', full_fn), 
                         intern = T)
  vcf_colnames <- gsub('#', '', vcf_colnames)
  vcf_colnames <- strsplit(x = vcf_colnames, split = '\\t')
  vcf_colnames <- tryCatch(vcf_colnames[[1]], error = function(e) { print(e); browser() }) 

  browser(text = 'Unexpected column length of VCF (1)', 
          expr = length(colnames(vcf)) != length(vcf_colnames))
  setnames(vcf, vcf_colnames)
  browser(text = 'Unexpected column length of VCF (2)', 
          expr = length(colnames(vcf)) != 11)
  vcf[, 'tlod' := as.numeric(gsub('.*TLOD=(\\d*\\.*\\d*);.*', '\\1', INFO))]
  vcf[, 'nlod' := as.numeric(gsub('.*NLOD=(\\d*\\.*\\d*);.*', '\\1', INFO))]
  vcf <- vcf[tlod >= tlod_threshold & nlod >= nlod_threshold]
  vcf <- vcf[FILTER == 'PASS']
  # com <- sprintf("grep '##INFO' %s", full_fn)
  # td <- fread(com, fill = T)
  # info_field_tags <- gsub('.*\\=(\\w*)', '\\1', unlist(td[, 1, with = F])) %>%
  #   setNames(NULL)
  # browser()
  # less(full_fn)
  vcf[, 'ANN' := gsub('ANN\\=(.*)\\=*', '\\1', INFO)]
  # vcf[, info]
  # sapply(strsplit(vcf[1, INFO], ';')[[1]], function(x) { 
  #   setNames(gsub('(.*)\\=(.*)', '\\1', x)
  # })
  # vcf[, (info_field_tags) := tstrsplit(INFO, '\\;')]

  # vcf[, 'ensg' := gsub('(ENSG=\\d*)', '\\1', INFO)]
  return(vcf)
}


read_VCF_less <- function(patient = 'pat_1', timepoint = 'Baseline') {
  full_fn <- lookup_WES_fn(patient = patient, timepoint = timepoint)
  if (length(full_fn) == 0) { return(NULL) }
  less(full_fn)
}


read_sequenza <- function(patient = 'pat_69', timepoint = 'Baseline') {
  l_cf_number <- lookup_DNA_cf(patient = patient, timepoint = timepoint)
  if (is.null(l_cf_number)) return(NULL)
  full_fn <- tryCatch(wes_table[tumor_cf == l_cf_number, sequenza_fn], 
                      error = function(e) { print(e); browser() }) 
  if (is.null(full_fn) || length(full_fn) == 0) 
    return(NULL)
  fh <- tryCatch(fread(full_fn), error = function(e) { print(e); browser() }) 
  if (null_dat(fh)) return(NULL)
  setnames(fh, gsub('\\.', '_', colnames(fh)))
  fh[, 'seg_length' := end_pos - start_pos]
  fh[, 'LOH' := A == 0 | B == 0]
  return(fh)
}


#' View a Sequenza file by looking up the corresponding table file and altering
#' it to represent the desired file
#'
#'
view_sequenza <- function(patient = 'pat_69', timepoint = 'Baseline',
                          ft_rep = 'genome_view.pdf',
                          open_directly = T) {
  l_cf_number <- lookup_DNA_cf(patient = patient, timepoint = timepoint)
  full_fn <- wes_table[tumor_cf == l_cf_number, sequenza_fn]
  if (is.null(full_fn)) return(NULL)

  new_fn <- gsub('segments.txt$', ft_rep, full_fn)
  if (open_directly) {
    sys_file_open(new_fn)
  }
  return(new_fn)
}


intersect_sequenza <- function(patient = 'pat_69', 
                               timepoint = 'Baseline',
                               coordinates) {
  fh <- read_sequenza(patient = patient, timepoint = timepoint)
  if (null_dat(fh)) return(NULL)
  seg_names <- c('seg_start', 'seg_end')
  setnames(fh, c('start_pos', 'end_pos'), seg_names)
  fh[, 'donor_id' := patient]
  coordinates[, donor_id := patient]
  fh <- tryCatch(fasanalysis::merge_muts_SNP6(muts = coordinates, fh,
                                              pos_col = seg_names,
                                              merge_cols = colnames(fh)), 
                 error = function(e) { print(e); browser() }) 
  setnames(fh, 'variant_classification', 'hugo_symbol')
  fh[, 'gain_or_loss' := ifelse(CNt > 2, 'gain', ifelse(CNt == 2, '', 'loss'))]
  fh[, 'variant_classification' := sprintf('%s %s', gain_or_loss, 
                                           ifelse(LOH, 'LOH', ''))]
  fh[, variant_classification := gsub('^ ', '', variant_classification)]
  fh[, variant_classification := gsub('\\s{2,}', ' ', variant_classification)]
  fh[, variant_classification := ifelse(variant_classification == '', '',
                            sprintf('chrom %s', variant_classification))]
  setnames(fh, 'donor_id', 'patient')
  setcolorder(fh, c('patient', setdiff(colnames(fh), 'patient')))
  return(fh)
}


#' Annotate variants with most deleterious effect
#'
#' @param fh, object of class \code{data.table} returned by filter_VCF
annotate_effects <- function(fh) {
  ##INFO=<ID=ANN,Number=.,Type=String,Description="Functional
  ## annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name |
  ## Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank |
  ## HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length |
  ## AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ">
  # fh[1, ANN]
  for (i in 1:nrow(fh)) {
    all_effects <- fh[i, strsplit(ANN, ';')][[1]]
    ann_IDX <- grep('\\|', all_effects)
    all_effects <- strsplit(strsplit(all_effects[ann_IDX], ',')[[1]], '\\|')
    all_effects <- as.data.frame(t(sapply(all_effects, function(x) 
                                          setNames(x[c(4, 2, 3)],
                                                   c('hugo_symbol',
                                                     'variant_classification',
                                                     'severity')))))
    severity_levs <- c('LOW', 'MODIFIER', 'MODERATE', 'HIGH')
    browser(expr = !all(all_effects$severity %in% severity_levs))
    all_effects$severity <- factor(all_effects$severity, levels = severity_levs)
    all_effects <- all_effects %>% 
      dplyr::arrange(desc(severity)) %>% 
      unique %>%
      dplyr::group_by(hugo_symbol) %>%
      dplyr::top_n(1)
    fh[i, 'effect' := all_effects$severity[1]]
    fh[i, 'hugo_symbol' := all_effects$hugo_symbol[1]]
    fh[i, 'variant_classification' := all_effects$variant_classification[1]]
  }
  return(fh)
}


#' Filter VCF list down to genes of interest
#'
#'
filter_var_list <- function(l_patient, 
                            # gs_search = paste(positions[, unique(gene_symbol)], 
                                     # collapse = '|')) {
                            gs_search = paste(positions[, unique(ensg)], 
                                              collapse = '|')) {
  messagef('Filtering VCF for %s', l_patient)
  fh <- filter_VCF(patient = l_patient, read_annotated = T)
  if (is.null(fh)) { return(NULL) }
  ## Pre filter variants using regex
  fh <- fh[grepl(gs_search, INFO)]
  if (null_dat(fh)) { return(NULL) }
  fh <- annotate_effects(fh)
  fh[, 'patient' := l_patient]
  fh[, INFO := NULL]
  fh[, ANN := NULL]
  fh <- controlled_merge(fh, 
                         patient_labels[, .(patient, arm, clinical_response)])
  return(fh)
}


#' Combine all aberrations for all patients in a user specified set of genes
#'
#'
combine_all_aberrations <- function(positions) {
  message('Filtering CNAs...')
  str_pres <- rbindlist(lapply(patient_labels[, unique(patient)], 
                                function(l_patient) {
      intersect_sequenza(patient = l_patient, 
                         coordinates = positions)
    })) %>%
    controlled_merge(patient_labels[, .(patient, clinical_response, arm)]) %>%
    { .[naturalorder(patient), ] }

  message('Filtering VCFs...')
  l_gs_search <- paste(positions[, unique(ensg)], collapse = '|')
  mut_pres <- rbindlist(lapply(patient_labels[, unique(patient)],
                               filter_var_list, gs_search = l_gs_search))
  ## Pat 6, 9, 13, 15, 19 no VCFs

  gene_dat <- rbind(
    str_pres[hugo_symbol %in% positions$gene_symbol, 
            .(patient, arm, clinical_response, hugo_symbol, 
              variant_classification)],
    mut_pres[hugo_symbol %in% positions$gene_symbol, 
             .(patient, arm, clinical_response, hugo_symbol, 
               variant_classification)]
  ) %>% { .[naturalorder(patient)] }

  setkey(gene_dat, patient)
  ## Make sure all patients are represented
  gene_dat <- gene_dat[patient_labels[, naturalsort(unique(patient))]]
  gene_dat[, 'mutated' := any(variant_classification != ''), by = patient]
  return(gene_dat)
}


compute_dnaseq_stats <- function(patient = 'pat_33') {
  fh <- read_sequenza(patient = patient)
  # view_sequenza(patient = 'pat_69', timepoint = 'Baseline',
  #                          ft_rep = 'genome_view.pdf',
                           # open_directly = T)
  if (null_dat(fh)) return(NULL)

  ## Check ploidy estimate against the one I compute
  sols_fn <- view_sequenza(patient = patient, timepoint = 'Baseline',
                           ft_rep = 'alternative_solutions.txt',
                           open_directly = F)
  ploidy_est <- fread(sols_fn)[1, ploidy]
  if (is.null(ploidy_est) || length(ploidy_est) == 0) return(NULL)
  browser(expr = !eps(ploidy_est, fh[, weighted.mean(CNt, seg_length)], 
                      epsilon = .25 * ploidy_est))

  fh[, 'CN_score_segment' := abs(CNt - 2)]

  vcf_fh <- filter_VCF(patient = patient)
  
  list('patient' = patient, 
       'mut_load' = vcf_fh[, .N],
       'gen_SCNA_score' = fh[, mean(CN_score_segment)],
       'weighted_gen_SCNA_score' = fh[, weighted.mean(CN_score_segment,
                                                      seg_length)],
       'perc_LOH' = fh[, mean(as.integer(LOH))],
       'weighted_perc_LOH' = fh[, weighted.mean(as.integer(LOH), seg_length)])
}
