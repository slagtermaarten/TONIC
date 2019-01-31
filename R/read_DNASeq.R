read_CONTRA <- function(patient = 'pat_1') {
  l_patient <- patient
  rel_subs <-
    patient_labels[patient == l_patient & timepoint == 'Baseline' &
                   !is.na(cf_number)]
  if (nrow(rel_subs) == 0) return(NULL)
  cf_number <- rel_subs[, cf_number]
  full_fn <- wes_table[tolower(tumor_cf) == tolower(cf_number),
                       file.path(forge_mirror, 'contra', 'CNATable',
                                 contra_fn)]
  if (length(full_fn) == 0) {
    return(NULL)
  }
  contra <- suppressWarnings(tryCatch(fread(full_fn),
              error = function(e) { print(e); browser() }))
  contra <- normalize_colnames(contra)
}


#' Find DNA CF by patient and timepoint
#'
#'
lookup_DNA_cf <- function(patient = 'pat_69', timepoint = 'Baseline') {
  l_patient <- patient
  l_timepoint <- timepoint
  rel_subs <-
    patient_labels[tolower(patient) == tolower(l_patient) &
                   timepoint == l_timepoint &
                   !is.na(cf_number)]
  if (nrow(rel_subs) == 0) return(NULL)
  cf_number <- rel_subs[, last(unique(cf_number))]
  browser(expr = length(cf_number) > 1)
  return(cf_number)
}


lookup_WES_fn <- function(patient = 'pat_65', timepoint = 'Baseline') {
  cf_number <- lookup_DNA_cf(patient = patient, timepoint = timepoint)
  if (is.null(cf_number)) return(NULL)
  # browser(expr = wes_tabpe[tumor_cf == cf_number, .N == 2])
  ## Pat 65 inadvertently got a technical replication which should be slightly
  ## better than the original. Both files are included in this overview.
  ## Hence the last() function call.
  fn <- wes_table[tolower(tumor_cf) == tolower(cf_number), unique(vcf_fn)]
  if (is.null(fn) || length(fn) == 0 || is.na(fn))
    return(NULL)
  return(file.path(forge_mirror, 'calls', last(fn)))
}


filter_VCF <- function(patient = 'pat_1', timepoint = 'Baseline',
                       read_annotated = T,
                       tumor = T,
                       tlod_threshold = 40, nlod_threshold = 10) {
  if (tumor == T) {
    full_fn <- lookup_WES_fn(patient = patient, timepoint = timepoint)
  } else {
    ## Look up tumor CF first
    l_tumor_cf <- lookup_DNA_cf(patient = patient, timepoint = timepoint)
    if (is.null(l_tumor_cf)) return(NULL)
    l_normal_cf <- wes_table[tolower(tumor_cf) == tolower(l_tumor_cf), normal_cf]
    full_fn <- file.path(forge_mirror, 'haplotypecaller-q100',
                         wes_table[tolower(tumor_cf) == tolower(l_tumor_cf),
                                   last(germline_vcf)])
  }

  if (is.null(full_fn) || length(full_fn) == 0) {
    return(NULL)
  }

  if (read_annotated) {
    full_fn <- gsub('\\.vcf$', '.annotated.vcf', full_fn)
  }

  if (!file.exists(full_fn)) {
    warningf('%s; %s does not exist', patient, full_fn)
    return(NULL)
  }
  com <- sprintf("grep -v '##' %s ", full_fn)
  vcf <- suppressWarnings(
           tryCatch(fread(cmd = com, fill = T, sep = '\t'),
           error = function(e) { print(e); return(NULL) }))
  if (null_dat(vcf)) return(NULL)

  setnames(vcf, gsub('#', '', colnames(vcf)))
  if (tumor) {
    vcf <- vcf[grepl('SOMATIC', INFO)]
    vcf[, 'tlod' := as.numeric(gsub('.*TLOD=(\\d*\\.*\\d*);.*', '\\1', INFO))]
    vcf[, 'nlod' := as.numeric(gsub('.*NLOD=(\\d*\\.*\\d*);.*', '\\1', INFO))]
    vcf <- vcf[tlod >= tlod_threshold & nlod >= nlod_threshold]
    vcf <- vcf[FILTER == 'PASS']
  }
  vcf[, 'ANN' := gsub('ANN\\=(.*)\\=*', '\\1', INFO)]
  return(vcf)
}


read_VCF_less <- function(patient = 'pat_1', timepoint = 'Baseline') {
  full_fn <- lookup_WES_fn(patient = patient, timepoint = timepoint)
  if (length(full_fn) == 0) { return(NULL) }
  less(full_fn)
}


#' Filter germline calls for genes of interest
#'
#'
filter_germline_gene_symbol <- function(patient = 'pat_1',
                                        timepoint = 'Baseline',
                                        gene_symbols = c('BRCA1', 'BRCA2')) {
  germline_calls <- filter_VCF(patient = patient, timepoint = timepoint,
                               read_annotated = T, tumor = F)
  if (null_dat(germline_calls)) return(NULL)
  # sr(entrez_table)
  entrez_table <- readRDS(file.path(rds_dir, 'entrez_table.rds'))
  coordinates <- entrez_table[hgnc_symbol %in% gene_symbols]
  germline_calls[, 'end_position' := POS + stringr::str_length(ALT)]
  setnames(germline_calls, 'POS', 'start_position')
  setnames(germline_calls, 'CHROM', 'chromosome')
  germline_calls[, 'donor_id' := patient]
  coordinates[, 'donor_id' := patient]
  setnames(coordinates, 'chromosome_name', 'chromosome')
  germline_calls[, 'variant_classification' := 'missense_mutation']
  germline_calls[, chromosome := as.character(chromosome)]
  germline_calls[, start_position := as.integer(start_position)]
  germline_calls[, end_position := as.integer(end_position)]
  coordinates[, chromosome := as.character(chromosome)]
  coordinates[, start_position := as.integer(start_position)]
  coordinates[, end_position := as.integer(end_position)]
  germline_calls <- germline_calls[chromosome %in%
                                   coordinates[, unique(chromosome)]]
  setkey(coordinates, chromosome, start_position, end_position)
  setkey(germline_calls, chromosome, start_position, end_position)
  germline_calls <- foverlaps(coordinates, germline_calls,
                              type = 'any', which = F)
  germline_calls %<>% annotate_effects
  return(germline_calls[order(effect)])
}


read_sequenza <- function(patient = 'pat_69', timepoint = 'Baseline') {
  l_cf_number <- lookup_DNA_cf(patient = patient, timepoint = timepoint)
  wes_table[, tumor_cf == l_cf_number]
  # l_cf_number <- lookup_DNA_cf(patient = patient, timepoint = 'Baseline')
  # l_cf_number <- lookup_DNA_cf(patient = patient, timepoint = 'Post-induction')
  # l_cf_number <- lookup_DNA_cf(patient = patient, timepoint = 'On nivo')
  if (is.null(l_cf_number)) return(NULL)
  # wes_table[, tolower(normal_cf) == tolower(l_cf_number)]
  full_fn <- tryCatch(wes_table[tolower(tumor_cf) == tolower(l_cf_number),
                      unique(sequenza_fn)],
                      error = function(e) { print(e); browser() })
  if (is.null(full_fn) || length(full_fn) == 0)
    return(NULL)
  if (length(full_fn) > 1) browser()
  fh <- tryCatch(fread(full_fn), error = function(e) { print(e); browser() })
  if (null_dat(fh)) return(NULL)
  setnames(fh, gsub('\\.', '_', colnames(fh)))
  fh[, 'seg_length' := end_pos - start_pos]
  fh[, 'LOH' := A == 0 | B == 0]
  sol_fn <- view_sequenza(patient = patient,
                          timepoint = timepoint,
                          ft_rep = 'alternative_solutions.txt',
                          open_directly = F)[1]
  if (is.null(sol_fn) || is.na(sol_fn) || length(sol_fn) != 1) browser()
  fh %<>% cbind(fread(sol_fn)[1])
  # fh[, c('CNt_n', 'A_n', 'B_n') := list(CNt / ploidy, A / ploidy, B / ploidy)]
  fh[, 'CNt_n' := CNt / ploidy]

  ## Merge in arm and band information
  fh <- merge_karyotyping(fh)
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
  full_fn <- wes_table[tolower(tumor_cf) == tolower(l_cf_number), sequenza_fn]
  if (is.null(full_fn) || is.na(full_fn) || length(full_fn) == 0) return(NULL)
  new_fn <- gsub('segments.txt$', ft_rep, full_fn)
  if (open_directly) sys_file_open(new_fn)
  return(new_fn)
}


intersect_sequenza <- function(patient = 'pat_69',
                               timepoint = 'Baseline',
                               coordinates) {
  fh <- tryCatch(read_sequenza(patient = patient, timepoint = timepoint),
                 warning = function(e) { print(e); browser() })
  if (null_dat(fh)) return(NULL)
  seg_names <- c('seg_start', 'seg_end')
  setnames(fh, c('start_pos', 'end_pos'), seg_names)
  fh[, 'donor_id' := patient]
  coordinates[, donor_id := patient]
  # setnames(fh, 'hugo_symbol', 'variant_classification')
  fh[, variant_classification := 'missense_mutation']

  fh <- tryCatch(fasanalysis::merge_muts_SNP6(muts = coordinates, fh,
                                              pos_col = seg_names,
                                              merge_cols = colnames(fh)),
                 error = function(e) { print(e); browser() },
                 warning = function(w) { print(w); browser() })
  # setnames(fh, 'variant_classification', 'hugo_symbol')
  fh[, 'gain_or_loss' := ifelse(CNt > 2, 'gain', ifelse(CNt == 2, '', 'loss'))]
  fh[CNt == 0, 'gain_or_loss' := 'homozygous loss']
  fh[CNt > 0, 'variant_classification' := sprintf('%s %s', gain_or_loss,
                                           ifelse(LOH, 'LOH', ''))]
  fh[CNt == 0, 'variant_classification' := gain_or_loss]
  fh[, variant_classification := gsub('^ ', '', variant_classification)]
  fh[, variant_classification := gsub('\\s{2,}', ' ', variant_classification)]
  fh[, variant_classification := ifelse(variant_classification == '', '',
                            sprintf('chrom %s', variant_classification))]
  setnames(fh, 'donor_id', 'patient')
  setcolorder(fh, c('patient', setdiff(colnames(fh), 'patient')))
  fh <- fh[!is.na(seg_start)]
  fh <- clean_columns('', fh, c('snp6_amp', 'coverage'))
  return(fh)
}


#' Annotate variants with most deleterious effect
#'
#' @param fh, object of class \code{data.table} returned by filter_VCF
annotate_effects <- function(fh) {
  if (null_dat(fh)) return(NULL)
  ## INFO=<ID=ANN,Number=.,Type=String,Description="Functional
  ## annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name |
  ## Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank |
  ## HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length |
  ## AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ">
  # fh[1, ANN]
  for (i in 1:nrow(fh)) {
    all_effects <- fh[i, strsplit(ANN, ';')][[1]]
    ann_IDX <- grep('\\|', all_effects)
    all_effects <- tryCatch(strsplit(strsplit(all_effects[ann_IDX], ',')[[1]],
                                     '\\|'),
                            error = function(e) { print(e); return(NULL) })
    if (is.null(all_effects)) next
    all_effects <- as.data.frame(t(sapply(all_effects, function(x)
                                          setNames(x[c(4, 2, 3, 10, 11)],
                                                   c('hugo_symbol',
                                                     'variant_classification',
                                                     'severity',
                                                     'bp_change',
                                                     'aa_change')))))
    ## Limit to first effect if multiple effects on transcript
    all_effects$variant_classification %<>% { gsub('&.*$', '', .) }
    severity_levs <- c('LOW', 'MODIFIER', 'MODERATE', 'HIGH')
    browser(expr = !all(all_effects$severity %in% severity_levs))
    all_effects$severity <- factor(all_effects$severity, levels = severity_levs)
    all_effects <- all_effects %>%
      dplyr::arrange(desc(severity)) %>%
      unique %>%
      dplyr::group_by(hugo_symbol) %>%
      { suppressMessages(dplyr::top_n(., 1)) }

    fh[i, 'effect' := all_effects$severity[1]]
    fh[i, 'bp_change' := all_effects$bp_change[1]]
    fh[i, 'aa_change' := all_effects$aa_change[1]]
    fh[i, 'hugo_symbol' := all_effects$hugo_symbol[1]]
    fh[i, 'variant_classification' := all_effects$variant_classification[1]]
  }
  fh[, variant_classification := simple_cap(variant_classification)]
  return(fh)
}


#' Filter VCF list down to genes of interest
#'
#'
filter_var_list <- function(l_patient,
                            # gs_search = paste(positions[, unique(gene_symbol)],
                                     # collapse = '|')) {
                            gs_search = paste(positions[, unique(gene_id)],
                                              collapse = '|')) {
  messagef('Filtering VCF for %s', l_patient)
  fh <- filter_VCF(patient = l_patient, read_annotated = T, tumor = T)
  if (is.null(fh)) { return(NULL) }
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


merge_karyotyping <- function(fh) {
  packrat::extlib('GenomicRanges')
  packrat::extlib('IRanges')
  chrom_dat <-
    fread(file.path(data_dir, 'ideogram_9606_GCF_000001305.14_850_V1')) %>%
    setnames('#chromosome', 'chromosome')
  # print(chrom_dat[, .(min(bp_start), max(bp_stop)), by = chromosome])

  fh_ranges <-
    GRanges(seqnames = Rle(fh[, chromosome]),
            ranges = fh[, IRanges(start = start_pos, end = end_pos)])
  cr_ranges <-
    GRanges(seqnames = Rle(chrom_dat[, chromosome]),
            ranges = chrom_dat[, IRanges(start = bp_start, end = bp_stop)],
            stain = chrom_dat[, stain],
            band = chrom_dat[, band],
            arm = Rle(chrom_dat[, arm]))

  IRanges::findOverlaps(fh_ranges, cr_ranges, type = 'within', select = 'all')
  i_m <- IRanges::mergeByOverlaps(fh_ranges, cr_ranges,
                                  type = 'any', select = 'all') %>%
    as.data.table(i_m) %>%
    setnames(c('fh_ranges.seqnames', 'fh_ranges.start', 'fh_ranges.end'),
             c('chromosome', 'start_pos', 'end_pos')) %>%
    .[, setdiff(colnames(.),
                c('cr_ranges.strand', 'fh_ranges.strand', 'cr_ranges.start',
                  'cr_ranges.seqnames', 'cr_ranges.width',
                  'cr_ranges.end', 'cr_ranges.width')),
      with = F]

  ## Identify segments overlapping a centromere; don't merge these in
  bl_segs <- i_m[, .('u_arms' = uniqueN(cr_ranges.arm)),
      by = .(chromosome, start_pos, end_pos)] %>%
        .[u_arms > 1, .(chromosome, start_pos, end_pos)]
  setkey(i_m, chromosome, start_pos, end_pos)
  fh %<>% controlled_merge(i_m[!bl_segs], all.x = T, all.y = F)
  return(fh)
}


#' Combine all aberrations for all patients in a user specified set of genes
#'
#'
combine_all_aberrations <- function(gene_symbols = 'B2M') {
  messagef('Quering gene locations...')
  ensg_mapping <- fread(file.path(p_root, 'salmon_rna',
                                  'ensembl86_ensg_to_symbol.tsv'))
  ensg_mapping[, 'variant_classification' := hugo_symbol]
  ensg_mapping <- ensg_mapping[hugo_symbol %in% gene_symbols, ]
  ensg_mapping[, c('start_position', 'end_position') :=
               list(min(start_position), max(end_position)), by = gene_id]
  positions <- unique(ensg_mapping, by = 'gene_id')
  browser(expr = all(!gene_symbols %in% positions[, unique(hugo_symbol)]),
          text = 'Not all genes found')

  message('Filtering CNAs...')
  str_pres <- rbindlist(lapply(patient_labels[, unique(patient)],
                               function(l_patient) {
      intersect_sequenza(patient = l_patient, coordinates = positions)
    }), fill = T) %>%
    controlled_merge(patient_labels[, .(patient, clinical_response, arm)]) %>%
    { .[naturalorder(patient), ] }

  message('Filtering VCFs...')
  # l_gs_search <- paste(positions[, unique(gene_id)], collapse = '|')
  l_gs_search <- paste(gene_symbols, collapse = '|')
  mut_pres <- rbindlist(lapply(patient_labels[, unique(naturalsort(patient))],
                               filter_var_list, gs_search = l_gs_search),
                        fill = T)
  if (nrow(mut_pres) > 0) {
    mut_pres <- mut_pres[hugo_symbol %in% positions$hugo_symbol,
                         .(patient, arm, clinical_response, hugo_symbol,
                           variant_classification, effect, ID)]
  } else {
    mut_pres <- NULL
  }

  gene_dat <- rbind(
    str_pres[hugo_symbol %in% positions$hugo_symbol,
            .(patient, arm, clinical_response, gene_id, hugo_symbol,
              variant_classification, CNt_n, CNt, A, B)],
    mut_pres, fill = T) %>% { .[naturalorder(patient)][order(hugo_symbol)] }

  setkey(gene_dat, patient)
  ## Make sure all patients are represented
  gene_dat <- gene_dat[patient_labels[, naturalsort(unique(patient))]]
  gene_dat[, 'mutated' := any(variant_classification != ''), by = patient]
  gene_dat <- unique(gene_dat)
  messagef('%f of patients aberrated for gene%s: %s',
           unique(gene_dat, by = 'patient')[, mean(mutated, na.rm = T)],
           ifelse(length(gene_symbols) > 1, 's', ''),
           paste(gene_symbols, collapse = ', '))

  rm_white <- function(x) gsub('^\\s*|\\s*$', '', x)
  messagef('%s',
           apply(gene_dat[, .N, by = variant_classification], 1,
                 function(x) sprintf('%s: %s', rm_white(x[[1]]),
                                     rm_white(x[[2]]))) %>%
           paste(collapse = ', '))
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
                           open_directly = F) %>% unique
  if (!is.null(sols_fn) && !is.na(sols_fn) && length(sols_fn) > 0) {
    ploidy_est <- fread(sols_fn[1])[1, .('purity' = cellularity, ploidy)]
    if (is.null(ploidy_est) || length(ploidy_est) == 0) return(NULL)
    # browser(expr = !eps(ploidy_est, fh[, weighted.mean(CNt, seg_length)],
    #                     epsilon = .25 * ploidy_est))
  } else {
    ploidy_est <- list('purity' = NA, 'ploidy' = NA)
  }

  chrom_scores <-
    fh[, .('chrom_gain' = pmin(quantile(CNt, .1) - 2),
           'chrom_loss' = 2 - pmax(quantile(CNt, .9))), by = chromosome] %>%
    .[, 'chrom_gain_b' := chrom_gain > 0] %>%
    .[, 'chrom_loss_b' := chrom_loss > 0]

  chrom_score <-
    chrom_scores[, .('gained_alleles' = .SD[chrom_gain > 0, sum(chrom_gain)],
                     'lost_alleles' = .SD[chrom_loss > 0, sum(chrom_loss)])] %>%
      unlist %>% sum

  arm_scores <-
    fh[, .('arm_gain' = pmin(quantile(CNt, .1) - 2),
           'arm_loss' = 2 - pmax(quantile(CNt, .9))),
       by = .(chromosome, arm)] %>%
    .[, 'arm_gain_b' := arm_gain > 0] %>%
    .[, 'arm_loss_b' := arm_loss > 0]
  arm_score <-
    arm_scores[, .('gained_alleles' = .SD[arm_gain > 0, sum(arm_gain)],
                   'lost_alleles' = .SD[arm_loss > 0, sum(arm_loss)])] %>%
      unlist %>% sum

  fh[, 'CN_score_segment' := pmin(abs(CNt - 2), 2)]
  fh[, 'CN_score_segment_unbounded' := abs(CNt - 2)]

  ## Identify chroms with only one whole chromosome event
  blacklist_chroms <- fh[, .N, by = chromosome][N == 1, sort(chromosome)]
  if (length(blacklist_chroms) > 0) {
    messagef('%s; following chroms have only one segment: %s', patient,
             paste(blacklist_chroms, collapse = ', '))
  }

  ## Try to compute LST, difficult since we only have exonic coverage
  ## TODO account for centromeres and chromosome length
  # fh[, 'padding_5' :=
  #    # c(.SD[, start_pos][1:(.N-1)], max(end_pos)) - c(0, .SD[, end_pos][2:.N]),
  #    .SD[, start_pos] - c(0, .SD[, end_pos][1:(.N-1)]),
  #    by = chromosome]
  # fh[, 'padding_3' := c(.SD[, padding_5][2:.N], 0), by = chromosome]
  # fh[, 'A_not_equal_previous_A' := .SD[, A] != c(0, .SD[, A][1:(.N-1)]),
  #    by = chromosome]
  # fh[, 'B_not_equal_previous_B' := .SD[, B] != c(0, .SD[, B][1:(.N-1)]),
  #    by = chromosome]
  # fh[chromosome == '20']
  # fh[chromosome == '20' & padding_5 <= 3e6]
  # fh[chromosome == '1']

  vcf_fh <- filter_VCF(patient = patient)
  if (null_dat(vcf_fh)) return(NULL)
  vcf_fh %<>% annotate_effects

  return(c(list('patient' = patient,
         'mut_load' = vcf_fh[variant_classification %in% mut_load_var_types, .N],
         'rs_frac' = vcf_fh[, length(grep('rs', ID))] / nrow(vcf_fh),
         'gen_SCNA_score' = fh[, sum(CN_score_segment)],
         'chrom_SCNA_score' = chrom_score,
         'arm_SCNA_score' = arm_score,
          ## Abkevich et al., (2012) doi:10.1038/bjc.2012.451
         'HRD_LOH' = fh[chromosome %nin% blacklist_chroms &
                        seg_length >= 15000000 & LOH == T, .N],
          ## Popova et al., (2012) doi:10.1158/0008-5472.CAN-12-1470
         'ntAI' = fh[, .SD[c(1, .N)], by = chromosome] %>%
            { .[CNt_n > 0 & A != B][, .N] },
         'gen_SCNA_score_unbounded' = fh[, sum(CN_score_segment_unbounded)],
         'weighted_gen_SCNA_score' = fh[, sum(CN_score_segment * seg_length /
                                              max(seg_length))],
         'perc_LOH' = fh[, mean(as.integer(LOH))],
         'weighted_perc_LOH' = fh[, weighted.mean(as.integer(LOH), seg_length)]),
    ploidy_est))
}


draw_lolliplot <- function(variant_list = brca_variants,
                           hugo_symbol = 'BRCA1') {
  l_hugo_symbol = hugo_symbol
  varlist <- brca_variants[hugo_symbol == l_hugo_symbol & aa_change != '',
                sprintf('%s@%d',
                        # gsub('p.[A-z]+(\\d+)[A-z]+$', '\\1', aa_change),
                        gsub('p.([A-z]+)(\\d+)([A-z]+)$', '\\1\\2\\3', aa_change),
                        .N),
                by = aa_change] %>% { paste(.[, V1], collapse = ' ') }
  lolly_loc <- '/Users/maartenslagter/libs/lollipops_1.4.0_mac64/lollipops'
  system(sprintf('%s -w=700 -legend -labels %s %s', lolly_loc,
                 hugo_symbol, varlist))
  invisible()
}
