#' Pick first non-NA entry from a vector
#'
#'
vec_first_non_NA <- function(vec) {
   # stopifnot(is.vector(vec))
   res <- vec[!is.na(vec)][1]
   if (length(res) == 0) return(NA) 
   else return(res)
}


#' Read adaptive seq data
#'
#'
read_adaptive_seqs <- function(force_reload = F, patient = 'pat_1', 
                               collapse_on_aa = T,
                               p_var = 'normalized_frequency') {
  ## Take most recent file if multiple are present
  fn <- list.files(file.path(p_root, 'data-raw'), pattern = 'AdvancedQuery',
                   full.names = T) %>% last

  messagef('Reading in TCR repertoires for %s', patient)
  ## Patient specific reading of files
  l_patient <- patient
  adaptive_samples <-
    sort(c(patient_labels[patient == l_patient,
                          unique(adaptive_sample_name)],
           blood_adaptive[patient == l_patient,
                          unique(adaptive_sample_name)]))
  grep_string <-
    sprintf('$1 == "%s"', adaptive_samples) %>%
    { paste(., collapse = ' || ') }
  read_com <- sprintf("awk '(%s || NR == 1) {print $0}' %s", grep_string, fn)
  arr <- fread(cmd = read_com)
  TCR_counts <- arr[, .('unique_TCRs' = .N), sample_name]
  print(TCR_counts)
  ## Ensure TCRs for all samples were found
  browser(expr = !all(TCR_counts[, sample_name] %in% adaptive_samples))

  arr[, productive_frequency := as.numeric(productive_frequency)]
  arr <- clean_columns(instance = '', fh = arr, col_names = "reads")
  arr[amino_acid == 'na', amino_acid := NA]
  arr[productive_frequency == 'na', productive_frequency := NA]
  arr <- arr[!is.na(amino_acid) & !is.na(productive_frequency)]
  setnames(arr, 'sample_name', 'adaptive_sample_name')

  arr1 <-
    controlled_merge(arr[grepl('T_', adaptive_sample_name)],
                     patient_labels[, .(patient, adaptive_sample_name,
                                        timepoint, arm, clinical_response,
                                        response)],
                     by_cols = 'adaptive_sample_name',
                     all.x = T, all.y = F)
  arr2 <-
    controlled_merge(arr[grepl('B_', adaptive_sample_name)],
                     blood_adaptive[, .(patient, adaptive_sample_name, arm,
                                        blood_timepoint,
                                        clinical_response, response)],
                     by_cols = 'adaptive_sample_name',
                     all.x = T, all.y = F)
  cond_setnames(arr2, 'blood_timepoint', 'timepoint')
  arr <- rbind(arr1, arr2, fill = T)

  if (collapse_on_aa) {
    arr[, productive_frequency := sum(productive_frequency, na.rm = T),
        by = .(amino_acid, timepoint)]
    arr <- unique(arr, by = c('timepoint', 'amino_acid'))
  }

  arr[, 'freq_rank' := (frank(productive_frequency, ties.method = 'max') - 1)
                        / (.N - 1), by = timepoint]
  arr[, timepoint := factor(as.character(timepoint),
                            levels = c(timepoints, blood_timepoints))]

  ## Merge in 'Fraction T-cells of nucleated cells'
  arr <-
    controlled_merge(arr,
                     patient_labels[, .(adaptive_sample_name,
                                        adaptive_t_cells)],
                     by_cols = 'adaptive_sample_name')
  arr <-
    controlled_merge(arr,
                     blood_adaptive[, .(adaptive_sample_name,
                                        adaptive_t_cells)],
                     by_cols = 'adaptive_sample_name')

  arr[, 'normalized_frequency' := productive_frequency * adaptive_t_cells]
  arr[, value := get(p_var)]
  arr <- arr[!is.na(productive_frequency)]
  arr[, timepoint := factor(timepoint,
                            levels = c(timepoints, blood_timepoints))]
  arr[, 'read_count' := 
      round(productive_frequency / min(productive_frequency)), 
      by = adaptive_sample_name]
  return(arr)
}


#' Further modify output of read_adaptive_seqs() for PBL analyses
#'
#'
read_annotated_bloodadaptive_seqs <- function(patient = 'pat_1',
                                              force_reload = F,
                                              collapse_on_aa = T,
                                              p_var = 'normalized_frequency',
                                              compute_pbl_it_timepoints = T) {
  arr <- read_adaptive_seqs(patient = patient, p_var = p_var, 
                            collapse_on_aa = collapse_on_aa)
  if (null_dat(arr)) return(NULL)
  av_tps <- arr[, unique(timepoint)]
  ## We need at least one of peripheral and intratumoral samples
  if (!any(blood_timepoints %in% av_tps) && any(timepoints %in% av_tps)) {
    messagef('Only one compartment available for %s', patient)
    return(NULL)
  }

  ## To 'annotate' CDR3s found in the blood with their presence in the tumor,
  ## we merge the two tables and allow multiple joins per blood observation.
  t1 <- arr[timepoint %in% blood_timepoints]
  t2 <- arr[timepoint %in% timepoints,
            .(amino_acid, timepoint, normalized_frequency, freq_rank)]
  setnames(t2,
           c('timepoint', 'normalized_frequency', 'freq_rank'),
           ## The first denotes the timepoint at which the blood TCR was
           ## detected intatumorally, the second the magnitude with which this
           ## happened
           c('it_timepoint', 'it_normalized_frequency', 'it_freq_rank'))
  if (F && any(blood_timepoints %in% av_tps) && any(timepoints %in% av_tps)) {
    browser()
  }
  ## We could get the following result for some CDR3s:
  ## b1 t1
  ## b1 t2
  ## b1 t3
  ## i.e., a blood observation that's duplicated because it's found at multiple
  ## tumor timepoints. We need to ensure this CDR3 is not represented multiple
  ## times in the final output.
  arr <- merge(t1, t2, all.x = T, all.y = F, allow.cartesian = T,
               by = 'amino_acid')

  suppressWarnings(arr[, shared_it_timepoints := NULL])
  arr[!is.na(it_timepoint), 
      'shared_it_timepoints' := as.integer(sum(value > 0)),
      by = .(timepoint, amino_acid)]
  # arr[, table(timepoint)]
  # arr[, table(it_timepoint)]
  # arr[, .SD[any(!is.na(it_timepoint))], by = amino_acid]
  # arr[!is.na(it_timepoint), unique(amino_acid)] %>%
  #   { arr[amino_acid %in% .][order(amino_acid)] }
  arr[, shared_it_timepoints := vec_first_non_NA(shared_it_timepoints), 
      by = amino_acid]
  arr[is.na(shared_it_timepoints), shared_it_timepoints := 0]
  print(arr[, table(shared_it_timepoints)])

  suppressWarnings(arr[, shared_pbl_timepoints := NULL])
  arr[!is.na(timepoint), 
      'shared_pbl_timepoints' := as.integer(sum(value > 0)),
      by = .(timepoint, amino_acid)]
  arr[, shared_pbl_timepoints := vec_first_non_NA(shared_pbl_timepoints), 
      by = amino_acid]
  arr[is.na(shared_pbl_timepoints), shared_pbl_timepoints := 0]
  print(arr[, table(shared_pbl_timepoints)])

  arr[timepoint %in% blood_timepoints, 
      'shared_pbl_timepoints' := factor(as.integer(sum(value > 0))),
      by = amino_acid]
  arr[, shared_pbl_timepoints := vec_first_non_NA(shared_pbl_timepoints), 
      by = amino_acid]
  arr[is.na(shared_pbl_timepoints), shared_pbl_timepoints := 0]

  ## Summarize all intratumoral timepoints into one string
  if (compute_pbl_it_timepoints) {
    suppressWarnings(invisible(arr[, it_timepoints := NULL]))
    if (F) {
      ## Faster but harder to parse
      system.time(arr[, it_timepoints := lapply(seq.int(.N),
                          function(x) list(as.list(unique(it_timepoint)))),
                      by = .(patient, amino_acid)])
    } else {
      system.time(arr[,
                  'it_timepoints' := paste(sort(unique(it_timepoint)),
                                           collapse = ' - '),
                  by = .(patient, amino_acid)])
    }
    arr[it_timepoints == 'NA', it_timepoints := NA]
    arr[is.na(it_timepoints) | it_timepoints == '', it_timepoints := 'None']
    it_timepoints_levs <- c('None', 'Baseline', 'Baseline - Post-induction',
                            'Baseline - On nivo',
                            'Baseline - Post-induction - On nivo',
                            'Post-induction', 'Post-induction - On nivo',
                            'On nivo')
    arr[, it_timepoints := factor(it_timepoints,
                                  levels = it_timepoints_levs)]
  }
  return(arr)
}


#' Finalize TCR data for parallel coordinate plotting
#'
#'
prepare_TCR_chrono <- function(patient = 'pat_11',
                               timepoint_v = 'timepoint',
                               facet_var = NULL,
                               colour_var = 'shared_timepoints',
                               sub_missing_with_zeros = T,
                               collapse_on_aa = T,
                               x_axis_tps = timepoints,
      highlight_function = function(arr) {
        arr[it_timepoint == 'On nivo' & freq_rank >= .8, hl_var := T]
      },
                               cluster_tcrs = T,
                               p_var = 'normalized_frequency',
                               compartment = 'tumor') {
  l_patient <- patient
  if (compartment == 'tumor') {
    arr <- read_adaptive_seqs(force_reload = F, patient = l_patient, 
                              collapse_on_aa = collapse_on_aa)
    if (null_dat(arr)) return(NULL)
  } else if (compartment == 'blood') {
    arr <- read_annotated_bloodadaptive_seqs(force_reload = F, 
                                             patient = l_patient,
                                             collapse_on_aa = collapse_on_aa)
  }

  if (null_dat(arr)) 
    return(NULL)
  if (arr[, uniqueN(timepoint)] == 1 && length(x_axis_tps) > 1) 
    return(NULL)

  ## Replace unobserved TCR/timepoint combinations with zero values
  if (sub_missing_with_zeros) {
    setkey(arr, amino_acid, timepoint)
    subs <- expand.grid('amino_acid' = arr[, unique(amino_acid)],
                        'timepoint' = arr[, unique(timepoint)]) %>%
      as.data.table
    # subs[, .N, by = amino_acid][N != length(timepoints)]
    ## Ensure each timepoints-TCR combination is represented in the data
    # arr[amino_acid %in% arr[, .N, by = amino_acid][N != 3, amino_acid]] %>%
    #   { .[order(amino_acid, timepoint)] } %>%
    #   { .[, .N == 0] } %>%
    #   stopifnot()
    arr <- arr[subs, ]
    arr[is.na(value), value := 0]
  }

  ## Unknown function of this code block
  if (F) {
    ## Turn this into resuable function when things are quiet
    comp_cols <- c('adaptive_sample_name')
    anchor_cols <- c('timepoint')
    sub_df <- unique(arr[!is.na(get(comp_cols))], by = anchor_cols) %>%
      { .[, c(comp_cols, anchor_cols), with = F] }
    setkeyv(sub_df, anchor_cols)
    comp_col_vals <- sub_df[arr[, get(anchor_cols)], get(comp_cols)]
    arr[is.na(get(comp_cols)), (comp_cols) := comp_col_vals]
  }

  ## Highlight TCRs of interest
  if (T || !is.null(highlight_function)) {
    arr[, 'hl_var' := F]
    arr <- highlight_function(arr)
    # arr[hl_var == T, .N, by = .(amino_acid, timepoint)][N > 1, amino_acid[1]] %>%
    #   { arr[amino_acid %in% .] }
    # arr[, .N, by = .(amino_acid, timepoint)][N > 1, amino_acid[1]] %>%
    #   { arr_u[amino_acid %in% .] }
    # arr_u <- unique(arr[order(-frank(hl_var))], 
    #        by = c('amino_acid', 'timepoint', 'hl_var'))
    arr <- unique(arr[order(-frank(hl_var))], 
                  by = c('amino_acid', 'timepoint'))
    # arr[hl_var == T, .N, by = .(amino_acid, timepoint)][N == 1, amino_acid[1]] %>%
    #   { arr[amino_acid %in% .] }
    # if (check_assumptions) {
    #   ## Select a recurrent amino_acid
    # }
    # arr <- arr[order(amino_acid, timepoint, hl_var)]
    # arr <- unique(arr, by = c('amino_acid', 'timepoint'))
  }

  ## Cluster TCRs such that identical TCRs (w.r.t. relevant variables) are
  ## represented by a single row
  if (cluster_tcrs) {
    ## Determine how to cluster TCRs into 'similar dynamic behaviour'
    if (compartment == 'tumor' ||
        (compartment == 'blood' && colour_var == 'shared_timepoints')) {
      wide_dat <- dcast(arr, amino_acid ~ timepoint, value.var = 'value')
    } else if (compartment == 'blood' && colour_var == 'it_timepoints') {
      wide_dat <- dcast(arr, amino_acid ~ timepoint + it_timepoints,
                        value.var = 'value')
    } else if (compartment == 'blood' && colour_var == 'hl_var') {
      arr <- highlight_function(arr)
      wide_dat <- dcast(arr, amino_acid ~ timepoint + hl_var, 
                        value.var = 'value')
    } else {
      stopf('Not implemented')
    }
    timepoints_present <- unique(arr[, unique(get(timepoint_v))])
    wide_dat_cols <- setdiff(colnames(wide_dat), 'amino_acid')
    wide_dat[, 'cluster_N' := .N, by = wide_dat_cols]
    cluster_assigns <- unique(wide_dat, by = wide_dat_cols) %>%
      { .[, 'cluster_ID' := 1:.N] }
    wide_dat <-
      controlled_merge(wide_dat,
                       cluster_assigns[, c('cluster_ID', wide_dat_cols),
                                       with = F])
    arr <- controlled_merge(arr,
                            wide_dat[, .(amino_acid, cluster_ID, cluster_N)])
    if (compartment == 'tumor') {
      arr <- unique(arr, by = c('cluster_ID', 'timepoint'))
    } else if (compartment == 'blood' && colour_var == 'it_timepoints') {
      arr <- unique(arr, by = c('cluster_ID', 'timepoint', 'it_timepoints'))
    } else if (compartment == 'blood' && colour_var == 'hl_var') {
      arr <- unique(arr, by = c('cluster_ID', 'timepoint', 'hl_var'))
    }
  }
  return(arr)
}
