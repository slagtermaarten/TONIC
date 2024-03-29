knitr::opts_chunk$set(message = FALSE, cache = T, cache.lazy = F,
                      cache.comments = F, autodep = T, warning = FALSE,
                      results = 'hide', fig.keep = 'high',
                      echo = TRUE, error = FALSE)
knitr::opts_knit$set(progress = TRUE, verbose = TRUE)
# knitr::clean_cache()
# knitr::opts_chunk$get("cache.path")
# knitr::knit_global()
# setwd('~/TONIC')
p_root <- '~/TONIC'
source(file.path(p_root, 'R', 'init.R'))
# read_adaptive_seqs()
source(file.path('~/libs', 'ALICE', 'ALICE.R'))
ALICE_input <- file.path(p_root, 'ALICE', 'ALICE_input')
dir.create(ALICE_input, showWarnings = F)

ALICE_lib <- file.path('~/libs', 'ALICE')
load(file.path(ALICE_lib, 'VDJT.rda'))
# print(naturalsort(unique(VDJT[, 1])))
# print(naturalsort(unique(VDJT[, 2])))

N_cores <- 24
library('doParallel')
doParallel::registerDoParallel(cores = N_cores)

format_VJ <- function(name = 'TCRBJ02-01') {
  vapply(name, function(n) {
    if (grepl('unresolved', n)) return('NA')
    pre <- gsub('^TCRB(V|J)(\\d{1,2})-(\\d{1,2})', 'TRB\\1', n)
    ## Get rid of leading zeros
    major <- as.integer(gsub('^TCRB(V|J)(\\d{1,2})-(\\d{1,2})', '\\2', n))
    minor <- as.integer(gsub('^TCRB(V|J)(\\d{1,2})-(\\d{1,2})', '\\3', n))
    ## Reconstruct TCR
    res <- sprintf('%s%d-%d', pre, major, minor)
    if (res %nin% c(unique(VDJT[, 1]), unique(VDJT[, 2]))) {
      res_c <- sprintf('%s%d', pre, major)
      if (res_c %nin% c(unique(VDJT[, 1]), unique(VDJT[, 2]))) {
        # warningf('%s really not in VDJT', res)
        res <- 'NA'
      } else {
        res <- res_c
      }
    }
    return(res)
  }, character(1))
}


run_ALICE <- function(patient = 'pat_1', redo_input = F, redo_output = F,
                      tps = blood_timepoints, analysis_name = '',
                      cores = 2,
                      combine_tps = F) {
  output_dir <- file.path(p_root, 'ALICE',
                          sprintf('%s%s', patient, analysis_name))
  o_fn <- file.path(output_dir, sprintf('%s_output%s.rds', patient,
                                        analysis_name))
  if (!redo_output && file.exists(o_fn))
     return(readRDS(o_fn))

  i_fn <- file.path(ALICE_input, sprintf('%s_i_list%s.rds', patient,
                                         analysis_name))

  if (redo_input || !file.exists(i_fn)) {
    arr <- read_adaptive_seqs(force_reload = F, patient = patient,
                              collapse_on_aa = T,
                              p_var = 'normalized_frequency')
    av_tps <- arr[, auto_name(intersect(unique(timepoint), tps))]
    av_tps %<>%
      {
        .[names(.) %in% blood_timepoints] <- blood_timepoint_names[.];
        .
      }
    av_tps_cn <- variabilize_character(av_tps)

    i_list <-
      plyr::llply(names(av_tps), function(tp) {
        arr[timepoint == tp, .('CDR3.amino.acid.sequence' = amino_acid,
                               timepoint = tp,
                               v_gene, j_gene,
                               'bestVGene' = format_VJ(v_gene),
                               'bestJGene' = format_VJ(j_gene),
                               'Read.count' = read_count)] %>%
        { .[!is.na(bestVGene) & !is.na(bestJGene)] }
      }, .parallel = T)

    
    if (combine_tps) {
      i_list_c <- i_list[[1]]
      setnames(i_list_c, 'Read.count', av_tps_cn[1])
      i_list_c[, (av_tps_cn[1]) := replace_na(get(av_tps_cn[1]), 0)]
      if (length(i_list) > 1) {
        for (i in 2:length(i_list)) {
          i_list_c <- dt_m <- merge(i_list_c, i_list[[i]], all = TRUE,
                                    by = c('CDR3.amino.acid.sequence', 'v_gene', 
                                           'j_gene', 'bestVGene', 'bestJGene'))
          i_list_c <- i_list_c[, -grep('timepoint', colnames(i_list_c)), 
                               with = F]
          setnames(i_list_c, 'Read.count', av_tps_cn[i])
        }
      }
      ## Replace NA CDR3 counts with 0
      for (tp in av_tps_cn) {
        i_list_c[, (tp) := replace_na(get(tp), 0)]
      }
      ## Compute total read counts to be used in clustering
      i_list_c$Read.count <- apply(i_list_c[, av_tps_cn, with = F], 1, sum, 
                                   na.rm = T)
      i_list <- list(i_list_c)
    }
    saveRDS(i_list, i_fn)
  } else {
    i_list <- readRDS(i_fn)
  }

  messagef('Starting ALICE for %s %s', patient, analysis_name)
  res <- ALICE_pipeline(DTlist = i_list, folder = output_dir,
                        cores = cores, iter = 10, nrec = 5e5)
  messagef('Finished ALICE for %s %s', patient, analysis_name)
  saveRDS(res, o_fn)
}

plyr::l_ply(patient_labels[!is.na(adaptive_sample_name),
                                  naturalsort(unique(patient))],
  function(patient) {
    tryCatch(run_ALICE(patient, tps = c(blood_timepoints, it_timepoints),
                       combine_tps = T,
                       cores = 2,
                       redo_input = F,
                       redo_output = F,
                       analysis_name = '_all_timepoints'),
             error = function(e) {
               print(e);
               return(NULL)
             })
  }, .parallel = F)

# plyr::l_ply(blood_adaptive[, naturalsort(unique(patient))],
#   function(patient) {
#     tryCatch(run_ALICE(patient, tps = blood_timepoints, 
#                        cores = 1,
#                        analysis_name = ''),
#              error = function(e) {
#                print(e);
#                return(NULL)
#              })
#   }, .parallel = F)
 
plyr::l_ply(patient_labels[!is.na(adaptive_sample_name),
                                  naturalsort(unique(patient))],
  function(patient) {
    tryCatch(run_ALICE(patient, tps = it_timepoints, analysis_name = 'TIL',
                       cores = 2,
                       redo_output = F),
             error = function(e) {
               print(e);
               return(NULL)
             })
  }, .parallel = F)

