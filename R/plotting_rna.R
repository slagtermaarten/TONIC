compute_FC <- function(dtf) {
  if (nrow(dtf) < 2 || dtf[, any(is.na(value))]) return(NULL)
  return(dtf[timepoint == 'Post-induction', log2(value)] -
         dtf[timepoint == 'Baseline', log2(value)])
}


compute_intensity <- function(dtf) {
  if (nrow(dtf) < 2 || dtf[, any(is.na(value))]) return(NULL)
  return(.5 * (dtf[timepoint == 'Post-induction', log2(value)] +
         dtf[timepoint == 'Baseline', log2(value)]))
}


compute_RNA_sigs <- function(sig = 'tis', genes = NULL) {
  if (!exists('rna_read_counts')) {
    source(file.path(p_root, 'R', 'load_rna_dat.R'))
  }

  if (!is.null(sig) && is.null(genes[1])) {
    gene_set_genes <- danaher_signature_genes[signature_name == sig,
                                              as.character(gene_name)]
  } else {
    gene_set_genes <- genes
  }

  rna_sigs <-
    rna_read_counts[external_gene_id %in% gene_set_genes,
                    .('value' = median(value, na.rm = T)),
                    by = c('patient', 'timepoint')]
  if (nrow(rna_sigs) == 0)
    return(NULL)

  rna_sigs[, 'FC' := compute_FC(.SD), by = c('patient')]
  rna_sigs[, 'A' := compute_intensity(.SD), by = c('patient')]

  return(rna_sigs)
}


prepare_rna_sig_MA <- function(sig = 'tis', genes = c('TP53'),
                                  colour_var = 'ca15_3',
                                  numeric_transform = 'log10') {
  rna_sigs <- compute_RNA_sigs(sig = sig, genes = genes)

  rna_sigs <-
    merge(rna_sigs, clinical_annotation[, c('patient', colour_var), with = F],
          by = 'patient')

  if (rna_sigs[, class(get(colour_var))] %in% c('integer', 'numeric', 'double')) {
    if (numeric_transform == 'Z') {
      rna_sigs[, (colour_var) := z_transform(get(colour_var))]
    } else if (numeric_transform == 'log10') {
      rna_sigs[, (colour_var) := log10(get(colour_var) + 1)]
    }
  }
  return(rna_sigs)
}


#' Bland-altman plot compare two timepoints
#'
#'
plot_RNA_sig_MA <- function(sig = 'tis', genes = c('TP53'),
                            colour_var = 'ca15_3',
                            numeric_transform = 'log10') {
  rna_sigs <- prepare_rna_sig_MA(sig = sig, genes = genes,
                                 colour_var = colour_var,
                                 numeric_transform = numeric_transform)
  ggplot(rna_sigs, aes(x = log10(A), y = FC)) +
    aes_string(colour = colour_var) +
    geom_point() +
    scale_x_continuous(name = 'Baseline value (log10)') +
    scale_y_continuous(name = 'FC Post-induction vs. Baseline (log2)') +
    ggtitle(sig) +
    theme(legend.position = 'right')
}


correlate_nano_rna <- function(sig = 'tis', genes = NULL) {
  RNA_sigs <- compute_RNA_sigs(sig = sig, genes = genes)
  if (null_dat(RNA_sigs)) return(NULL)
  RNA_sigs[, value := log2(value + 1)]
  RNA_sigs[, FC := NULL]
  RNA_sigs[, A := NULL]
  RNA_sigs[, 'rna_q50' := median(value, na.rm = T), by = c('timepoint')]
  setnames(RNA_sigs, 'value', 'rna_value')

  ## Correlate signature
  if (!is.null(sig) && is.null(genes[1])) {
    nano_dat <- danaher_scores.m[variable %in% gsub('_', '.', sig),
                                  .(patient, timepoint, value)]
    setnames(nano_dat, 'value', 'nano_value')
    if (null_dat(nano_dat)) return(NULL)
    ## Consider using only the samples for which both data types are available
    nano_dat[, 'nano_q50' := median(nano_value, na.rm = T), by = c('timepoint')]
  } else {
    nano_dat <- 
      melt(exp_levels[gene_symbol %in% genes], id.vars = 'gene_symbol',
           variable.name = 'filename', value.name = 'nano_value') %>%
      merge(patient_labels[, c('filename', 'patient', 'timepoint'), with = F],
            by = 'filename')
    if (null_dat(nano_dat)) return(NULL)
    nano_dat[, filename := NULL]
    if (!is.null(nano_dat[, nano_value])) {
      tryCatch(nano_dat[, 'nano_value' := log2(nano_value)], 
               error = function(e) { print(e); browser() }) 
    }
    nano_dat[, 'nano_q50' := median(nano_value, na.rm = T), by = c('timepoint')]
  }
  p_dat <- merge(RNA_sigs, nano_dat, by = c('patient', 'timepoint'))

  if (nrow(p_dat) == 0) return(NULL)
  return(p_dat)
}


#' Compare RNASeq and NanoString data with scatter plot
#'
#'
compare_nano_rna_sigs <- function(sig = 'tis', genes = NULL) {
  p_dat <- correlate_nano_rna(sig = sig, genes = genes)
  cor_v <- p_dat[, cor(rna_value, nano_value,
                       method = 'pearson', use = 'complete')]

  p <- ggplot(p_dat, aes(x = rna_value, y = nano_value)) + 
    geom_point() +
    geom_smooth(method = 'lm') +
    scale_x_continuous(name = 'RNASeq signature score') +
    scale_y_continuous(name = 'NanoString signature score') +
    ggtitle(as.formula(sprintf('%s~(gene~count==%d)', sig,
      danaher_signature_genes[signature_name == sig, uniqueN(gene_name)])))

  p <- p + annotate('text', label = sprintf('Pearson~italic(R)^2==%.2f', cor_v),
    x = interpolate_in_gg_range(plot = p, axis = 'x', degree = .1),
    y = interpolate_in_gg_range(plot = p, axis = 'y', degree = .9),
    colour = 'grey40',
    hjust = 0, vjust = 1, size = 3, parse = T)

  return(p)
}
