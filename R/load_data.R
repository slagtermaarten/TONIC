## Load response data {{{
harmonize_celllines <- function(vec) {
  gsub('.CL_Untreated|\\.|-', '', vec)
}

fh <- fread(file.path('data-raw', 'cell_surface_data.csv'))

setnames(fh, c('cell_line', 
               paste(rep(c('MHC', 'PDL1'), each = 4), 
                     c('US', '0', '1', '10'), sep = "_")))
fh <- fh[-1]
num_cols <- setdiff(colnames(fh), 'cell_line')
fh[, (num_cols) := lapply(.SD, function(x) as.numeric(gsub(',', '.', x))), 
   .SDcols = num_cols]
setkey(fh, cell_line)

fh[, 'PDL1_FC_log2' := log2(PDL1_10) - log2(PDL1_0)]
fh[, 'PDL1_FC' := 2^PDL1_FC_log2]
# fh[, 'geom_mean' := (PDL1_0 * PDL1_10)^(1/2)]
# fh[, 'harm_mean' := harmonic_mean(PDL1_0, PDL1_10)]
fh[, cell_line := harmonize_celllines(cell_line)]
## Global response variable to use in regression analyses
resp_var <- 'PDL1_FC'
## Load response data }}}

## Load_transcriptome
dat_path <- file.path('data-raw', 'RNAseq_data', 'DESeq2_normalized.RData')
if (!exists('Sample_annotation') || !exists('RC_data')) 
  load(dat_path)
idx <- (is.na(Sample_annotation$Hours) | Sample_annotation$Hours == 24) & 
       Sample_annotation$Treatment %in% c('Untreated')
soi <- as.character(Sample_annotation$File_name[idx])
Sample_names <- paste(Sample_annotation$Cell_line[idx], 
                      Sample_annotation$Treatment[idx], sep = '_')

## Subselect to samples of interest and aggregate expression levels
gene_var <- 'ensembl_gene_id'
gene_var <- 'external_gene_id'
RC_f <- RC_data[, c(gene_var, soi)]
colnames(RC_f) <- c(gene_var, Sample_names)
# rownames(RC_f) <- RC_f[, gene_var]
if (any(table(RC_f[, gene_var]) > 1)) {
  RC_f <- aggregate(formula = as.formula(sprintf('. ~ %s', gene_var)), 
                    data = RC_f, 
                    FUN = sum)
  ## Some rows should have been combined and some columns unselected
  stopifnot(all(dim(RC_f) < dim(RC_data)))
}
rownames(RC_f) <- RC_f[, gene_var]

## Normalize expression data 
if (F) {
  ## Manual method
  col_facs <- 1e-6 * apply(RC_f[, Sample_names], 2, sum)
  RC_f_n <- scale(RC_f[, Sample_names], center = F, scale = col_facs) 
  ## All column sums should equal 1e6
  stopifnot(all(eps(apply(RC_f_n[, Sample_names], 2, sum), 1e6)))
  ## No rows should have been lost
  stopifnot(nrow(RC_f_n) == nrow(RC_f))
  ## Name rows and columns
  rownames(RC_f_n) <- as.character(RC_f[, gene_var])
} else if (F) { 
  ## Normalize with DESeq
  ## Data already in CPM format?
  colnames(RC_f)
  rownames(RC_f)
  gexp_counts <- RC_f
  gexp_counts <- gexp_counts[, colSums(gexp_counts) > nrow(gexp_counts)]
  print(dim(gexp_counts))
} else {
  RC_f_n <- RC_f
}

colnames(RC_f_n) <- harmonize_celllines(colnames(RC_f_n))
if (F) {
  inspect_mat(RC_f_n)
  setdiff(harmonize_celllines(colnames(RC_f_n)[-1]),
                              harmonize_celllines(fh[, cell_line]))
}

## Remove genes for which variance is 0, this confuses ggsea code
idx <- which(!eps(apply(RC_f_n[, setdiff(colnames(RC_f_n), gene_var)], 
                        1, sd), 1e-6))
message(
  format_frac(msg = 'RNA exp data genes removed with zero std', 
              num = nrow(RC_f_n) - length(idx),
              denom = nrow(RC_f_n)))
RC_f_n <- RC_f_n[idx, ]
rm(idx)

## Handle technical duplicates
if (F) {
  ## Take mean of technical replicates
  sample_names <- gsub('(.+)\\.1', '\\1', colnames(RC_f_n))
  tech_means <- sapply(auto_name(unique(sample_names)), function(sn) {
    idx <- which(sample_names == sn)
    means <- apply(RC_f_n[, idx], 1, mean)
    return(means)
  })
  RC_f_n <- tech_means
  rm(sample_names)
  rm(tech_means)
} else {
  ## Restrict data to experiments conducted in TS lab
  sample_names <- gsub('(.+)\\.1', '\\1', colnames(RC_f_n))
  ## Invert the boolean to select DP experiments instead
  RC_f_n <- RC_f_n[, !duplicated(sample_names)]
}

## Select down to cell lines for which both expression and PD-L1 staining is
## available
shared_celllines <- intersect(harmonize_celllines(colnames(RC_f_n)[-1]),
                              harmonize_celllines(fh[, cell_line]))
RC_f_n <- RC_f_n[, colnames(RC_f_n) %in% shared_celllines]
setkey(fh, cell_line)
fh <- fh[shared_celllines]
stopifnot(all(colnames(RC_f_n) == fh[, cell_line]))