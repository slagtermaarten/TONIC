library(grid)

caplist_tonic <- c(caplist_def, 'Gu-Trantien', 'IFNy', 'IFNg', 'TNBC', 'TIL',
                   'mm2', 'Ichao1', 'cells', 'cell', 'TIS', 'score', 's', 'ng',
                   'LDH', 'IFN', 'Th1', 'CRP', 'TIL', 'CIBERSORT', 'GC',
                   'gamma', 'LumA', 'LumB', 'APM', 'FastQ',
                   'PD-L1', 'MMR', 'loss', 'SASP', 'RNA')
if (exists('rna_read_counts_salmon')) {
  caplist_tonic %<>% { c(., rownames(rna_read_counts_salmon)) }
}
# simple_cap('ifny', caplist = caplist_tonic)
tonic_cap <- function(x, ...) {
  stopifnot(is.character(x) || is.factor(x))
  ## First some preprocessing before final capitalisation
  x <- gsub('-', '_', x) %>%
       gsub('\\.', ' ', .) %>%
       gsub('IFNG', 'IFNy', .) %>%
       gsub('ifn gamma', 'IFNy', .) %>%
       gsub('fqreads', 'FastQ reads', .) %>%
       gsub('GU_TRANTIEN', 'Gu-Trantien', .) %>%
       gsub('PD_L1', 'PD-L1', .) %>%
       gsub('pdl1', 'PD-L1', .)
  return(simple_cap(tolower(x), caplist = caplist_tonic, ...))
}

prep_es <- function(gsea_res,
                    fdr_thresh = .25,
                    subset_var = 'arm',
                    combined_lev = 'all_arms',
                    N_pathways = uniqueN(gsea_res$GeneSet)) {
  p_dat <- gsea_res

  ## Put combined cohort first in line
  # subset_var %in% colnames(p_dat)
  # dplyr::select_(p_dat, subset_var)
  # unlist(dplyr::select_(p_dat, 'abs'))
  if (!is.null(subset_var)) {
    if (subset_var %in% colnames(p_dat)) {
      flevels <- unique(p_dat$subset_var)
      idx <- which(flevels == combined_lev)
      flevels <- flevels[c(idx, setdiff(seq_along(flevels), 1))]
      p_dat$subset_var <- factor(p_dat$subset_var, levels = flevels)
    }
  }

  p_dat$es[p_dat$fdr > fdr_thresh] <- NA
  # p_dat$nes[p_dat$fdr > fdr_thresh] <- NA
  p_dat$fdr[p_dat$fdr > fdr_thresh] <- NA
  p_dat$GeneSet <- tonic_cap(p_dat$GeneSet)
  ## Order by amount of significant hits
  if (F) {
    p_dat <- p_dat %>% group_by(GeneSet) %>%
      mutate(nes_tally = sum(!is.na(nes))) %>%
      arrange(-nes_tally, fdr)
    gs_ordering <- unique(p_dat$GeneSet)
  } else if (F) {
    p_dat <- p_dat %>% group_by(GeneSet) %>%
      mutate(nes_sum = sum(nes)) %>%
      arrange(nes_sum, fdr)
    gs_ordering <- unique(p_dat$GeneSet)
  } else {
    if (all(c('timepoint', 'arm') %in% colnames(p_dat) == c(T, T))) {
      wide_dat <- dcast(p_dat, GeneSet ~ timepoint, value.var = 'nes')
    } else if (all(c('timepoint', 'arm') %in% colnames(p_dat) == c(F, T))) {
      wide_dat <- dcast(p_dat, GeneSet ~ arm, value.var = 'nes')
    }
    d_obj <- dist(wide_dat[, c(2:ncol(wide_dat)), with = F], method = 'euclidean')
    clust <- hclust(d_obj, method = 'ward.D2')
    gs_ordering <- unlist(wide_dat$GeneSet[rev(clust$order)])
  }
  p_dat$GeneSet <- factor(p_dat$GeneSet, levels = gs_ordering)
  # idx <- match(unique(p_dat$resp_exp), hr_resp_exp)
  # p_dat$resp_exp <- factor(p_dat$resp_exp,
  #                          labels = hr_resp_exp[idx],
  #                          levels = names(hr_resp_exp)[idx])
  ## Only plot top N_pathways gene sets
  p_dat <- p_dat[p_dat$GeneSet %in% levels(p_dat$GeneSet)[1:N_pathways], ]
  return(p_dat)
}


plot_es <- function(p_dat,
                    N_pathways = length(unique(p_dat$GeneSet)),
                    ptitle = NULL,
                    fdr_thresh = .25,
                    # legend.position = c(.85, .55),
                    legend.position = 'bottom',
                    subset_var = 'arm',
                    fn = NULL,
                    x_var = 'resp_exp',
                    combined_lev = 'all_arms',
                    nes_range = NULL) {
  if (null_dat(p_dat)) return(NULL)
  p_dat <- prep_es(gsea_res = p_dat,
                   N_pathways = N_pathways,
                   subset_var = subset_var,
                   combined_lev = combined_lev,
                   fdr_thresh = fdr_thresh)
  if (is.null(nes_range)) {
    nes_range <- range(p_dat$nes)
    max_val <- max(ceiling(abs(nes_range)))
    nes_range <- c(-max_val, max_val)
  }

  p <- ggplot(p_dat, aes_string(y = 'GeneSet', x = x_var, fill = 'nes',
                                size = '1/fdr')) +
    geom_tile(size = 0) +
    geom_point(color = 'black') +
    rotate_x_labels(rotate_labels = 90) +
    scale_fill_gradient2(name = 'Normalized\nEnrichment\nScore',
                         low = muted('blue'),
                         high = muted('red'),
                         na.value = 'white',
                         # low = 'gold1',
                         # high = 'red',
                         limits = nes_range) +
    scale_y_discrete(name = '', expand = c(0, 0)) +
    scale_x_discrete(name = '', expand = c(0, 0))

  size_labelling <- function(x) round(1/x, 2)
  size_labelling(p_dat$fdr)

  ## Prevent errors this way
  if (!all(is.na(p_dat$fdr))) {
    p <- p + scale_size_continuous(name = 'FDR',
                                   range = c(0, 1),
                                   # trans = 'probability',
                                   # trans = 'log10',
                                   # limits = c(1e-5, fdr_thresh),
                                   # labels = function(x) x) +
                                   labels = size_labelling)
    # p <- p + guides(fill = guide_legend(order = 1), 
    #                 size = guide_legend(order = 2))
  }

  if (!is.null(subset_var) && subset_var %in% colnames(p_dat)) {
    p <- p + facet_wrap(as.formula(sprintf('~ %s', subset_var)),
                        ncol = uniqueN(dplyr::select_(p_dat, subset_var)))
  }

  legend_direction <- ifelse(legend.position == 'bottom', 
                             'horizontal', 'vertical')
  p <- p + theme(legend.position = legend.position,
                 legend.direction = legend_direction,
                 strip.text = element_text(size = 8),
                 legend.key.size = grid::unit(5, 'mm'),
                 axis.ticks = element_line(size = 0),
                 legend.background = element_blank())
  if (F) {
    p <- p + theme(axis.text.x = element_blank())
  } else {
    p <- p + rotate_x_labels(45)
  }

  if (!is.null(ptitle)) {
    ## Number of blocks per cm
    NR <- 3
    datef <- gsub('/', '_', format(Sys.time(), "%x"))
    date_dir <- file.path(img_dir, datef)
    dir.create(date_dir, showWarnings = F)

    fn <- sprintf('%s/GSEA_%s.pdf', date_dir, ptitle)
    mult <- ifelse(!is.null(subset_var) && 
                   subset_var %in% colnames(p_dat), 5, 1)
    pwidth <- uniqueN(p_dat$resp_exp)/NR
    pheight <- uniqueN(p_dat$GeneSet)/NR
    fw <- mult * pwidth + 10
    fh <- max(8, pheight + 4)
    if (F) {
      ## Correct width for amount of panels
      mult <- unit(mult, 'cm')
      pwidth <- unit(pwidth, 'cm')
      pheight <- unit(height, 'cm')
      p_grid <- set_panel_size(p, width = pwidth, height = pheight)
      fw <- unit(fw, 'cm')
      fh <- unit(fh, 'cm')
      pdf(file = fn, 
          width = convertUnit(fw, 'in'), 
          height = convertUnit(fh, 'in'))
      grid.draw(p_grid)
      dev.off()
    } else {
      if (!dir.exists(dirname(fn))) dir.create(dirname(filename),
                                                     recursive = T)
      ggsave(plot = p, 
             filename = fn, width = fw, height = fh, units = 'cm',
             fillOddEven = T,
             compress = F,
             useKerning = F, useDingbats = F)
      # ggsave(plot = p, 
      #        filename = sprintf('%s/GSEA_%s.eps', date_dir, ptitle), 
      #        width = fw, height = fh, units = 'cm')
    }
  }
  return(p)
}


plot_leading_edge_heatmap <- function(res, gsea_name = 'paired',
                                      timepoint = timepoints[1:2],
                                      fdr_thresh = .25, ...) {
  lapply(res[, unique(arm)], function(l_arm) {
    leading_edge <- res[arm == l_arm & fdr <= fdr_thresh,
                        setNames(leading_edge_genes, GeneSet)]
    if (is.null(leading_edge) || length(leading_edge) == 0) return(NULL)
    gs <- strsplit(x = leading_edge, split = '\\;')
    if (is.null(gs) || length(gs) == 0) return(NULL)
    sapply(seq_along(gs), function(idx) {
      plot_gene_set(gs_name = sprintf('%s_tp_comp-leading_edge_%s_%s',
                                      gsea_name, l_arm,
                                      names(gs)[idx]),
                    gene_symbols = gs[idx],
                    exp_mat = tpms_salmon,
                    # timepoints = 'Baseline')
                    # timepoints = 'On nivo')
                    timepoints = timepoint, ...)
    })
  })
}


#' Combine two timepoints and only plot all arms combined, not the individual
#' arms
#'
#'
plot_cr_gsea_res <- function(ptitle = 'unpaired',
                             bl_fn = 'rds2/GSEA_unpaired_baseline_TMM.rds',
                             pi_fn = 'rds2/GSEA_unpaired_postinduction_TMM.rds') {
  res_1 <- readRDS(file = bl_fn)
  res_2 <- readRDS(file = pi_fn)
  coln <- colnames(res_2)
  p_dat <- rbind(cbind(res_1[arm == 'All arms', coln, with = F], 
                       'timepoint' = 'Baseline'), 
                 cbind(res_2[arm == 'All arms', coln, with = F], 
                       'timepoint' = 'Post-induction'))
  plot_es(p_dat, legend.position = 'right', ptitle = ptitle,
          subset_var = NULL, x_var = 'timepoint')
  invisible()
}
