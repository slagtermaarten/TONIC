library(grid)
source('R/geom_rounded_tile.R')

gen_oncoplot <- function(dtf) {
  setkey(dtf, patient, hugo_symbol)
  dtf <- dtf[expand.grid(dtf[, unique(patient)], 
                         dtf[, unique(hugo_symbol)])]
  setkey(dtf, patient)
  dtf %<>% { .[as.character(oncop_dat$patient)] } %>%
    dplyr::mutate(patient = factor(patient, levels = levels(oncop_dat$patient)))

  ## Format chromosomal data
  chrom_data <- dtf[grepl('chrom', variant_classification)] 
  chrom_data[, LOH := grepl('LOH', variant_classification)]
  chrom_data[, variant_classification := gsub('chrom ', '', variant_classification)]
  chrom_data[, variant_classification := gsub('LOH', '', variant_classification)]
  chrom_data[, variant_classification := gsub('\\s*$', '', variant_classification)]
  chrom_data[, variant_classification := simple_cap(variant_classification)]
  chrom_data[is.na(variant_classification) | variant_classification == '', 
             variant_classification := 'Not altered']
  # chrom_data[variant_classification == '', variant_classification := NA]
  chrom_data <- 
    clean_columns('', chrom_data, c('arm', 'clinical_response', 'gene_id'))

  mut_data <- dtf[!grepl('chrom', variant_classification) &
                  variant_classification != ''] %>%
    dplyr::filter(!is.na(hugo_symbol))

  ## Preprocess shapes
  unique(mut_data$variant_classification)
  shapes <- c('Missense Variant' = 21, 'Frameshift Variant' = 22, 
              'Upstream Gene Variant' = 23)
  missing_shapes <- setdiff(dtf[, unique(variant_classification)], 
                            names(shapes)) %>% setdiff(NA)
  shapes <- c(shapes, setNames(22 + seq_along(missing_shapes), missing_shapes))
  shapes <- c(shapes, setNames(rep(NA, length(missing_shapes)), missing_shapes))
  
  p <- ggplot(chrom_data, aes(x = patient, y = hugo_symbol, 
              fill = variant_classification)) +
    # geom_tile(color = NA) +
    geom_rounded_tile(color = 'black') +
    geom_rounded_tile(data = chrom_data[LOH == T], alpha = .8, fill = NA, 
                      size = .3,
                      roundness = 1,
                      color = 'black') +
    # geom_point(data = chrom_data[LOH == T], alpha = .8, fill = 'white', 
    #            color = 'white', pch = 22, size = 3, show.legend = F) +
    geom_point(data = mut_data, 
               mapping = aes(x = patient, y = hugo_symbol,
                             shape = variant_classification),
               color = 'black', fill = 'black', size = 1, inherit.aes = F) +
    scale_shape_manual(name = 'SNV type', labels = simple_cap, 
                       values = shapes) +
    scale_x_discrete(name = '', drop = F, expand = c(0, 0)) +
    scale_y_discrete(name = '', drop = T, expand = c(0, 0)) +
    scale_fill_manual(name = 'Copy number status', 
      na.value = 'white',
      values = c('Gain' = tonic_color_palettes[['clinical_response']][['R']], 
                 'Loss' = tonic_color_palettes[['clinical_response']][['NR']],
                 'Not altered' = 'grey80')) +
    # rotate_x_labels(90) +
    track_style + 
    theme(axis.text.y = element_text(size = 6))

  return(p)
}

# B2M_track <- gen_oncoplot(rbind(MMR_dat, B2M_gene_dat, fill = T))
# MMR_dat[patient == 'pat_26']
# options(max.print = 2000)
# MMR_dat[hugo_symbol == 'POLE']
# print(B2M_track)

if (F) {
  data(iris)
  ggplot(iris, aes(x = Sepal.Length, y = Petal.Length)) + geom_rounded_square()

  if (F) {
    g <- ggplotGrob(B2M_track)
    str(g$grobs[[6]])
    g$grobs[[6]]$children[[3]]$gp$r <- unit(0.5, 'snpc')
    g$grobs[[6]]$children[[3]]$gp$r <- 0.5
    grid.draw(g)
  }
}
