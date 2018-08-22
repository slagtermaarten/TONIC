gen_oncoplot <- function(dtf) {
  setkey(dtf, patient)
  dtf %<>% { .[as.character(p_dat$patient)] }
  stopifnot(unique(dtf$patient) == p_dat$patient)
  p <- ggplot(dtf[!is.na(variant_classification)], 
                      aes(x = patient, y = hugo_symbol, fill = CNt)) +
    geom_tile(color = 'black') +
    geom_point(data = dtf[!is.na(variant_classification) & 
                        grepl('LOH', variant_classification)],
               color = 'white') +
    # scale_fill_gradient2(name = 'Normalized\ncopy number',
    #                      na.value = 'grey80',
    #                      midpoint = 1,
    #                      # high = scales::muted('green'),
    #                      high = 'forestgreen',
    #                      # high = 'gold1',
    #                      low = scales::muted('red')) +
    #                      # limits = range(MMR_dat$CNt_n, na.rm = T) %>%
    #                      # { c(floor(.[1]), ceiling(.[2])) }) +
    #                      # limits = range(MMR_dat$CNt_n, na.rm = T)) +
    track_style + 
    theme(axis.text.y = element_text(size = 6))
  return(p)
}
