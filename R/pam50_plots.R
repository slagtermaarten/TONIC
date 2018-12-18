if (F) {
  l_patient_ordering <- rna_patient_ordering
} else {
  l_patient_ordering <- patient_ordering
}

subtype_predictions <-
  readRDS('rds/subtype_predictions.rds') %>%
  { .$subtype.proba } %>%
  as.data.frame %>%
  rownames_to_column('cf_number') %>%
  melt(id.vars = 'cf_number') %>%
  controlled_merge(rna_sample_annotation[, .(cf_number, patient, timepoint)],
                   by_cols = 'cf_number') %>%
  dplyr::mutate(variable = factor(variable, levels = names(tonic_color_palettes[['pam50']]))) %>%
  dplyr::mutate(patient = factor(patient, levels = l_patient_ordering)) %>%
  { .[!is.na(patient)] } %>%
  setkey(patient, timepoint) %>%
  { .[expand.grid(l_patient_ordering, timepoints)] } %>%
  { .[is.na(variable), variable := "NA"] } %>%
  { .[is.na(value), value := 1] }


bar_border_size <- 0

pam50_track_BL <- subtype_predictions %>%
  dplyr::filter(timepoint == 'Baseline') %>%
  ggplot(aes(x = patient, fill = variable, y = value)) +
  geom_col(color = 'black', size = bar_border_size, width = .8) +
  scale_fill_manual(name = 'PAM50', values = tonic_color_palettes[['pam50']]) +
  # track_style +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0), drop = F)

pam50_track_PI <- subtype_predictions %>%
  dplyr::filter(timepoint == 'Post-induction') %>%
  ggplot(aes(x = patient, fill = variable, y = value)) +
  geom_col(color = 'black', size = bar_border_size, width = .8) +
  scale_fill_manual(name = 'PAM50', values = tonic_color_palettes[['pam50']]) +
  # track_style +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0), drop = F)

nanostring_pam50 <- danaher_scores.m %>%
  dplyr::filter(variable %in% c('luma.cor', 'lumb.cor', 'basal.cor',
                                'her2.cor')) %>%
  dplyr::filter(!is.na(value)) %>%
  dplyr::mutate(variable = plyr::mapvalues(variable,
                              c('luma.cor', 'lumb.cor',
                                'basal.cor', 'her2.cor'),
                              c('LumA', 'LumB',
                                'Basal', 'Her2'))) %>%
  dplyr::group_by(patient, timepoint) %>%
  dplyr::mutate(value = value / sum(value)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(patient = factor(patient, levels = l_patient_ordering)) %>%
  dplyr::arrange(patient) %>%
  { .[!is.na(patient)] } %>%
  as.data.table %>%
  setkey(patient, timepoint) %>%
  { .[expand.grid(l_patient_ordering, timepoints)] } %>%
  { .[is.na(variable), variable := "NA"] } %>%
  { .[is.na(value), value := 1] }

# nanostring_pam50[, .N, by = .(patient, timepoint)]
# table(subtype_predictions$patient)
# table(nanostring_pam50$patient)

pam50_nano_track_BL <- nanostring_pam50 %>%
  dplyr::filter(timepoint == 'Baseline') %>%
  ggplot(aes(x = patient, fill = variable, y = value)) +
  geom_col(color = 'black', size = bar_border_size, width = .8) +
  scale_fill_manual(name = 'PAM50', values = tonic_color_palettes[['pam50']]) +
  # track_style +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0), drop = F)

pam50_nano_track_PI <- nanostring_pam50 %>%
  dplyr::filter(timepoint == 'Post-induction') %>%
  ggplot(aes(x = patient, fill = variable, y = value)) +
  geom_col(color = 'black', size = bar_border_size, width = .8) +
  scale_fill_manual(name = 'PAM50', values = tonic_color_palettes[['pam50']]) +
  # track_style +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0), drop = F)

pam50_nano_track_ON <- nanostring_pam50 %>%
  dplyr::filter(timepoint == 'On nivo') %>%
  ggplot(aes(x = patient, fill = variable, y = value)) +
  geom_col(color = 'black', size = bar_border_size, width = .8) +
  scale_fill_manual(name = 'PAM50', values = tonic_color_palettes[['pam50']]) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0), drop = F)


pam50_nano_track_BL_single <- nanostring_pam50 %>%
  dplyr::filter(timepoint == 'Baseline') %>%
  copy %>%
  { .[, y_var := 'PAM50'] } %>%
  ggplot(aes(x = patient, y = y_var, fill = variable)) +
  geom_rounded_tile(color = 'black') +
  scale_x_discrete(expand = c(0, 0), drop = F) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_manual(name = 'PAM50',
                    values = tonic_color_palettes[['pam50']])

if (F) {
  ## Compare PAM50 subtypings, draw separate plot for that
  plots <- list(rna_cr_track + track_style,
                pam50_track_BL + track_style + theme(legend.position = 'none'),
                pam50_track_PI + track_style + theme(legend.position = 'none'),
                pam50_nano_track_BL + track_style + theme(legend.position = 'none'),
                pam50_nano_track_PI + track_style + theme(legend.position = 'none'),
                pam50_nano_track_ON + track_style + theme(legend.position = 'none'))

  sth <- .6
  stopifnot(length(heights <- c(sth, 3, 3, 3, 3, 3)) ==
            length(plots))

  pam50_main <- ggpubr::ggarrange(plotlist = plots,
                                  ncol = 1,
                                  nrow = length(plots),
                                  heights = heights,
                                  # align = "hv",
                                  legend = 'none')
  print(pam50_main)
}
# grid.draw(get_legend(pam50_track_BL))
