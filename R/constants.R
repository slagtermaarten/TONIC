mut_load_var_types <- 
 c('Missense Variant', 'Frameshift Variant', 'Stop Gained', 
   'Conservative Inframe Insertion', 'Conservative Inframe Deletion', 
   'Disruptive Inframe Insertion', 'Disruptive Inframe Deletion',
   'Structural Interaction Variant', 
   # 'TF Binding Site Variant', 
   'Stop Lost', 
   # 'Start Lost', 
   'Protein Protein Contact', 'Stop Retained Variant') %>% sort
# tolower(paste(mut_load_var_types, collapse = ', '))

timepoints <- c('baseline' = 'Baseline', 'post.induction' = 'Post-induction',
                'on.nivo' = 'On nivo')
it_timepoints <- timepoints
timepoints_inv <- setNames(names(timepoints), timepoints)

## Labels for timepoints 1, 2 and 3
timepoint_labels <- c('baseline', 'post_induction', 'on_nivo')

timepoint_colors <- maartenutils::gen_color_vector('Zissou1', 1) %>%
  darken(factor = rev(c(.75, 1.0, 1.25))) %>%
  setNames(timepoints)

blood_timepoints <- auto_name(c(-2, 0, 6, 10, 12))
blood_timepoints_subs <- setNames(blood_timepoints, 
                                  c(timepoints, 'On nivo 2', 'On nivo 3'))
blood_timepoints_subs_inv <- setNames(c(timepoints, 'On nivo 2', 'On nivo 3'),
                                      blood_timepoints)
blood_timepoint_names <-
  setNames(paste0('PBL ', c('Baseline', 'Post-induction', 'On nivo', 
                  'On nivo + 4', 'On nivo + 6')), blood_timepoints)


## Marleen's ingrained way of ordering the projects
treatment_arms <- c('Radiotherapy', 'Doxorubicin', 'Cyclophosphamide',
                    'Cisplatin', 'No induction')
## Ordering by clinical benefit and NanoString results
treatment_arms <- c('No induction', 'Radiotherapy', 'Cyclophosphamide',
                    'Cisplatin', 'Doxorubicin')
label_reps <- c('efron_thisted_estimator' = 'TCRSeq repertoire size', 
                'adaptive_t_cells' = 'TCRSeq %T-cells',
                'sample_clonality' = 'TCRSeq clonality',
                'ca15_3' = 'Baseline CA15.3',
                's_til' = 'Baseline stromal TIL [%]',
                'pd_l1_immunoinfiltrate' = 'Baseline PD-L1 [%]')

## More human readable labels for plotting purposes
axis_subs <- c('tp1' = 'Baseline', 
               'tp2' = 'Post-induction',
               'tp3' = 'On-nivo',
               'baseline' = 'Baseline', 
               'post.induction' = 'Post-induction',
               'post_induction' = 'Post-induction',
               'on.nivo' = 'On-nivo',
               'on_nivo' = 'On-nivo')

arm_colors <- 
  maartenutils::gen_color_vector('Royal1', 5) %>%
  # `[`(c(3,4,5)) %>%
  { c(darken(maartenutils::gen_color_vector('Zissou1', 1), 
             factor = c(.9, 1.3)), .) } %>%
  `[`(c(3,1,2,5,4)) %>%
  setNames(treatment_arms) %>%
  attr_pass('class', 'color_vector')

## Color scheme developed by Marleen and Leonie
arm_colors <- c('No induction' =  '#ABAAAC', 'Radiotherapy' = '#B51533',
                'Cyclophosphamide' = '#009343', 'Cisplatin' = '#982593', 
                'Doxorubicin' = '#0996D8') %>%
  attr_pass('class', 'color_vector')

resp_colors <- maartenutils::gen_color_vector('Zissou1', 2) %>%
  darken(factor = c(1.0, 1.2)) %>%
  setNames(c('R', 'NR')) %>%
  attr_pass('class', 'color_vector')

resp_colors <- 
  maartenutils::darken(c(rev(gen_color_vector(name = 'Zissou1', n = 2)),
                         'grey70'),
                       c(1.2, 1.0, 1.0)) %>%
  setNames(c('R', 'NR', 'NA')) %>%
  attr_pass('class', 'color_vector')
# plot(resp_colors)

recist_colors <-
  c(rev(gen_color_vector(name = 'Zissou1', n = 2)), 'grey70') %>%
  { c(.[1], .[1], .[1], .[2], .[3]) } %>%
  maartenutils::darken(c(.7, 1.1, 1.3, 1.0, 1.0)) %>%
  setNames(c('CR', 'PR', 'SD', 'PD', 'NA')) %>%
  attr_pass('class', 'color_vector')

gene_colors <- setNames(c('white', 'grey20', "orangered3", "orangered3",
                          muted('blue'), muted('blue')),
                        c('chrom LOH', 'hemizygous loss', 'chrom loss', 
                          'chrom loss LOH', 'chrom gain', 'chrom gain LOH'))

tissue_colors <- gen_color_vector(name = 'Darjeeling1', n = 8) %>%
  { c(.[1:7], 'grey80') } %>%
  darken(c(rep(1.2, 7), 1)) %>%
  setNames(c('Lymphnode', 'Mamma', 'Muscle', 'Liver', 'Stomach', 'Skin', 'Lung',
             'NA')) %>%
  attr_pass('class', 'color_vector')

pam50_colors <- rev(gen_color_vector(name = 'Zissou1', n = 5)) %>%
  { c(.[1], arm_colors[3], arm_colors[3], .[4], .[5], 'grey80') } %>%
  setNames(c('Basal', 'LumA', 'LumB', 'Normal', 'Her2', 'NA')) %>%
  darken(c(1.2, .9, 1.2, 1.0, 1.4, 1)) %>%
  { .[c(1, 2, 3, 4, 5, 6)] } %>%
  attr_pass('class', 'color_vector')

darkred <- '#AC1200'

var_type_colors <- rev(gen_color_vector(name = 'Zissou1', n = 5)) %>%
  { c(.[1], .[1], 'grey50', .[4], .[5]) } %>%
  setNames(c('Frameshift Indel', 'Inframe Indel', 'Missense Variant', 
             'Stop Gained', 'Stop Lost')) %>%
  darken(c(1.1, 1.4, 1.0, 1.0, 1.4)) %>%
  # { .[c(1, 4, 5, 2, 3)] } %>%
  attr_pass('class', 'color_vector')
# structure(var_type_colors)

brca_pathology_colors <- c('white', darkred, darkred, 'grey80') %>%
  darken(c(1, 1/1.4, 1.0, 1)) %>%
  setNames(c('WT', 'BRCA1-mut', 'BRCA2-mut', 'NA')) %>%
  attr_pass('class', 'color_vector')

brcaness_colors <- { c('white', darkred, 'grey80') } %>%
  setNames(c('WT', 'BRCA1-like', 'NA')) %>%
  attr_pass('class', 'color_vector')

mut_effect_colors <- rev(gen_color_vector('Spectral', n = 4)) %>%
  setNames(c('NA', 'MODIFIER', 'MODERATE', 'HIGH')) %>%
  attr_pass('class', 'color_vector')

invisible(patient_labels[, 'pfs_binned' := cut(pfs_w, breaks = 4)])
patient_labels[, pfs_w := .SD[!is.na(pfs_w), pfs_w], by = patient]

tonic_color_palettes <- list('clinical_response' = resp_colors,
                             'recist' = recist_colors,
                             'pam50' = pam50_colors,
                             'tissue' = tissue_colors,
                             'brca_pathology' = brca_pathology_colors,
                             'var_type_colors' = var_type_colors,
                             'brcaness' = brcaness_colors,
                             'mut_effect' = mut_effect_colors,
                             'arm' = arm_colors,
                             'gene_classification' = gene_colors,
                             'pfs_binned' = rev(gen_color_vector('Spectral', 
                                 n = patient_labels[, uniqueN(pfs_binned)])),
                             'timepoint' = timepoint_colors)
