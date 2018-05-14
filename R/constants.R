timepoints <- c('baseline' = 'Baseline', 'post.induction' = 'Post-induction',
                'on.nivo' = 'On nivo')
timepoints_inv <- setNames(names(timepoints), timepoints)
blood_timepoints <- auto_name(c(-2, 0, 6, 10, 12))
## Marleen's ingrained way of ordering the projects
treatment_arms <- c('Radiotherapy', 'Doxorubicin', 'Cyclophosphamide',
                    'Cisplatin', 'No induction')
## Ordering by clinical benefit and NanoString results
treatment_arms <- c('No induction', 'Radiotherapy', 'Cyclophosphamide',
                    'Cisplatin', 'Doxorubicin')
label_reps <- c('efron_thisted_estimator' = 'Repertoire size', 
                'adaptive_t_cells' = 'Adaptive %T-cells',
                'sample_clonality' = 'TCR clonality')

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
  setNames(c('R', 'NR'))

timepoint_colors <- maartenutils::gen_color_vector('Zissou1', 1) %>%
  darken(factor = rev(c(.75, 1.0, 1.25))) %>%
  setNames(timepoints)

tonic_color_palettes <- list('clinical_response' = resp_colors,
                             'arm' = arm_colors,
                             'timepoint' = timepoint_colors)
