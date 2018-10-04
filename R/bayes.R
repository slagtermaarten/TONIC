pacman::p_load(shinystan)
pacman::p_load("rstan")
if (!require(brms)) {
  devtools::install_github("paul-buerkner/brms")
}
library(brms)

brms_settings <- list('family' = bernoulli(),
                      'warmup' = 1000, 'iter' = 3e3, 'chains' = 4,
                      'control' = list(adapt_delta = 0.95))

brms_fc_settings <- list('family' = gaussian(),
                         'warmup' = 1000, 'iter' = 3e3, 'chains' = 4,
                         'control' = list(adapt_delta = 0.95))

run_brm <- function(p_dat, instance_name, save_RDS = T) {
  brm_tonic <- brm(formula = value ~ 1 + (timepoint | arm),
                   data = p_dat, family = gaussian(),
                   prior = NULL,
                   warmup = 1000, iter = 4000, chains = 5,
                   control = list(adapt_delta = 0.95))
  messagef('finished: %s, size: %s', instance_name,
           format(object.size(brm_tonic), units = 'Mb'))
  if (save_RDS) {
    if (!dir.exists(file.path(p_root, 'rds')))
      dir.create(file.path(p_root, 'rds'))
    saveRDS(brm_tonic,
            file.path(p_root, 'rds', sprintf('%s.rds', instance_name)))
  }
  return(brm_tonic)
}


run_brm_gene <- function(gene_symbol = 'CD274', save_RDS = T) {
  p_dat <- prep_gene_parallel_coords(gene = gene_symbol)[['p_dat']]
  run_brm(p_dat = p_dat, instance_name = gene_symbol, save_RDS = save_RDS)
}


run_brm_geneset <- function(geneset = 'apm', save_RDS = T) {
  p_dat <- danaher_scores.m[variable == geneset][!is.na(value)]
  run_brm(p_dat = p_dat, instance_name = geneset, save_RDS = save_RDS)
}


#' Extract parameters of interest from object
#'
#'
extract_params <- function(brms_object = brm_tonic,
                           params, burnin = 1000) {
  stopifnot(class(brms_object) == 'brmsfit')

  res <- lapply(auto_name(params), function(param) {
    ## Concatenate values of all chains, discarding burn-in samples
    Reduce(c, lapply(brms_object$fit@sim$samples, function(sim) {
      vec <- sim[[param]]
      return(vec[burnin:length(vec)])
    }))
  })

  if (!is.null(names(params)) && any(names(params) != params)) {
    names(res) <- names(params)
  }

  return(res)
}


compare_params <- function(brms_object = brm_tonic,
                           p1 = 'r_induction_therapy[Doxorubicin,Intercept]',
                           p2 = 'r_induction_therapy[Cisplatin,Intercept]') {
  print(sprintf('%s vs %s', p1, p2))
  extract_params(brms_object = brms_object, params = p1) -
  extract_params(brms_object = brms_object, params = p2)
}


#' Perform some comparisons of interest on brms object
#'
#'
perform_model_comps <- function(brms_object = brm_tonic) {
  print(brms_object$fit)
  compare_params(brms_object = brms_object,
                 p1 = 'r_induction_therapy[Doxorubicin,Intercept]',
                 p2 = 'r_induction_therapy[Cisplatin,Intercept]') %>%
   quantile(., probs = c(0.05, .5, .95)) %>%
   exp %>%
   print(digits = 5)

  compare_params(brms_object = brms_object,
                 p1 = 'r_induction_therapy[Doxorubicin,Intercept]',
                 p2 = 'r_induction_therapy[Cisplatin,Intercept]') %>%
   { mean(. > 0) } %>%
   print(digits = 5)

  compare_params(brms_object = brms_object,
                 p1 = 'r_induction_therapy[Cisplatin,Intercept]',
                 p2 = 'r_induction_therapy[No.induction,Intercept]') %>%
   quantile(., probs = c(0.05, .5, .95)) %>%
   exp %>%
   print(digits = 5)

  compare_params(brms_object = brms_object,
                 p1 = 'r_induction_therapy[Cisplatin,Intercept]',
                 p2 = 'r_induction_therapy[Cisplatin,Intercept]') %>%
   { mean(. > 0) } %>%
   print(digits = 5)

  # compare_params(brms_object = brm_tonic,
  #                p1 = 'r_induction_therapy[Doxorubicin,Intercept]',
  #                p2 = 'r_induction_therapy[Cisplatin,Intercept]')
}


#' Plot posterios distributions of offset parameters
#'
#'
posterior_histogram <- function(brms_object = brm_tonic,
                                name = '',
                                offset_val = 1,
                                group_name = 'induction_therapy',
                                param_capture = '.*\\[(.*),Intercept\\]',
                                binwidth = .1,
                                x_lab = 'Contribution to response probability',
                                cols = NULL,
                                variables = c('Doxorubicin', 'Cisplatin')) {
  ## Subselect group variables
  p_names <- names(brms_object$fit@sim$samples[[1]])
  p_names <- grep(param_capture, p_names, value = T)
  p_names <- setNames(p_names,
                      tolower(gsub(param_capture, '\\1', p_names)))
  if (all(sapply(variables, is.null)))
    variables <- names(p_names)
  if (is.null(cols)) {
    cols <- setNames(gen_color_vector(name = 'Royal1',
                                      n = length(p_names)),
                     names(p_names))
  }
  cols <- setNames(cols, gsub(' ', '\\.', tolower(names(cols))))

  p_dat <- as.data.frame(extract_params(brms_object,
                                        p_names[tolower(variables)]))
  # comp_val <- setdiff(1:length(p_names), offset_val)[1]
  # prop_low <- mean(p_dat[, offset_val] < p_dat[, comp_val])
  prop_low <- NULL
  N_samples <- max(1e3, nrow(p_dat))
  p_dat <- melt(p_dat[sample(1:nrow(p_dat))[1:N_samples], ], formula = ~ .)
  medians <- p_dat %>% group_by(variable) %>% summarise('med' = median(value))
  p_dat <- as.data.table(p_dat)
  # p_dat <- remove_outliers(p_dat, by_cols = 'variable', test_cols = 'value')

  p <- ggplot(dplyr::filter(p_dat, variable %in% tolower(variables)),
              aes(x = value, fill = variable, y=..density..)) +
    geom_histogram(alpha = .3, position = 'identity', binwidth = binwidth) +
    scale_x_continuous(trans = 'identity', name = x_lab) +
    theme(legend.position = 'right', legend.direction = 'vertical')

  p <- p + geom_vline(data = medians,
                      mapping = aes(xintercept = med, colour = variable),
                      size = 1.5, alpha = .8)
  p <- p + scale_colour_manual(values = cols, guide = F)
  p <- p + scale_fill_manual(values = cols, name = '')
  text_size = 2.5
  if (!is.null(prop_low)) {
    p <- p + annotate(geom = 'text',
                      y = interpolate_in_gg_range(plot = p, axis = 'y',
                                                  degree = .9),
                      x = interpolate_in_gg_range(plot = p, axis = 'x',
                                                  degree = .05),
                      parse = T,
                      hjust = 0,
                      size = text_size,
                      # label = parse(text=sprintf('P(%s<%s)==%.2f',
                      #                 variables[1], variables[2], prop_low)))
                      label = sprintf('italic(P)(%s<%s)==%.2f',
                                      gsub(' ', '~', variables[offset_val]),
                                      gsub(' ', '~', variables[comp_val]),
                                      prop_low))
  }

  p <- p + annotate(geom = 'text',
                    y = interpolate_in_gg_range(plot = p, axis = 'y',
                                                degree = .98),
                    x = interpolate_in_gg_range(plot = p, axis = 'x',
                                                degree = .05),
                    parse = F,
                    hjust = 0,
                    size = text_size,
                    # label = parse(text=sprintf('P(%s<%s)==%.2f',
                    #                 variables[1], variables[2], prop_low)))
                    label = name)
  return(p)
}


prep_comp_data <- function(tp1 = 'Baseline', tp2 = 'On nivo') {
  dtf <- prepare_test_gene_set_difference(gene_set = sig_gene_sets,
                                          tp1 = tp1,
                                          tp2 = tp2)
  invisible(dtf[, 'FC' := get(tp2) - get(tp1)])
  dtf <- controlled_merge(dtf, patient_labels[, .(patient, clinical_response,
                                                  s_til, ca15_3)],
                          by_cols = 'patient') %>%
    { .[complete.cases(.), ] }
  return(dtf)
}


explore_confounder <- function(conf_var = 's_til') {
  if (!require(ggExtra)) {
    devtools::install_github('daattali/ggExtra')
    library(ggExtra)
  }
  p1 <- ggplot(dtf, aes_string(x = conf_var, y = 'FC', colour = 'arm')) +
    scale_colour_manual(name = '', values = tonic_color_palettes[['arm']]) +
    geom_point() +
    scale_y_continuous(name = sprintf('log2 Fold change between %s and %s\nin the median of response-associated genesets scores', tp1, tp2)) +
    scale_x_continuous(name = var_to_label(conf_var, label_reps)) +
    ggplot2::theme(legend.position = 'left', legend.direction = 'vertical',
                   aspect.ratio = 1)
  ggMarginal(p1, groupColour = TRUE, groupFill = TRUE, type = 'histogram')
}


plot_params <- function(p_dat = hlm, mn = 'on_nivo', plot_normalized = T,
                        save_plot = T) {
  if (plot_normalized == FALSE) {
    par_n <- intersect(names(p_dat),
                       c('sigma', 'bc[1]', sprintf('Y_arm[%d]', 5:1)))
    arm_means <- grep('Y_arm', par_n, value = T)
    par_n_print <- c(setdiff(par_n, arm_means), 
                     paste0('Y_', treatment_arms)) %>%
      { gsub('bc\\[1\\]', 'Baseline', .) }
    p <- stan_plot(p_dat, pars = par_n)
    p <- p + theme_ms()
    p <- p + scale_x_continuous(name = 'Credible interval')
    p <- p + scale_y_discrete(name = '', limits = rev(par_n),
                              labels = rev(par_n_print))
  } else {
    par_n <- intersect(names(p_dat),
                       c(sprintf('Y_arm_normalized[%d]', 1:4)))
    arm_means <- grep('Y_arm_normalized', par_n, value = T)
    par_n_print <- c(setdiff(par_n, arm_means), treatment_arms[2:5]) %>%
      { gsub('bc\\[1\\]', 'Baseline', .) }

    probs <- sapply(2:5, function(idx) {
      rstan::extract(p_dat, pars = c('Y_arm[1]', sprintf('Y_arm[%d]', idx))) %>%
        as.data.table %>%
        set_colnames(c('ref', 'test')) %>%
        { .[, .(mean(test >= ref))][[1]] }
    }) %>% setNames(par_n_print)

    posterior2 <- extract(p_dat, inc_warmup = F, permuted = FALSE)
    pacman::p_load(bayesplot)
    p <- mcmc_areas(posterior2, pars = par_n)
    p <- p + theme_ms()
    p <- p + scale_y_discrete(name = '', limits = rev(par_n),
                              labels = rev(sprintf('%s (%s)',
                                               names(probs),
                                               scales::percent(probs))))
    p <- p + scale_x_continuous(name = 'Mean logFC compared to no induction')
  }

  extra_pars <- c('rc[1]', 'bc[1]') %>% intersect(names(p_dat))
  rel_heights <- { .2 * (1.5^(length(extra_pars))) } %>% 
    { c((1 - as.numeric(.)), .) }

  if (length(extra_pars) > 0) {
    p_ann <- plot_extra_params(p_dat, params = extra_pars)
    # pacman::p_load(ggpubr)
    # pacman::p_load(gridExtra)
    pacman::p_load(cowplot)
    # p_comb <- ggpubr::ggarrange(plotlist = list(p, p_ann), ncol = 1, align = 'v',
    #                             heights = c(.8, .2))
    # p_comb <- grid.arrange(p, p_ann, ncol = 1)
    p_comb <- cowplot::plot_grid(plotlist = list(p, p_ann), 
                                 align = 'v',
                                 ncol = 1, 
                                 rel_heights = rel_heights)
  } else {
    p_comb <- p
  }

  if (save_plot) {
    fn <- file.path(img_dir, sprintf('bayesian_fc_%s.pdf', mn))
    ggsave(p_comb, filename = fn,
           width = 9, height = 7, units = 'cm')
  }
  return(p_comb)
}


simulate_data <- function(stanmodel = hlm_pi, nsims = 5) {
  params <- rstan::extract(stanmodel, pars = names(stanmodel), 
                           inc_warmup = F) %>%
    as.tibble %>%
    { .[1001:nrow(.), ] } %>%
    dplyr::arrange(-lp__) %>%
    { .[1:nsims, ] } %>%
    dplyr::select(c(sprintf('Y_arm[%d]', 1:5), 'sigma')) %>%
    set_colnames(c(treatment_arms, 'sigma'))

  rnorm_vec <- function(means, sd) {
    sapply(means, function(x) rnorm(1, mean = unlist(x), sd = unlist(sd)))
  }

  all_dat <- rbind(purrr::map_dfr(1:nsims, function(iter) {
    data.table(dplyr::select(data_pi, -FC), 
               'FC' = data_pi[, rnorm_vec(means = params[iter, as.integer(arm)], 
                                          sd = params[iter, 'sigma'])],
               'type' = sprintf('sim_%d', iter))
  }), cbind(data_pi, 'type' = 'observed'))

  ggplot(all_dat, aes(x = arm, y = FC)) + 
   ggbeeswarm::geom_quasirandom(width = .1, alpha = .8) + 
   geom_boxplot(outlier.size = 0, alpha = .5, fill = 'lightblue') + 
   facet_wrap(~type) +
   rotate_x_labels(45) + 
   labs(list(x = '', y = 'FC'))
}


plot_extra_params <- function(p_dat = hlm_pi, params) {
  p <- stan_plot(p_dat, pars = params)
  p <- p + theme_ms()
  p <- p + scale_x_continuous(name = 'Credible interval')
  p <- p + scale_y_discrete(name = '', 
                            limits = params,
                            labels = c('bc[1]' = 'Baseline contribution', 
                                       'rc[1]' = 'Clinical response contribution') %>%
                                       { .[params] } %>% rev)
  return(p)
}
