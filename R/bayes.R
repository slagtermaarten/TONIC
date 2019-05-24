library(shinystan)
library("rstan")
# if (!require(brms)) {
#   devtools::install_github("paul-buerkner/brms")
# }
# library(brms)


prep_comp_data <- function(tp1 = 'Baseline', tp2 = 'On nivo') {
  dtf <- prepare_test_gene_set_difference(gene_set = sig_gene_sets,
                                          tp1 = tp1,
                                          tp2 = tp2)
  invisible(dtf[, 'FC' := get(tp2) - get(tp1)])
  dtf <- controlled_merge(dtf,
    patient_labels[, .(patient, clinical_response, s_til, ca15_3,
      lines_of_therapy_for_metastatic_disease, lymphnode_only_disease)],
    by_cols = 'patient') %>%
    { .[complete.cases(.), ] }
  dtf[patient == 'pat_5', lymphnode_only_disease := F]
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


plot_params <- function(
  p_dat = hlm, mn = 'post_induction', plot_normalized = T, save_plot = F,
  ncol = 1) {
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
      c(sprintf('FC_arm_normalized[%d]', 1:4)))
    arm_means <- grep('FC_arm_normalized', par_n, value = T)
    par_n_print <- c(setdiff(par_n, arm_means), treatment_arms[2:5]) %>%
      { gsub('bc\\[1\\]', 'Baseline', .) }

    probs <- sapply(2:5, function(idx) {
      rstan::extract(p_dat,
        pars = c('FC_arm[1]', sprintf('FC_arm[%d]', idx))) %>%
        as.data.table %>%
        set_colnames(c('ref', 'test')) %>%
        { .[, .(mean(test >= ref))][[1]] }
    }) %>% setNames(par_n_print)

    library(bayesplot)
    posterior2 <- rstan::extract(p_dat, inc_warmup = F, permuted = FALSE)
    p <- bayesplot::mcmc_areas(posterior2, pars = par_n)
    p <- p + theme_ms()
    p <- p + scale_y_discrete(name = '', limits = rev(par_n),
      labels = rev(sprintf('%s (%s)', names(probs), scales::percent(probs))))
    p <- p + xlab(expression(mu['arm']-mu['no induction']))
  }

  extra_pars <- c('bc[1]', 'rc[1]', 'tc[1]', 'lc[1]') %>%
    intersect(names(p_dat))
  rel_heights <- max(.15 * length(extra_pars), .2)
  ## 2019-02-02 11:36 Pipe forwarding of numerics gives unexpected behaviour, so
  ## I'm opting for the old-fashioned way of doing it instead
  rel_heights <- c(1 - rel_heights, rel_heights)

  if (length(extra_pars) > 0) {
    p_ann <- plot_extra_params(p_dat, params = extra_pars)
    # library(ggpubr)
    # library(gridExtra)
    library(cowplot)
    # p_comb <- ggpubr::ggarrange(plotlist = list(p, p_ann), ncol = 1, align = 'v',
    #                             heights = c(.8, .2))
    # p_comb <- grid.arrange(p, p_ann, ncol = 1)
    if (ncol == 1) {
      p_comb <- cowplot::plot_grid(plotlist = list(p, p_ann),
        align = 'v', ncol = ncol, rel_heights = rel_heights)
    } else if (ncol == 2) {
      p_comb <- cowplot::plot_grid(plotlist = list(p, p_ann),
        align = 'h', ncol = ncol, rel_widths = c(.5, .5))
    }
  } else {
    p_comb <- p
  }

  if (save_plot) {
    fn <- file.path(img_dir, sprintf('bayesian_fc_%s.pdf', mn))
    ggsave(p_comb, filename = fn,
           width = 9, height = 7, units = 'cm')
    cat(sprintf('Printed to %s\n', fn))
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
                            # labels = c('bc[1]' = 'Baseline expression',
                            #            'lc[1]' = 'Lymph-node only metastasis',
                            #            'tc[1]' = '# Prior treatment lines',
                            #            'rc[1]' = 'Clinical response') %>%
                            labels = c('bc[1]' = 'b',
                                       'lc[1]' = 'l',
                                       'tc[1]' = 't',
                                       'rc[1]' = 'c') %>%
                                       { .[params] } %>% rev %>% setNames(NULL))
  return(p)
}


#' R wrapper around the Stan model defined in 'rmd/stan_model_cutting_edge.stan'
#'
#'
run_sim <- function(data,
                    df_sigma_FC = 3,
                    df_sigma_arm = 3,
                    df_b = 3,
                    df_response_contrib = 3,
                    df_FC_mu = 3,
                    use_prior_only = F,
                    use_clinical_response = F,
                    use_baseline = F,
                    use_LN = F,
                    use_TL = F,
                    TL_thresh = 1,
                    iter = 5e5,
                    chains = 10,
                    stan_model = 'rmd/stan_model_t_version.stan',
                    adapt_delta = .9) {
  stan_dat <- list(
    prior_only = as.integer(use_prior_only),
    use_baseline = as.integer(use_baseline),
    use_clinical_response = as.integer(use_clinical_response),
    N_arms = data[, uniqueN(arm)],
    N_obs = nrow(data),
    arm = as.integer(data[, arm]),
    clinical_response = as.integer(data$clinical_response == 'R'),
    TL = data$lines_of_therapy_for_metastatic_disease,
    LN = as.integer(data$lymphnode_only_disease),
    'df_FC_mu' = df_FC_mu,
    'df_sigma_arm' = df_sigma_arm,
    'df_sigma_FC' = df_sigma_FC,
    'df_b' = df_b,
    'df_response_contrib' = df_response_contrib,
    baseline = data[, Baseline],
    FC = data[, FC]
  )

  if (!is.null(TL_thresh)) {
    stan_dat$TL <- stan_dat$TL < TL_thresh
  }

  if (!use_baseline) {
    stan_dat$baseline <- rep(0, nrow(data))
  }

  if (!use_clinical_response) {
    stan_dat$clinical_response <- rep(0, nrow(data))
  }

  if (!use_TL) {
    stan_dat$TL <- rep(0, nrow(data))
  }

  if (!use_LN) {
    stan_dat$LN <- rep(0, nrow(data))
  }

  hlm <- stan(model_name = 'Hierarchical Model',
    file = stan_model,
    data = stan_dat,
    iter = iter,
    warmup = ceiling(iter / 4),
    control = list(adapt_delta = adapt_delta),
    chains = chains, verbose = FALSE)
  return(hlm)
}


gen_all_model_fits <- function(
  dtf = prep_comp_data(tp1 = 'Baseline', tp2 = 'Post-induction'),
  stan_model = 'rmd/stan_model_t_version.stan',
  iter = 1e5) {
  sims <- list(
    hlm = run_sim(data = dtf,
      stan_model = stan_model,
      iter = iter,
      use_baseline = F, use_clinical_response = F),
    hlm_b = run_sim(data = dtf,
      stan_model = stan_model,
      iter = iter,
      use_baseline = T, use_clinical_response = F),
    hlm_r = run_sim(data = dtf,
      stan_model = stan_model,
      iter = iter,
      use_baseline = F, use_clinical_response = T),
    hlm_b_r = run_sim(data = dtf,
      stan_model = stan_model,
      iter = iter,
      use_baseline = T, use_clinical_response = T),
    hlm_b_r_d = run_sim(data = dtf,
      stan_model = stan_model,
      iter = iter,
      use_baseline = T, use_clinical_response = T,
      use_TL = T, use_LN = F),
    hlm_b_r_l = run_sim(data = dtf,
      stan_model = stan_model,
      iter = iter,
      use_baseline = T, use_clinical_response = T,
      use_TL = F, use_LN = T),
    hlm_b_r_d_l = run_sim(data = dtf,
      stan_model = stan_model,
      iter = iter,
      use_baseline = T, use_clinical_response = T,
      use_TL = T, use_LN = T)
  )
  return(sims)
}


study_covar_values <- function(object_list, name = '') {
  purrr::map(object_list, function(obj) {
    i_pars <- c('FC_arm_normalized', intersect(fitted_coefs(obj), 
      c('bc[1]', 'rc[1]', 'tc[1]', 'lc[1]')))
    if (length(i_pars) == 0) return(NULL)
    rstan::extract(obj, pars = i_pars) %>%
      as.data.table %>%
      { .[, lapply(.SD, function(x) 
        sprintf('%.2f [%.2f, %.2f]', 
          quantile(x, .5), quantile(x, .1), quantile(x, .9))
        )] } %>%
      c(list('modelname' = name))
  })
}


# rbindlist(purrr::imap(pi_BM_sims, function(obj, n) {
#   i_pars <- intersect(fitted_coefs(obj), 
#     c('FC_arm_normalized[1]', 'FC_arm_normalized[2]', 
#       'FC_arm_normalized[3]', 'FC_arm_normalized[4]'))
#   if (length(i_pars) == 0) return(NULL)
#   rstan::extract(obj, pars = i_pars) %>%
#     as.data.table %>%
#     { .[, lapply(.SD, mean)] } %>%
#     cbind(., 'name' = n)
# }))
