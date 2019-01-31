library(randomForest)
library(mlbench)
library(DEoptimR)
# library(caret)
# devtools::install_github('topepo/caret/pkg/caret')
library(caret)
library(glmnet)

prepare_ml_dat <- function(tp = 'Baseline', type = 'gene_sets') {
  type <- match.arg(type, choices = c('genes', 'gene_sets'), several.ok = F)

  if (tolower(tp) == 'all') {
    # tp <- c('Baseline', 'On-nivo', 'Post-induction')
    tp <- patient_labels[, levels(timepoint)]
  }

  ## Training data info
  t_dat_info <- patient_labels[timepoint %in% tp,
                               .(filename, arm, response)] %>%
    { .[!is.na(response)] } %>%
    { .[response != 'SD'] }

  if (type == 'genes') {
    # t_dat_info %>% .[, table(response)]
    # rownames(exp_levels) <- exp_levels[, gene_symbol]
    ## Subselect relevant columns from expression matrix
    t_dat <- exp_levels[, .SD, .SDcols = t_dat_info[, as.character(filename)]] %>%
      { log2(. + 1) } %>%
      ## Scale data
      # base::scale(.)
      t
    colnames(t_dat) <- exp_levels[, gene_symbol]
  } else if (type == 'gene_sets') {
    browser()
  }

  ## Get rid of some genes with NA values for some patients
  t_dat <- t_dat[, apply(t_dat, 2, function(x) all(!is.na(x)))]

  ## Look up responses
  setkey(t_dat_info, filename)
  responses <- t_dat_info[rownames(t_dat),
                 setNames(response %in% c('PR', 'CR'), rownames(t_dat))]
  return(list('t_dat' = t_dat, 'responses' = responses))
}


glmnet_class_test <- function(alphas = seq(0.0, 1, by = .1),
                              ncol = 4,
                              nrow = 3,
                              tp = 'Baseline',
                              nfolds = nrow(t_dat)) {
  prep <- prepare_ml_dat(tp = tp)
  t_dat <- prep[['t_dat']]
  responses <- prep[['responses']]

  print(table(responses))

  if (sum(responses) == 0) {
    stop('no responders among patients with this timepoint')
  }

  if (length(alphas) > (ncol * nrow)) {
    stop('adjust amount of plots to number of alphas')
  }

  indices <- 1:(ceiling(length(alphas) / ncol / nrow) * ncol * nrow)
  layout(matrix(indices, ncol = ncol, nrow = nrow, byrow = T))

  for (a in alphas) {
    res <- cv.glmnet(x = t_dat,
                     y = responses,
                     alpha = a,
                     family = 'binomial',
                     type.measure = 'class',
                     nfolds = nfolds)
    plot(res)
    abline(a = mean(responses), b = 0, type = 'dashed')
    title(parse(text=sprintf('alpha==%s', a)))
  }
  return(res)
}


perform_rf <- function(tp = 'Baseline') {
  prep <- prepare_ml_dat_rnaseq(tp = tp)
  rf_default <- caret::train(x = prep[[1]], y = prep[[2]],
                             method="rf",
                             metric="Kappa",
                             tuneGrid = expand.grid(.mtry=ceiling(sqrt(ncol(t_dat)))),
                             trControl = trainControl(method="repeatedcv",
                                                      number=nrow(t_dat),
                                                      repeats=10))

  return(rf_default)
}


#' Plot variable importance for RF object
#'
#'
plot_variable_importance <- function(fit) {
  vec <- importance(fit)
  p_dat <- data.frame(gini = vec,
                      gene = rownames(vec), row.names = rownames(vec))
  p_dat <- p_dat %>% arrange(-MeanDecreaseGini) %>%
    as.data.table %>%
    transform(perc_rank = frank(-MeanDecreaseGini, ties.method = 'min') - 1) %>%
    transform(perc_rank = perc_rank / max(perc_rank))

  library(ggrepel)

  ggplot(p_dat, aes(x = perc_rank, y = MeanDecreaseGini, label = gene)) +
    geom_point(alpha = .5, size = .5) +
    scale_x_continuous(name = 'Rank percentile') +
    geom_label_repel(data = p_dat[perc_rank <= .25], label.size = .1,
                     size = 1.5,
                     segment.alpha = .5)
}


#' Plot variable importance for RF object (easier)
#'
#'
plot_variable_importance <- function(fit) {
  varImpPlot(fit, type = 2,
             n.var = min(100, nrow(fit$importance, main = tp)))
}
