## extract stability paths from refits of <model> on resampled <data>.
## arguments:
##    model: a fitted model object of class "regsubsets"
##    data: the data used for the model fit
##    method: resampling method, either "subsample" (without replacement) or
##      "bootstrap" (with replacement)
##    strata: for subsampling, vector (length = nrow(data)) defining the strata
##      for stratified sampling.
##    fraction: subsampling fraction
## return: a <max. subsetsize + 1> x <covariates> matrix with relative selection
##    frequencies. first row is the null model (no covariates,i.e., all 0s)
## dependencies: {checkmate}, {leaps}
get_stability_paths <- function(model, data, reps = 100,
                                method = c("subsample", "bootstrap"),
                                strata = NULL, fraction = 0.5) {
  checkmate::assert_class(model, "regsubsets")
  checkmate::assert_data_frame(data)
  checkmate::assert_count(reps)
  method <- match.arg(method)
  checkmate::assert_vector(strata, any.missing = FALSE,
                           len = NROW(data), null.ok = TRUE)
  checkmate::assert_number(fraction, lower = 0, upper = 1)

  selected <- vector("list", reps)
  for (i in seq_len(reps)) {
    new_data <- resample(data, method = method, strata = strata,
                         fraction = fraction)
    new_model <- refit(model, new_data)
    selected[[i]] <- get_selected(new_model)
  }
  stability_paths <- make_paths(selected)
  stability_paths
}

############## resample ########################################################

resample <- function(data, method = c("subsample", "bootstrap"),
                     strata = NULL, fraction = 0.5) {
  nrows <- nrow(data)
  rows <- resample_rows(nrows, method, strata, fraction)
  data[rows, ]
}

resample_rows <- function(nrows, method, strata = NULL, fraction = 0.5) {
  switch(method,
         "bootstrap" = sample_with_replacement(nrows, strata),
         "subsample" = sample_without_replacement(nrows, strata,
                                                  fraction = fraction)
  )
}

sample_with_replacement <- function(nrows, strata = NULL) {
  if (is.null(strata)) {
    return(sample(nrows, replace = TRUE)) # --> early exit!
  }
  rows <- tapply(
    X = seq_len(nrows), INDEX = strata, FUN = sample, replace = TRUE
  )
  as.vector(rows)
}

############## refit ###########################################################

# redo subset selection <model> on <new_data>1
refit <- function(model, new_data) {
  # works by overwriting the data argument of the original model
  # and then re-doing the function call that produced the original model
  modelcall <- model$call
  modelcall$data <- new_data
  # use regsubsets-generic here instead of regsubsets.formula or other method as
  # these methods are not exported by {leaps}
  # (quote s.t. just the name of the function is handed over, not the
  # function code itself...)
  modelcall[[1]] <- quote(leaps::regsubsets)
  eval(modelcall)
}

# ------------------------------------------------------------------------------
# BOTTOM-LEVEL FUNCTIONS FOR STABILITY PATH COMPUTATION
# ------------------------------------------------------------------------------

# SUBSAMPLING ------------------------------------------------------------------

sample_without_replacement <- function(nrows, 
                                       strata = NULL, 
                                       fraction = 0.5) {
  
  # In absence of strata, just sample from rows w/o replacement
  
  if (is.null(strata)) {
    
    return(sample(
      seq_len(nrows), 
      size = floor(fraction * nrows), # recommendation Meinshausen/Buehlmann
      replace = FALSE))
    
  }
  
  # In presence of strata, do the same per stratum
  # Choose number via ceiling and floor uniformly at random, so that, in 
  # expectation, n / 2 rows are sampled
  
  row_indices <- lapply(
    unique(strata),
    function(i) {
      sample(
        seq_len(nrows)[strata == i], 
        size = ifelse(
          runif(1) > 0.5, 
          ceiling(sum(strata == i) * fraction),
          floor(sum(strata == i) * fraction)),
        replace = FALSE)})
  
  unlist(row_indices)
  
}

# FINDING SELECTED COVARIATES --------------------------------------------------

get_selected <- function(model) {
  
  # Retrieve selected variables, remove intercept and add 0 row for constant
  # model (corresponds to maximum shrinkage)
  
  selected <- summary(model)$which
  selected_effects <- selected[, colnames(selected) != "(Intercept)"]
  rbind("0" = rep(FALSE, ncol(selected_effects)), selected_effects)
  
}

# ESTIMATING SELECTION PROBABILITY ---------------------------------------------

make_paths <- function(selected) {
  
  Reduce("+", selected) / length(selected)
  
}

# PLOTTING STABILITY PATHS -----------------------------------------------------

plot_stability_paths <- function(stability_paths) {
  
  as.data.frame(stability_paths) %>% 
    rownames_to_column("n_covariates") %>%
    gather("covariate", "estimated_prob", -n_covariates) %>% 
    ggplot(aes(x = n_covariates, y = estimated_prob, group = covariate)) +
    geom_point(aes(col = covariate)) +
    geom_line(aes(col = covariate)) +
    xlab("# covariates") +
    ylab(expr(Pi))
  
}