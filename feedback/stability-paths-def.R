sample_without_replacement <- function(nrows, strata = NULL, fraction = 0.5) {
  if (is.null(strata)) {
    return(sample(1:nrows, size = ceiling(nrows * fraction), replace = FALSE))
  }
  # vectorize subsampling over:
  # <rows belonging to the different strata>,
  strata_index <- split(seq_len(nrows), strata)
  # <subsample sizes in different strata>
  # round up so that subsample has at least 1 member from each strata
  strata_sizes <- ceiling(table(strata) * fraction)
  rows <- mapply(sample,
    x = strata_index, size = strata_sizes,
    replace = FALSE
  )
  unlist(rows)
}

# extract <# subset sizes> x <# covariates> inclusion indicator matrix
# from fitted subset selection model <new_model>
get_selected <- function(new_model) {
  selected <- summary(new_model)$which
  # add (named) row for null model / maximal regularization: intercept only
  selected <- rbind("0" = rep(FALSE, ncol(selected)), selected)
  # remove intercept column & return
  selected[, colnames(selected) != "(Intercept)"]
}

# collapse list of selection indicator matrices <selected> into
# a single matrix containing frequencies of inclusion
make_paths <- function(selected) {
  # add all matrices elementwise and divide by their number to
  # get elementwise means
  Reduce(`+`, selected) / length(selected)
}

