plot_stability_paths <- function(stability_paths, path_pars = list(),
                                 label_pars = list(), ...) {
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))

  checkmate::assert_matrix(stability_paths, mode = "numeric",
    any.missing = FALSE, col.names = "named")

  new_pars <- list(...)
  default_pars <- list(mar = opar$mar + c(0, 0, 0, 1))
  do.call(par, modifyList(default_pars, new_pars))

  n_covariates <- ncol(stability_paths)
  n_step <- nrow(stability_paths)

  # default qualitative palette:
  col <- hcl(
    h = seq(0, 360, l = n_covariates + 1) + 15, c = 100,
    l = 65
  )[-(n_covariates + 1)]

  default_path_pars <- list(
    lty = 1, type = "b", bty = "n", col = col,
    x = 0:(n_step - 1), pch = 19, xlab = "# covariates",
    ylab = expression(Pi), ylim = c(0, 1)
  )
  path_pars <- c(
    list(y = stability_paths),
    modifyList(default_path_pars, path_pars)
  )
  do.call(matplot, path_pars)

  # arrange labels according to mean inclusion probability
  rank_covariate <- rank(colMeans(stability_paths), ties.method = "first")
  label_pos <- seq(0, 1, l = n_covariates)[rank_covariate]
  ## alternatively:  vertical position = mean inclusion probability
  ## # label_pos <- colMeans(stability_paths)
  ## ## scale to [0, 1] for better spacing:
  ## # label_pos <- (label_pos - min(label_pos))/diff(range(label_pos))

  default_label_pars <- list(
    text = colnames(stability_paths), side = 4,
    at = label_pos, col = col, las = 1
  )
  label_pars <- modifyList(default_label_pars, label_pars)

  do.call(mtext, label_pars)
}
