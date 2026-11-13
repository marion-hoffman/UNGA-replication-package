# Rubin's rules pooling for a list of model results
# Each result should be a list with fields:
# - coef: named numeric vector
# - vcov: variance-covariance matrix with matching names
pool_rubin <- function(results_list) {
  m <- length(results_list)
  if (m == 0) stop("No results to pool")
  # filter out NULLs
  results_list <- Filter(Negate(is.null), results_list)
  m <- length(results_list)
  if (m == 0) stop("No non-null results to pool")

  # Find common parameter names
  param_sets <- lapply(results_list, function(r) names(r$coef))
  common_params <- Reduce(intersect, param_sets)
  if (length(common_params) == 0) stop("No common parameters to pool")

  # Build matrix of estimates (m x p)
  betas <- do.call(rbind, lapply(results_list, function(r) as.numeric(r$coef[common_params])))
  colnames(betas) <- common_params

  # Pooled mean
  beta_bar <- colMeans(betas, na.rm = TRUE)

  # Within-imputation variance W (average vcov)
  W <- Reduce(`+`, lapply(results_list, function(r) {
    vc <- r$vcov
    if (is.null(dim(vc)) || length(vc) == 0) {
      # No variance provided; approximate with zeros matrix
      matrix(0, nrow = length(common_params), ncol = length(common_params), dimnames = list(common_params, common_params))
    } else {
      vc[common_params, common_params, drop = FALSE]
    }
  })) / m

  # Between-imputation variance B (sample covariance of coef)
  centered <- sweep(betas, 2, beta_bar, FUN = "-")
  B <- t(centered) %*% centered / (m - 1)

  # Total variance T
  Tvar <- W + (1 + 1/m) * B

  list(coef = beta_bar, vcov = Tvar, m = m)
}

# Summarize pooled result from pool_rubin
# pooled: list with fields coef (named numeric), vcov (matrix), m (integer)
# conf.level: confidence level for Wald intervals
summary_pooled <- function(pooled, conf.level = 0.95) {
  if (is.null(pooled$coef) || is.null(pooled$vcov)) stop("Invalid pooled object: missing coef/vcov")
  beta <- pooled$coef
  V <- pooled$vcov
  se <- sqrt(diag(V))
  z <- beta / se
  p <- 2 * stats::pnorm(-abs(z))
  alpha <- 1 - conf.level
  zcrit <- stats::qnorm(1 - alpha/2)
  lower <- beta - zcrit * se
  upper <- beta + zcrit * se
  out <- data.frame(
    term = names(beta),
    estimate = as.numeric(beta),
    std.error = as.numeric(se),
    statistic = as.numeric(z),
    p.value = as.numeric(p),
    conf.low = as.numeric(lower),
    conf.high = as.numeric(upper),
    row.names = NULL,
    check.names = FALSE
  )
  attr(out, "m") <- pooled$m
  out
}

# Convenience: summarize one model from pooled_models list
# pooled_models: list as saved in pooled_models.RDS (named by model)
# model_name: one of names(pooled_models)
# conf.level: confidence level
summary_pooled_models <- function(pooled_models, model_name, conf.level = 0.95) {
  if (!(model_name %in% names(pooled_models))) stop("model_name not found in pooled_models")
  summary_pooled(pooled_models[[model_name]], conf.level = conf.level)
}
