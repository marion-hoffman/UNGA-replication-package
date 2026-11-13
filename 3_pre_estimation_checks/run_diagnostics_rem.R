# Diagnostics for the REM (conditional logit / clogit) fitted on
# eventlist_for_rem.csv. Each risk set = observed pair + sampled non-events
# identified by ORDER. Run AFTER all covariates have been added (full
# dataset, one random order). Uses the first imputed sequence (rep_01)
# of the parallel pipeline.

set.seed(123)
source(file.path(Find(function(d) file.exists(file.path(d, "paths.R")),
                      c(".", "..", "../..")), "paths.R"), chdir = TRUE)

library(survival)

all_events <- read.csv(REP01_EVENTLIST, stringsAsFactors = FALSE)

predictor_cols <- setdiff(colnames(all_events), c("IS_OBSERVED", "SOURCE", "TARGET", "TIME", "ORDER",
                                                  "num_possible", "num_sampled", "weight", "offset"))

model_data <- all_events[, c("TIME", "IS_OBSERVED", "ORDER", predictor_cols,
                             "offset", "num_possible", "num_sampled")]

# Full model (random-order effects) -- the reported specification.
# Dyad attributes enter as proportion/average summaries (trade as the
# count above the yearly mean); the 'all' unanimity dummies are excluded
# (see section 3).
fml_full <- IS_OBSERVED ~ ego_log_gdp + ego_conflict + ego_hr +
  ego_staff + ego_group + ego_sponsor1year +
  alter_topic_hr + alter_topic_coop + alter_topic_ispa + alter_nrefs +
  ego_hr:alter_topic_hr +
  alter_recur + alter_cons + alter_sponsor1year +
  dyad_group_prop + dyad_regime_prop +
  dyad_diplo_av + dyad_trade_gtmean + dyad_recur + dyad_4cycle_sponsor1year +
  offset(offset) + strata(ORDER)

# Keep rows complete on the full-model variables so that fitted values,
# residuals and the data stay row-aligned in every section below.
keep <- stats::complete.cases(model_data[, intersect(all.vars(fml_full), names(model_data))])
cat("Dropping", sum(!keep), "rows with missing full-model covariates\n")
model_data <- model_data[keep, ]
all_events <- all_events[keep, ]

cat("Fitting full model...\n")
fit_full <- clogit(fml_full, data = model_data)
print(summary(fit_full))

# All plots go to a single PDF for inspection
pdf(file.path(PRE_CHECKS_OUT, "diagnostics_rem.pdf"), width = 8, height = 6)

## ===== 1. Stability across number of sampled non-events =====
# Sub-sample k of the stored non-events per risk set (per ORDER) and recompute
# the offset as log(num_possible / k_eff), k_eff = min(k, num_sampled).
cat("\n===== 1. Stability across number of non-events =====\n")
obs_idx <- which(model_data$IS_OBSERVED == 1)
non_idx <- which(model_data$IS_OBSERVED == 0)
non_by_order <- split(non_idx, model_data$ORDER[non_idx])

# Pre-shuffle non-events within each risk set so that smaller k are NESTED
# subsets of the k=10 sample -- shared Monte-Carlo noise makes convergence
# across k visible instead of confounded by independent resampling.
shuffled <- lapply(non_by_order, sample)
k_grid <- c(1, 2, 4, 6, 8, 10)
coefs_by_k <- list()
for (k in k_grid) {
  keep <- unlist(lapply(shuffled, function(ix) ix[seq_len(min(k, length(ix)))]))
  sub <- model_data[c(obs_idx, keep), ]
  k_eff <- pmin(k, sub$num_sampled)
  sub$offset <- ifelse(sub$IS_OBSERVED == 0, log(sub$num_possible / k_eff), 0)
  fit_k <- clogit(fml_full, data = sub)
  coefs_by_k[[as.character(k)]] <- coef(fit_k)
  cat(sprintf("  fitted k=%d (n=%d)\n", k, nrow(sub)))
}
coef_tab <- do.call(cbind, coefs_by_k)
colnames(coef_tab) <- paste0("k=", k_grid)
cat("\nCoefficients by number of non-events:\n")
print(round(coef_tab, 4))

# Plot each coefficient as a function of k; dashed line = k=10 estimate.
# Interpretation: flat line approaching the reference = estimates have
# converged in the number of sampled controls.
npar <- nrow(coef_tab)
nc <- ceiling(sqrt(npar)); nr <- ceiling(npar / nc)
par(mfrow = c(nr, nc), mar = c(4, 4, 2, 1))
for (i in seq_len(npar)) {
  plot(k_grid, coef_tab[i, ], type = "b", pch = 19,
       xlab = "k sampled non-events", ylab = "coefficient",
       main = rownames(coef_tab)[i], cex.main = .8)
  abline(h = coef_tab[i, ncol(coef_tab)], lty = 2, col = "grey40")
}
par(mfrow = c(1, 1))

## ===== 2. Trade parametrization: count-above-mean vs raw mean =====
cat("\n===== 2. Trade parametrization: dyad_trade_gtmean vs dyad_trade_av =====\n")
# Same full spec, but trade coded as the raw mean importance over all
# previous sponsors (zeros included) instead of the count above the
# yearly matrix mean. Non-nested comparison on identical data -> AIC /
# pseudo-R2 / concordance decide which parametrization to keep. This
# runs BEFORE the non-linearity checks so the spline diagnostics apply
# to the winner.
fml_trade_raw <- update(fml_full, . ~ . - dyad_trade_gtmean + dyad_trade_av)
fit_trade_raw <- clogit(fml_trade_raw, data = model_data)

r2 <- function(f) 1 - f$loglik[2] / f$loglik[1]
cat(sprintf("AIC:                main %.2f | raw %.2f (Delta %+ .2f, >0 favours main)\n",
            AIC(fit_full), AIC(fit_trade_raw),
            AIC(fit_trade_raw) - AIC(fit_full)))
cat(sprintf("McFadden pseudo-R2: main %.4f | raw %.4f\n",
            r2(fit_full), r2(fit_trade_raw)))
cat(sprintf("Concordance:        main %.4f | raw %.4f\n",
            summary(fit_full)$concordance[1],
            summary(fit_trade_raw)$concordance[1]))

# Do the two parametrizations rank the candidates the same way?
lp_main <- predict(fit_full, type = "lp")       # includes the offset
lp_raw  <- predict(fit_trade_raw, type = "lp")
cat("Correlation of linear predictors:", round(cor(lp_main, lp_raw), 4), "\n")
pick_main <- ave(lp_main, model_data$ORDER, FUN = function(x) x == max(x)) > 0
pick_raw  <- ave(lp_raw,  model_data$ORDER, FUN = function(x) x == max(x)) > 0
cat(sprintf("Risk sets where argmax agrees:  %.1f%%\n",
            100 * mean(pick_main == pick_raw)))
cat(sprintf("Observed sponsor is argmax:     main %.1f%% | raw %.1f%%\n",
            100 * mean(pick_main[model_data$IS_OBSERVED == 1]),
            100 * mean(pick_raw[model_data$IS_OBSERVED == 1])))

# Shared-coefficient stability: if the rest of the model barely moves,
# the trade parametrization doesn't affect any substantive conclusions.
shared <- setdiff(names(coef(fit_full)), "dyad_trade_gtmean")
b1 <- coef(fit_full)[shared]; b2 <- coef(fit_trade_raw)[shared]
cmp <- data.frame(main = b1, raw_mean = b2,
                  rel_change = (b2 - b1) / pmax(abs(b1), 1e-9))
cat("Shared coefficients (sorted by |relative change|):\n")
print(round(cmp[order(-abs(cmp$rel_change)), ], 3))
cat(sprintf("Trade terms: dyad_trade_gtmean = %+.3f | dyad_trade_av = %+.3f\n",
            coef(fit_full)["dyad_trade_gtmean"],
            coef(fit_trade_raw)["dyad_trade_av"]))

## ===== 3. Non-linearity of continuous effects =====
# NB: the 'all' unanimity dummies (dyad_group_all, dyad_regime_all) are
# permanently dropped: mechanically redundant with 'prop', no equivalents
# for diplo/trade, not of theoretical interest, and excluding them made the
# non-linearity fits below easier.
# pspline(df = 4) splits each term into a linear + non-linear component;
# the printed table gives a Chi-square test for the non-linear part.
cat("\n===== 3. Non-linearity (pspline on continuous predictors) =====\n")
fml_nl <- IS_OBSERVED ~ ego_conflict + ego_staff + ego_group +
  alter_topic_hr + alter_topic_coop + alter_topic_ispa +
  ego_hr:alter_topic_hr + alter_recur + alter_cons +
  dyad_recur +
  pspline(ego_log_gdp, df = 4) + pspline(ego_hr, df = 4) +
  pspline(ego_sponsor1year, df = 4) + pspline(alter_nrefs, df = 4) +
  pspline(alter_sponsor1year, df = 4) + pspline(dyad_group_prop, df = 4) +
  pspline(dyad_regime_prop, df = 4) +
  pspline(dyad_trade_gtmean, df = 4) +
  pspline(dyad_diplo_av, df = 4) +
  pspline(dyad_4cycle_sponsor1year, df = 4) +
  offset(offset) + strata(ORDER)
fit_nl <- clogit(fml_nl, data = model_data)
print(fit_nl)   # check 'nonlin' Chisq p-values per term

# Fitted spline curves vs the LINEAR assumption of the reported model.
# Red dashed = purely linear contribution (both centered). Divergence over
# the data mass = the non-linearities a parametric fix must capture.
term_mat <- predict(fit_nl, type = "terms")
ps_cols <- grep("^pspline", colnames(term_mat))
ps_vars <- sub("pspline\\(([^,]+),.*", "\\1", colnames(term_mat)[ps_cols])
rows_used <- rownames(term_mat)          # rows actually used in the fit
idx_s <- sample(length(rows_used), min(20000, length(rows_used)))
b_lin <- coef(fit_full)
nc2 <- ceiling(sqrt(length(ps_cols))); nr2 <- ceiling(length(ps_cols) / nc2)
par(mfrow = c(nr2, nc2), mar = c(4, 4, 2, 1))
for (j in seq_along(ps_cols)) {
  v <- ps_vars[j]
  x <- model_data[rows_used[idx_s], v]
  y <- term_mat[idx_s, ps_cols[j]]; y <- y - mean(y)
  o <- order(x)
  plot(x[o], y[o], pch = ".", col = "grey50",
       xlab = v, ylab = "partial contribution (centered)",
       main = paste("pspline:", v), cex.main = .8)
  if (v %in% names(b_lin)) {
    xv <- model_data[[v]]
    xx <- quantile(xv, c(.02, .98), na.rm = TRUE)
    lines(xx, b_lin[v] * (xx - mean(xv, na.rm = TRUE)),
          col = "red", lty = 2, lwd = 2)
  }
}
par(mfrow = c(1, 1))

## ----- 3b. Candidate parametric specification -----------------------------
# Fixes suggested by the spline plots: cubic terms for ego_log_gdp
# (up-down-up) and alter_sponsor1year (mild S-shape); quadratic terms for
# ego_hr (accelerating), dyad_group_prop, dyad_regime_prop (inverted U),
# dyad_diplo_av and dyad_trade_gtmean; squared HR-score interaction
# ego_hr_sq:alter_topic_hr. NB: alter_nrefs is log1p-transformed at source;
# trade enters as dyad_trade_gtmean (count of previous sponsors above
# the yearly matrix-mean trade importance) -- the raw-mean dyad_trade_av
# parametrization is only used in the trade_raw robustness model.
cat("\n===== 3b. Alternative spec: do the parameters change? =====\n")
model_data$ego_log_gdp_sq <- model_data$ego_log_gdp^2
model_data$ego_log_gdp_cb <- model_data$ego_log_gdp^3
model_data$ego_hr_sq <- model_data$ego_hr^2
model_data$alter_sponsor1year_sq <- model_data$alter_sponsor1year^2
model_data$alter_sponsor1year_cb <- model_data$alter_sponsor1year^3
model_data$dyad_group_prop_sq <- model_data$dyad_group_prop^2
model_data$dyad_regime_prop_sq <- model_data$dyad_regime_prop^2
model_data$dyad_diplo_av_sq <- model_data$dyad_diplo_av^2
model_data$dyad_trade_gtmean_sq <- model_data$dyad_trade_gtmean^2

fml_alt <- IS_OBSERVED ~ ego_log_gdp + ego_log_gdp_sq + ego_log_gdp_cb +
  ego_conflict + ego_hr + ego_hr_sq +
  ego_staff + ego_group + ego_sponsor1year +
  alter_topic_hr + alter_topic_coop + alter_topic_ispa + alter_nrefs +
  ego_hr:alter_topic_hr + ego_hr_sq:alter_topic_hr +
  alter_recur + alter_cons +
  alter_sponsor1year + alter_sponsor1year_sq + alter_sponsor1year_cb +
  dyad_group_prop + dyad_group_prop_sq + 
  dyad_regime_prop + dyad_regime_prop_sq +
  dyad_diplo_av + dyad_diplo_av_sq + 
  dyad_trade_gtmean + dyad_trade_gtmean_sq +
  dyad_recur +
  dyad_4cycle_sponsor1year +
  offset(offset) + strata(ORDER)
fit_alt <- clogit(fml_alt, data = model_data)
print(summary(fit_alt))
cat("LR test, reported spec vs alternative (approximate: specs are not\n")
cat("strictly nested -- linear terms vs their polynomial extensions):\n")
print(anova(fit_full, fit_alt))
cat("AIC reported spec:   ", AIC(fit_full), "\n")
cat("AIC alternative spec:", AIC(fit_alt), "\n")

# Do the shared coefficients move? (terms with identical meaning in both
# specs -- excludes reshaped terms)
common <- intersect(names(coef(fit_full)), names(coef(fit_alt)))
cmp <- data.frame(b_reported = coef(fit_full)[common],
                  b_alternative = coef(fit_alt)[common],
                  pct_change = 100 * (coef(fit_alt)[common] - coef(fit_full)[common]) /
                    abs(coef(fit_full)[common]))
print(round(cmp[order(-abs(cmp$pct_change)), ], 3))

# For the reshaped terms the coefficients are not directly comparable;
# compare the IQR-contrast contribution f(Q75) - f(Q25) on the linear
# predictor instead (same quantity the multiplier exp() is computed on).
b_alt <- coef(fit_alt)
# Contribution of raw-scale variable v under the alternative spec, picking
# up whichever of v, v^2, v^3, log1p(v) terms exist in the fit.
contrib_alt <- function(v, x) {
  out <- 0
  sq <- paste0(v, "_sq"); cb <- paste0(v, "_cb"); lg <- paste0("log_", v)
  if (v  %in% names(b_alt)) out <- out + b_alt[v]  * x
  if (sq %in% names(b_alt)) out <- out + b_alt[sq] * x^2
  if (cb %in% names(b_alt)) out <- out + b_alt[cb] * x^3
  if (lg %in% names(b_alt)) out <- out + b_alt[lg] * log1p(x)
  out
}
reshaped <- c("ego_log_gdp", "ego_hr", "alter_nrefs", "alter_sponsor1year",
              "dyad_group_prop", "dyad_regime_prop", "dyad_diplo_av",
              "dyad_trade_gtmean")
cat("\nIQR-contrast contribution f(Q75)-f(Q25), reported vs alternative:\n")
for (v in reshaped) {
  q <- quantile(model_data[[v]], c(.25, .75), na.rm = TRUE)
  c_old <- unname(b_lin[v] * diff(q))
  c_new <- unname(contrib_alt(v, q[2]) - contrib_alt(v, q[1]))
  cat(sprintf("  %-18s reported: %+.3f   alternative: %+.3f   (ratio %.2f)\n",
              v, c_old, c_new, c_new / c_old))
}

## ----- 3c. Is the non-linearity actually solved? ---------------------------
# Refit splines ON THE MODELING SCALE of the alternative spec: if the
# parametric form captured the shape, the residual spline should be flat
# (large nonlin p-values) and the grey cloud should hug the fitted
# parametric contribution (red line).
cat("\n===== 3c. Residual curvature under the alternative spec =====\n")
fml_nl_alt <- IS_OBSERVED ~ ego_conflict + ego_staff + ego_group +
  alter_topic_hr + alter_topic_coop + alter_topic_ispa +
  ego_hr:alter_topic_hr + ego_hr_sq:alter_topic_hr + alter_recur + alter_cons +
  dyad_recur +
  pspline(ego_log_gdp, df = 4) + pspline(ego_hr, df = 4) +
  pspline(ego_sponsor1year, df = 4) + pspline(alter_nrefs, df = 4) +
  pspline(alter_sponsor1year, df = 4) + pspline(dyad_group_prop, df = 4) +
  pspline(dyad_regime_prop, df = 4) +
  pspline(dyad_trade_gtmean, df = 4) +
  pspline(dyad_diplo_av, df = 4) +
  pspline(dyad_4cycle_sponsor1year, df = 4) +
  offset(offset) + strata(ORDER)
fit_nl_alt <- clogit(fml_nl_alt, data = model_data)
print(fit_nl_alt)   # nonlin p-values should be much larger than before

term_mat <- predict(fit_nl_alt, type = "terms")
ps_cols <- grep("^pspline", colnames(term_mat))
ps_vars <- sub("pspline\\(([^,]+),.*", "\\1", colnames(term_mat)[ps_cols])
rows_used <- rownames(term_mat)
idx_s <- sample(length(rows_used), min(20000, length(rows_used)))
# Overlay = contribution implied by the alternative parametric spec;
# contrib_alt auto-detects linear/quadratic/cubic/log terms per variable.
has_alt_term <- function(v)
  any(c(v, paste0(v, "_sq"), paste0(v, "_cb"), paste0("log_", v)) %in% names(b_alt))
nc2 <- ceiling(sqrt(length(ps_cols))); nr2 <- ceiling(length(ps_cols) / nc2)
par(mfrow = c(nr2, nc2), mar = c(4, 4, 2, 1))
for (j in seq_along(ps_cols)) {
  v <- ps_vars[j]
  x <- model_data[rows_used[idx_s], v]
  y <- term_mat[idx_s, ps_cols[j]]; y <- y - mean(y)
  o <- order(x)
  plot(x[o], y[o], pch = ".", col = "grey50",
       xlab = v, ylab = "partial contribution (centered)",
       main = paste("alt spec, spline:", v), cex.main = .8)
  if (has_alt_term(v)) {
    xv <- model_data[[v]]
    xx <- seq(quantile(xv, .02, na.rm = TRUE), quantile(xv, .98, na.rm = TRUE),
              length.out = 200)
    # contrib_alt is defined on the raw scale; for log-splined vars the
    # modeling-scale column is the raw name prefixed with log_
    rv <- sub("^log_", "", v)
    xx_raw <- if (v != rv) expm1(xx) else xx
    yy <- contrib_alt(rv, xx_raw); yy <- yy - mean(contrib_alt(rv, model_data[[rv]]))
    lines(xx, yy, col = "red", lty = 2, lwd = 2)
  }
}
par(mfrow = c(1, 1))

## ===== 4. Choice-model fit =====
# In the conditional-multinomial-logit reading, each risk set is a discrete
# choice among the observed pair + sampled controls.
cat("\n===== 4. Choice-model fit =====\n")
# McFadden pseudo-R^2: 1 - LL_full/LL_null (fit$loglik = c(null, full))
cat("McFadden pseudo-R2:", round(1 - fit_full$loglik[2] / fit_full$loglik[1], 4), "\n")

# Hit rate / rank of the OBSERVED alternative within each risk set (ORDER).
# lp scores are used for ranking; the offset is constant within a risk set
# (same num_possible/num_sampled) so it does not affect ranks.
lp <- predict(fit_full, type = "lp")
md_fit <- model_data[names(lp), ]   # rows used by the fit (NAs dropped)
rk <- tapply(seq_along(lp), md_fit$ORDER, function(ix)
  rank(-lp[ix], ties.method = "average")[md_fit$IS_OBSERVED[ix] == 1])
obs_rank <- unlist(rk)
K_mean <- mean(model_data$num_sampled[!is.na(model_data$num_sampled)]) + 1
cat("Hit rate (observed pair ranked 1st):", round(mean(obs_rank == 1), 4),
    "   (chance =", round(1 / K_mean, 4), ")\n")
cat("Mean rank of observed pair:", round(mean(obs_rank), 2),
    "   (chance =", round((K_mean + 1) / 2, 2), ")\n")
cat("Share of observed pairs in top 3:", round(mean(obs_rank <= 3), 4),
    "   (chance =", round(3 / K_mean, 4), ")\n")

## ===== 5. Influential observations (dfbetas) =====
cat("\n===== 5. Influence diagnostics (dfbetas) =====\n")
# clogit uses the exact partial likelihood, for which dfbetas are not
# implemented. Refit the SAME model via coxph(Surv(1, status) ~ ... +
# strata(ORDER), ties = "efron"): with exactly one event per risk set the
# tie method is irrelevant, so coefficients are identical, and dfbetas are
# available for the Efron method.
model_data$one <- 1
fml_efron <- update(fml_full, Surv(one, IS_OBSERVED) ~ .)
fit_efron <- coxph(fml_efron, data = model_data, ties = "efron")
cat("Max |coef difference| vs clogit fit:",
    max(abs(coef(fit_efron) - coef(fit_full))), "\n")
dfb <- residuals(fit_efron, type = "dfbetas")  # scaled: change in coef / SE
max_dfb <- apply(abs(dfb), 2, max)
cat("Max |dfbetas| per coefficient:\n")
print(sort(round(max_dfb, 4), decreasing = TRUE))

# most affected coefficient = 20 -> dyad_4cycle_sponsor1year (max |dfbeta| = 0.1997)
# dfbeta sign: + means the row pushes the estimate up, - pulls it down
top_coef <- which.max(max_dfb)
top_rows <- order(abs(dfb[, top_coef]), decreasing = TRUE)[1:10]
# rownames of residuals = rownames of model_data = rownames of all_events
cat("Most influential observations for coef", top_coef, ":\n")
print(data.frame(
  all_events[top_rows, c("SOURCE", "TARGET", "TIME", "ORDER", "IS_OBSERVED")],
  dfbeta = round(unname(dfb[top_rows, top_coef]), 3),
  x      = round(all_events$dyad_4cycle_sponsor1year[top_rows], 3))) # careful, hard coded
summary(all_events$dyad_4cycle_sponsor1year) # not particularly very high values

# second most affected coefficient = 7 -> alter_topic_hr (max |dfbeta| = 0.1811)
max_dbf2 <- max_dfb
max_dbf2[top_coef] <- 0
top_coef2 <- which.max(max_dbf2)
top_rows2 <- order(abs(dfb[, top_coef2]), decreasing = TRUE)[1:10]
cat("Most influential observations for coef", top_coef2, ":\n")
print(data.frame(
  all_events[top_rows2, c("SOURCE", "TARGET", "TIME", "ORDER", "IS_OBSERVED")],
  dfbeta = round(unname(dfb[top_rows2, top_coef2]), 3),
  x      = round(all_events$alter_topic_hr[top_rows2], 3)))  # careful, hard coded
# here influential values are all for hr topic drafts

## ===== 6. Multicollinearity (VIF) =====
cat("\n===== 6. VIF on the full-model design matrix =====\n")
# NB: model.matrix on the full formula would expand strata(ORDER) into one
# dummy column per risk set (~50+ Gb). Build the design matrix from the
# predictor term labels only (interactions preserved, offset/strata dropped).
pred_labels <- attr(terms(fml_full), "term.labels")
pred_labels <- pred_labels[!grepl("^(offset|strata)", pred_labels)]
X <- model.matrix(reformulate(pred_labels, intercept = FALSE), data = model_data)
vifs <- diag(solve(cor(X)))
names(vifs) <- colnames(X)
print(sort(round(vifs, 2), decreasing = TRUE)) # all below 5, probably OK

dev.off()
cat("\nDiagnostics done. Plots written to", file.path(PRE_CHECKS_OUT, "diagnostics_rem.pdf"), "\n")
