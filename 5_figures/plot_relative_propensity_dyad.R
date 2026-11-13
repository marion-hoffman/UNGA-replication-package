# plot_relative_propensity_dyad.R
#
# Relative propensity curves R(x) = exp{ eta(x) - eta(x0) } for the four
# H3 dyadic statistics (x0 = 0) -- predicted propensity of a dyad with
# statistic x vs an otherwise identical dyad with no ties/overlap (all
# four stats are proportions in [0, 1] with substantial mass at 0).
# Sponsoring only: the initiating model has no dyadic terms.
#
# Three curves per panel:
#   No controls             -- univariate h3 spec (single linear term)
#   Full model              -- model_all_effects (single linear term)
#   Full model (non-linear) -- model_nl_effects (linear + squared term)
#
# All four stats use the same variable in every tier (the non-linear
# spec now carries dyad_trade_gtmean + its square, matching the main
# spec). rescale stays available if codings ever diverge again.
#
# Run:
#   Rscript plot_relative_propensity_dyad.R

suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
})

source(file.path(Find(function(d) file.exists(file.path(d, "paths.R")),
                      c(".", "..", "../..")), "paths.R"), chdir = TRUE)
results_dir <- SPONSORING_OUT
pooled_dir <- file.path(results_dir, "pooled")

rep1 <- file.path(results_dir, "rep_01", "eventlist_for_rem.csv")
stopifnot(file.exists(rep1))
spon_ev <- utils::read.csv(rep1)

# coef/vcov extractor for pooled lists and coxph fits
get_cv <- function(obj) {
  if (inherits(obj, c("coxph", "clogit")))
    return(list(cf = stats::coef(obj), V = stats::vcov(obj)))
  list(cf = obj$coef, V = obj$vcov)   # pooled list(coef, vcov, m)
}

# R(x) curve for a linear (+optional uncentred quadratic) term.
# xgrid/x0 are on the variable's own scale; x_plot rescales for display.
curve_df <- function(cf, V, xgrid, x0, lin, sq, tier, stat, sx = identity) {
  nm <- names(cf)
  miss <- c(lin, sq)[!c(lin, sq) %in% nm]
  if (length(miss))
    warning("curve_df: term(s) ", paste(miss, collapse = ", "),
            " absent from coef names (", tier, " / ", stat,
            ") -- curve is flat; pooled file is likely stale")
  out <- lapply(xgrid, function(x) {
    w <- setNames(rep(0, length(nm)), nm)
    if (lin %in% nm) w[lin] <- x - x0
    if (!is.null(sq) && sq %in% nm) w[sq] <- x^2 - x0^2
    d <- sum(w * cf)
    s <- sqrt(drop(t(w) %*% V %*% w))
    data.frame(x = sx(x), R = exp(d),
               lo = exp(d - 1.96 * s), hi = exp(d + 1.96 * s))
  })
  cbind(do.call(rbind, out), tier = tier, stat = stat)
}

# Per-stat configuration: linear term in univ/full specs, the non-linear
# spec's linear term (nl_var) and its square, and the pooled univ file.
stat_specs <- list(
  group  = list(lab = "Shared groups", var = "dyad_group_prop",
                sq = "dyad_group_prop_sq",
                univ = "pooled_univariate_h3_group.RDS"),
  regime = list(lab = "Same regime", var = "dyad_regime_prop",
                sq = "dyad_regime_prop_sq",
                univ = "pooled_univariate_h3_regime.RDS"),
  diplo  = list(lab = "Diplomatic ties", var = "dyad_diplo_av",
                sq = "dyad_diplo_av_sq",
                univ = "pooled_univariate_h3_diplo.RDS"),
  trade  = list(lab = "Trade ties", var = "dyad_trade_gtmean",
                nl_var = "dyad_trade_gtmean", sq = "dyad_trade_gtmean_sq",
                univ = "pooled_univariate_h3_trade.RDS")
)

minmax01 <- function(v) {
  r <- range(v, na.rm = TRUE)
  function(x) (x - r[1]) / diff(r)   # rescale to observed [0,1]
}

tier_files <- function(sp) list(
  "No controls" = list(file = sp$univ,
                       var = sp$var, sq = NULL),
  "Full model" = list(file = "pooled_model_all_effects.RDS",
                      var = sp$var, sq = NULL),
  "Full model (non-linear)" = list(file = "pooled_model_nl_effects.RDS",
                                   var = if (is.null(sp$nl_var)) sp$var
                                         else sp$nl_var,
                                   sq = sp$sq)
)

df <- NULL; rug_df <- NULL
for (st in names(stat_specs)) {
  sp  <- stat_specs[[st]]
  tfs <- tier_files(sp)
  # Rug: the linear spec's variable (normalized like the curves below).
  lin_v <- spon_ev[[sp$var]]
  rescale <- any(vapply(tfs, function(t) t$var != sp$var, logical(1)))
  sx_lin <- if (rescale) minmax01(lin_v) else identity
  rug_df <- rbind(rug_df,
    data.frame(x = sx_lin(lin_v[spon_ev$IS_OBSERVED == 1]), stat = st))
  for (tier in names(tfs)) {
    tf <- tfs[[tier]]
    f  <- file.path(pooled_dir, tf$file)
    if (!file.exists(f)) { warning("missing pooled file: ", tf$file)
                           next }
    cv <- get_cv(readRDS(f))
    v  <- spon_ev[[tf$var]]
    x0 <- 0
    gr <- seq(min(v, na.rm = TRUE), max(v, na.rm = TRUE), length.out = 200)
    sx <- if (rescale) minmax01(v) else identity
    df <- rbind(df,
      curve_df(cv$cf, cv$V, gr, x0, tf$var, tf$sq, tier, st, sx))
  }
}

df$stat <- factor(df$stat, levels = names(stat_specs),
                  labels = vapply(stat_specs, function(s) s$lab, ""))
rug_df$stat <- factor(rug_df$stat, levels = names(stat_specs),
                      labels = vapply(stat_specs, function(s) s$lab, ""))
df$tier <- factor(df$tier, levels = names(tier_files(stat_specs$group)))

tier_cols <- c("No controls" = "grey65", "Full model" = "grey30",
               "Full model (non-linear)" = "black")
tier_lwd  <- c("No controls" = 0.6, "Full model" = 0.9,
               "Full model (non-linear)" = 1.2)
tier_lt   <- c("No controls" = "dashed", "Full model" = "solid",
               "Full model (non-linear)" = "solid")

p <- ggplot(df, aes(x = x, y = R, colour = tier, fill = tier,
                    linewidth = tier, linetype = tier)) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey40") +
  geom_ribbon(aes(ymin = lo, ymax = hi, colour = NULL), alpha = 0.12) +
  geom_line() +
  geom_rug(data = rug_df, aes(x = x),
           inherit.aes = FALSE, sides = "b", alpha = 0.25,
           colour = "grey30") +
  facet_wrap(~ stat, nrow = 2, scales = "free") +
  scale_colour_manual(values = tier_cols, name = NULL) +
  scale_fill_manual(values = tier_cols, name = NULL) +
  scale_linewidth_manual(values = tier_lwd, name = NULL) +
  scale_linetype_manual(values = tier_lt, name = NULL) +
  scale_y_log10() +
  labs(x = "Value of the dyadic statistic",
       y = "Relative propensity (vs no ties/overlap)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank())

fig_dir <- FIGURES_OUT
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(fig_dir, "relative_propensity_dyad.pdf"), p,
       width = 8, height = 6)
ggsave(file.path(fig_dir, "relative_propensity_dyad.png"), p,
       width = 8, height = 6, dpi = 300)
cat("Saved", file.path(fig_dir, "relative_propensity_dyad.pdf"), "\n")
