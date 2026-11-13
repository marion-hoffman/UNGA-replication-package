# plot_relative_propensity.R
#
# Relative propensity curves R(x,z) = exp{ eta(x,z) - eta(x0,z) },
# x0 = median ego_hr -- predicted propensity of a country with HR score
# x vs an otherwise identical median-score country (NOT a probability).
# All non-HR terms cancel in the contrast; alter_topic_hr's main effect
# is constant within a column and cancels too.
#
# Figure: rows = outcome (draft creation | sponsoring),
#         cols = draft scope (all drafts | HR-condemnation only),
#         lines = no controls / full model / full model (non-linear),
#         ribbons = 95% band from (pooled) vcov, rug = ego_hr of
#         observed events, log-scaled y so 1 is the reference.
#
# Sponsoring tiers come from the Rubin-pooled objects; the initiating
# tiers are fitted here on the initiating eventlist (that pipeline has
# no replicate pooling). All specs use RAW squares (center = 0).
#
# Run:
#   Rscript plot_relative_propensity_hr.R
# Optionally override the sponsoring results folder:
#   RESULTS_DIR="path/to/results" Rscript plot_relative_propensity_hr.R

suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
})

source(file.path(Find(function(d) file.exists(file.path(d, "paths.R")),
                      c(".", "..", "../..")), "paths.R"), chdir = TRUE)
results_dir <- SPONSORING_OUT
pooled_dir <- POOLED_DIR
init_csv <- CREATION_EVENTLIST
stopifnot(dir.exists(pooled_dir), file.exists(init_csv))

# ---- generic contrast -------------------------------------------------
# Delta(x,z) weights for the HR block of a coef vector. 'center' is the
# centering constant of the squared term (0 everywhere now: both models
# use raw squares). Terms absent from the spec contribute nothing.
hr_weights <- function(nm, x, x0, z, center = 0) {
  find <- function(a, b) {
    if (paste0(a, ":", b) %in% nm) return(paste0(a, ":", b))
    if (paste0(b, ":", a) %in% nm) return(paste0(b, ":", a))
    NA_character_
  }
  w <- setNames(rep(0, length(nm)), nm)
  d_lin <- x - x0
  d_sq  <- (x - center)^2 - (x0 - center)^2
  if ("ego_hr" %in% nm)    w["ego_hr"] <- d_lin
  if ("ego_hr_sq" %in% nm) w["ego_hr_sq"] <- d_sq
  if (z == 1) {
    i1 <- find("ego_hr", "alter_topic_hr")
    i2 <- find("ego_hr_sq", "alter_topic_hr")
    if (!is.na(i1)) w[i1] <- d_lin
    if (!is.na(i2)) w[i2] <- d_sq
  }
  w
}

# Curves over an x grid for one (coef, vcov) pair and topic z
curve_df <- function(cf, V, xgrid, z, x0, center, tier, outcome,
                     topic = z) {
  nm <- names(cf)
  out <- lapply(xgrid, function(x) {
    w <- hr_weights(nm, x, x0, z, center)
    d <- sum(w * cf)
    s <- sqrt(drop(t(w) %*% V %*% w))
    data.frame(x = x, eta = d, se = s,
               R = exp(d),
               lo = exp(d - 1.96 * s), hi = exp(d + 1.96 * s))
  })
  cbind(do.call(rbind, out),
        topic = topic, tier = tier, outcome = outcome)
}

# coef/vcov extractor for pooled lists and coxph fits
get_cv <- function(obj) {
  if (inherits(obj, c("coxph", "clogit")))
    return(list(cf = stats::coef(obj), V = stats::vcov(obj)))
  list(cf = obj$coef, V = obj$vcov)   # pooled list(coef, vcov, m)
}

# ================= sponsoring row (pooled objects) ====================
# NOTE: main HR measure is CIRI -- pooled 'ego_hr' coefs hold CIRI
# coefficients, so the x axis uses the ego_hr_ciri column.
rep1 <- file.path(results_dir, "rep_01", "eventlist_for_rem.csv")
spon_ev <- utils::read.csv(if (file.exists(rep1)) rep1 else init_csv)
x0_spon <- stats::median(spon_ev$ego_hr_ciri, na.rm = TRUE)
spon_grid <- seq(min(spon_ev$ego_hr_ciri, na.rm = TRUE),
                 max(spon_ev$ego_hr_ciri, na.rm = TRUE), length.out = 200)

# Column 1 = ALL drafts: propensity over every draft, i.e. the two HR
# main parameters only (z = 0). Column 2 = HR-condemnation drafts: same
# spec at z = 1, adding the topic-HR interactions. Descriptive tier uses
# the h2 univariate spec (HR block + alter_topic_hr only).
spon_files <- list(
  descriptive = "pooled_univariate_h2_hr_x_topic.RDS",
  full        = "pooled_model_all_effects.RDS",
  nl          = "pooled_model_nl_effects.RDS")

# Score landmarks for the console report (same distribution the grid spans)
xq_spon <- stats::quantile(spon_ev$ego_hr_ciri, c(0, .25, .75, 1),
                           na.rm = TRUE, names = FALSE)
report_df <- NULL
add_report <- function(cf, V, xs, z, x0, center, tier, outcome, topic) {
  d <- curve_df(cf, V, xs, z, x0, center, tier, outcome, topic = topic)
  d$q <- c("min", "Q25", "Q75", "max")[seq_len(nrow(d))]
  report_df <<- rbind(report_df, d)
}

spon_df <- NULL
for (tier in names(spon_files)) {
  f <- file.path(pooled_dir, spon_files[[tier]])
  f <- f[file.exists(f)]
  if (length(f) == 0) { warning("no pooled file for sponsoring ", tier)
                        next }
  cv <- get_cv(readRDS(f[1]))
  spon_df <- rbind(spon_df,
    curve_df(cv$cf, cv$V, spon_grid, 0, x0_spon, 0, tier, "Sponsoring",
             topic = "all"),
    curve_df(cv$cf, cv$V, spon_grid, 1, x0_spon, 0, tier, "Sponsoring",
             topic = "hr"))
  add_report(cv$cf, cv$V, xq_spon, 0, x0_spon, 0, tier, "Sponsoring", "all")
  add_report(cv$cf, cv$V, xq_spon, 1, x0_spon, 0, tier, "Sponsoring", "hr")
}

# ================= initiating row (refit, single eventlist) ===========
init_ev <- utils::read.csv(init_csv)
# MAIN HR measure = CIRI (same swap as estimate_model.R)
init_ev$ego_hr <- init_ev$ego_hr_ciri
init_ev$ego_hr_sq <- init_ev$ego_hr^2
# polynomial controls of the non-linear spec (f_nl in estimate_model.R)
init_ev$ego_log_gdp_sq      <- init_ev$ego_log_gdp^2
init_ev$ego_log_gdp_cb      <- init_ev$ego_log_gdp^3
init_ev$ego_sponsor1year_sq <- init_ev$ego_sponsor1year^2
init_ev$ego_sponsor1year_cb <- init_ev$ego_sponsor1year^3
init_center <- 0   # raw square (same as sponsoring + initiating specs)
x0_init <- stats::median(init_ev$ego_hr, na.rm = TRUE)
init_grid <- seq(min(init_ev$ego_hr, na.rm = TRUE),
                 max(init_ev$ego_hr, na.rm = TRUE), length.out = 200)

# Same structure as sponsoring: one interacted spec per tier, evaluated
# at z = 0 (all drafts: two HR main params) and z = 1 (HR-condemnation).
f_init <- list(
  # descriptive: HR block + topic + interactions only (h2 analogue)
  descriptive = Surv(TIME, IS_OBSERVED) ~ ego_hr + ego_hr_sq +
    alter_topic_hr + ego_hr:alter_topic_hr + ego_hr_sq:alter_topic_hr,
  # full: attributes + endogenous count
  full = Surv(TIME, IS_OBSERVED) ~ ego_log_gdp + ego_conflict +
    ego_hr + ego_hr_sq + ego_staff + ego_group +
    ego_sponsor1year +
    alter_topic_hr + ego_hr:alter_topic_hr + ego_hr_sq:alter_topic_hr,
  # non-linear: f_nl from initiating-analysis/estimate_model.R
  nl = Surv(TIME, IS_OBSERVED) ~
    ego_log_gdp + ego_log_gdp_sq + ego_log_gdp_cb + ego_conflict +
    ego_hr + ego_hr_sq + ego_staff + ego_group +
    ego_sponsor1year + ego_sponsor1year_sq + ego_sponsor1year_cb +
    alter_topic_hr + ego_hr:alter_topic_hr + ego_hr_sq:alter_topic_hr)

xq_init <- stats::quantile(init_ev$ego_hr, c(0, .25, .75, 1),
                           na.rm = TRUE, names = FALSE)

init_df <- NULL
for (tier in names(f_init)) {
  fit <- tryCatch(coxph(f_init[[tier]], data = init_ev, ties = "efron"),
                  error = function(e) NULL)
  if (is.null(fit)) { warning("initiating fit failed: ", tier); next }
  cv <- get_cv(fit)
  init_df <- rbind(init_df,
    curve_df(cv$cf, cv$V, init_grid, 0, x0_init, init_center,
             tier, "Draft creation", topic = "all"),
    curve_df(cv$cf, cv$V, init_grid, 1, x0_init, init_center,
             tier, "Draft creation", topic = "hr"))
  add_report(cv$cf, cv$V, xq_init, 0, x0_init, init_center,
             tier, "Draft creation", "all")
  add_report(cv$cf, cv$V, xq_init, 1, x0_init, init_center,
             tier, "Draft creation", "hr")
}

# ================= assemble + plot ====================================
df <- rbind(spon_df, init_df)
tier_levels <- c("descriptive", "full", "nl")
tier_labels <- c("No controls", "Full model", "Full model (non-linear)")
df$tier    <- factor(df$tier, levels = tier_levels, labels = tier_labels)
df$topic   <- factor(df$topic, levels = c("all", "hr"),
                     labels = c("All drafts", "HR-condemnation drafts"))
df$outcome <- factor(df$outcome,
                     levels = c("Draft creation", "Sponsoring"),
                     labels = c("Draft creation", "Draft sponsoring"))

# rug: ego_hr of OBSERVED events (same in both topic columns)
rug_df <- rbind(
  data.frame(x = spon_ev$ego_hr_ciri[spon_ev$IS_OBSERVED == 1],
             outcome = "Sponsoring"),
  data.frame(x = init_ev$ego_hr[init_ev$IS_OBSERVED == 1],
             outcome = "Draft creation"))
rug_df$outcome <- factor(rug_df$outcome,
                         levels = c("Draft creation", "Sponsoring"),
                         labels = c("Draft creation", "Draft sponsoring"))

tier_cols <- c("No controls" = "grey65", "Full model" = "grey30",
               "Full model (non-linear)" = "black")
tier_lwd  <- c("No controls" = 0.6, "Full model" = 0.9,
               "Full model (non-linear)" = 1.2)
tier_lt   <- c("No controls" = "dashed", "Full model" = "solid",
               "Full model (non-linear)" = "solid")

p <- ggplot(df, aes(x = x, y = R, colour = tier, fill = tier,
                    linewidth = tier, linetype = tier)) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey40") +
  geom_ribbon(aes(ymin = lo, ymax = hi, colour = NULL),
              alpha = 0.12) +
  geom_line() +
  geom_rug(data = rug_df, aes(x = x),
           inherit.aes = FALSE, sides = "b", alpha = 0.25,
           colour = "grey30") +
  facet_grid(outcome ~ topic, scales = "free_y") +  # per-row y range
  scale_colour_manual(values = tier_cols, name = NULL) +
  scale_fill_manual(values = tier_cols, name = NULL) +
  scale_linewidth_manual(values = tier_lwd, name = NULL) +
  scale_linetype_manual(values = tier_lt, name = NULL) +
  scale_y_log10() +
  labs(x = "Human-rights score (CIRI)",
       y = "Relative propensity (vs median score)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank())

fig_dir <- FIGURES_OUT
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(fig_dir, "relative_propensity_hr.pdf"), p,
       width = 8, height = 6)
ggsave(file.path(fig_dir, "relative_propensity_hr.png"), p,
       width = 8, height = 6, dpi = 300)
cat("Saved", file.path(fig_dir, "relative_propensity_hr.pdf"), "\n")

# ---- console report: R at score extremes/quartiles (for the text) ----
# R(x, z) vs the median score x0, 95% CI from the same delta-method vcov
# as the ribbons. Quartiles of the plotted score distribution per row.
if (!is.null(report_df)) {
  report_df$topic <- factor(report_df$topic, levels = c("all", "hr"),
                            labels = c("All drafts",
                                       "HR-condemnation drafts"))
  report_df$tier  <- factor(report_df$tier,
                            levels = tier_levels, labels = tier_labels)
  cat("\nRelative propensity vs median score (95% CI):\n")
  for (oc in unique(report_df$outcome)) {
    cat(sprintf("\n== %s ==\n", oc))
    for (tp in levels(report_df$topic)) {
      for (tr in levels(report_df$tier)) {
        sub <- report_df[report_df$outcome == oc &
                         report_df$topic == tp &
                         report_df$tier == tr, ]
        if (!nrow(sub)) next
        cat(sprintf("  %s | %s\n", tp, tr))
        for (i in seq_len(nrow(sub)))
          cat(sprintf("    %-3s (x = %5.2f): R = %7.3f  [%6.3f, %6.3f]\n",
                      sub$q[i], sub$x[i], sub$R[i],
                      sub$lo[i], sub$hi[i]))
      }
    }
  }
}

