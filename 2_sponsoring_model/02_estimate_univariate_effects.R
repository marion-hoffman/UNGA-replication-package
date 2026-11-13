# Minimal hypothesis specs on the SAME case-control risk sets as the
# main models: one stratified coxph per hypothesis block, same
# strata(ORDER), Efron ties, robust SEs, then Rubin-pooled across reps.
#
# This is the paper's "step 1" descriptive evidence: the unconditional
# association of each hypothesized effect before attribute / endogeneity
# controls enter. Comparing these to model_all_effects shows how much
# of each signal is absorbed by the conditioning terms.
#
# Blocks keep jointly-needed terms together:
#   H1 (HR shape)        : ego_hr + ego_hr_sq
#   H2 (HR x topic)      : both HR terms + alter_topic_hr + both
#                          interactions (hierarchy: the squared
#                          interaction needs the linear one present)
#   H3 (dyad affinity)   : one spec per measure; the two extra trade
#                          parametrizations are H3 robustness variants
#   H1/H2 robustness     : same blocks with HR quartile dummies
#                          (<=25% / >=75%, mid-range reference) and
#                          with the Fariss measure in place of CIRI
#                          (main HR measure is now CIRI)
#
# Output:
#   results_dir/pooled/pooled_univariate_<spec>.RDS  (one per spec,
#     same list(coef, vcov, m) format -- make_latex_tables.R renders
#     them automatically since they match the pooled_*.RDS pattern)
#   console table grouped by spec
#
# Run:
#   Rscript estimate_univariate_effects.R
# Optionally override the results folder:
#   RESULTS_DIR=/path/to/results Rscript estimate_univariate_effects.R

suppressPackageStartupMessages(library(survival))

source(file.path(Find(function(d) file.exists(file.path(d, "paths.R")),
                      c(".", "..", "../..")), "paths.R"), chdir = TRUE)
results_dir <- Sys.getenv(
  "RESULTS_DIR",
  SPONSORING_OUT)

source(file.path(SPONSORING_FUNCTIONS_DIR, "functions_pool_rubin.R"))

Surv <- survival::Surv
coxph <- survival::coxph
strata <- survival::strata

rep_dirs <- list.dirs(results_dir, full.names = TRUE, recursive = FALSE)
rep_dirs <- rep_dirs[grepl("^rep_\\d+$", basename(rep_dirs))]
if (length(rep_dirs) == 0)
  stop("No rep_XX folders found in ", results_dir)

# One minimal spec per hypothesis (names become the output labels).
univ_specs <- list(
  h1_hr_quadratic = "ego_hr + ego_hr_sq",
  h2_hr_x_topic   = "ego_hr + ego_hr_sq + alter_topic_hr +
                     ego_hr:alter_topic_hr + ego_hr_sq:alter_topic_hr",
  h3_group        = "dyad_group_prop",
  h3_regime       = "dyad_regime_prop",
  h3_diplo        = "dyad_diplo_av",
  h3_trade        = "dyad_trade_gtmean",
  h3_trade_raw    = "dyad_trade_av",   # robustness: raw mean importance
  h3_trade_pos    = "dyad_trade_pos",  # robustness: any positive tie
  # Robustness: HR as low/high quartile dummies (discrete shape check)
  h1_hr_dummies     = "ego_hr_low + ego_hr_high",
  h2_hr_dum_x_topic = "ego_hr_low + ego_hr_high + alter_topic_hr +
                       ego_hr_low:alter_topic_hr + ego_hr_high:alter_topic_hr",
  # Robustness: Fariss human-rights measure in place of CIRI
  h1_hr_fariss        = "ego_hr_fariss + ego_hr_fariss_sq",
  h2_hr_fariss_x_topic = "ego_hr_fariss + ego_hr_fariss_sq + alter_topic_hr +
                          ego_hr_fariss:alter_topic_hr + ego_hr_fariss_sq:alter_topic_hr"
)

fit_one <- function(model_data, rhs) {
  fml <- stats::as.formula(
    paste("Surv(TIME, IS_OBSERVED) ~", rhs, "+ strata(ORDER)"))
  fit <- tryCatch(
    coxph(fml, data = model_data, ties = "efron", robust = TRUE),
    error = function(e) e)
  if (inherits(fit, "error")) return(NULL)
  cf <- tryCatch(stats::coef(fit), error = function(e) NULL)
  vc <- tryCatch(stats::vcov(fit), error = function(e) NULL)
  if (is.null(cf) || is.null(vc) || length(cf) == 0) return(NULL)
  list(coef = cf, vcov = vc)
}

# spec -> list of {coef, vcov} per rep
by_spec <- setNames(vector("list", length(univ_specs)), names(univ_specs))

for (rep_dir in sort(rep_dirs)) {
  csv <- file.path(rep_dir, "eventlist_for_rem.csv")
  if (!file.exists(csv)) next
  cat(sprintf("[%s] minimal-spec fits...\n", basename(rep_dir)))
  model_data <- utils::read.csv(csv)
  # MAIN HR measure = CIRI: keep Fariss under ego_hr_fariss, let ego_hr
  # carry the CIRI score (sq/dummy terms below then all use CIRI).
  if (all(c("ego_hr", "ego_hr_ciri") %in% names(model_data))) {
    model_data$ego_hr_fariss <- model_data$ego_hr
    model_data$ego_hr <- model_data$ego_hr_ciri
  }
  # Derived terms needed by the hypothesis blocks (same construction as
  # functions_estimate_models.R)
  if (!"ego_hr_sq" %in% names(model_data))
    model_data$ego_hr_sq <- model_data$ego_hr^2
  if (!"ego_hr_fariss_sq" %in% names(model_data) &&
      "ego_hr_fariss" %in% names(model_data))
    model_data$ego_hr_fariss_sq <- model_data$ego_hr_fariss^2
  if (!all(c("ego_hr_low", "ego_hr_high") %in% names(model_data))) {
    q_hr <- stats::quantile(model_data$ego_hr, c(.25, .75),
                            na.rm = TRUE, names = FALSE)
    model_data$ego_hr_low  <- as.integer(model_data$ego_hr <= q_hr[1])
    model_data$ego_hr_high <- as.integer(model_data$ego_hr >= q_hr[2])
  }
  for (spec in names(univ_specs)) {
    # skip specs whose columns are absent from this eventlist
    vars_needed <- all.vars(stats::as.formula(paste("~", univ_specs[[spec]])))
    if (!all(vars_needed %in% names(model_data))) next
    r <- fit_one(model_data, univ_specs[[spec]])
    if (!is.null(r)) by_spec[[spec]] <- c(by_spec[[spec]], list(r))
  }
}

pooled_dir <- file.path(results_dir, "pooled")
dir.create(pooled_dir, showWarnings = FALSE, recursive = TRUE)

stars <- function(p) ifelse(p < .001, "***", ifelse(p < .01, "**",
                     ifelse(p < .05, "*", ifelse(p < .1, ".", ""))))

tab <- data.frame()
for (spec in names(univ_specs)) {
  rlist <- by_spec[[spec]]
  if (length(rlist) < 2) {
    cat(sprintf("  %-18s skipped (only %d non-null fits)\n",
                spec, length(rlist)))
    next
  }
  pooled <- pool_rubin(rlist)
  s <- summary_pooled(pooled)
  saveRDS(pooled,
          file.path(pooled_dir, paste0("pooled_univariate_", spec, ".RDS")))
  tab <- rbind(tab, data.frame(
    spec = spec, term = s$term, estimate = s$estimate,
    std.error = s$std.error, z = s$statistic, p.value = s$p.value,
    sig = stars(s$p.value), m = pooled$m, row.names = NULL))
}

cat("\nMinimal-spec pooled estimates (same risk sets as main specs):\n")
tab$estimate <- round(tab$estimate, 3)
tab$std.error <- round(tab$std.error, 3)
tab$z <- round(tab$z, 2)
tab$p.value <- signif(tab$p.value, 3)
print(tab, right = FALSE)
cat("\nSaved per-spec pooled objects to", pooled_dir,
    "-- run make_latex_tables.R to render them.\n")
