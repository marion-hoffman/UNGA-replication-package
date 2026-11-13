# Estimate the initiation model: which country creates a draft.
# Standard Cox formulation: at each event time, all states are in the risk
# set (full enumeration -> no sampling offset). Simultaneous creation events
# form a tied risk set handled by Efron (main) / Breslow (check). No
# per-event strata or randomization of concurrent events: every covariate
# is exogenous to event ordering.

source(file.path(Find(function(d) file.exists(file.path(d, "paths.R")),
                      c(".", "..", "../..")), "paths.R"), chdir = TRUE)

all_events <- read.csv(CREATION_EVENTLIST, stringsAsFactors = FALSE)

predictor_cols <- setdiff(colnames(all_events), c("IS_OBSERVED", "SOURCE", "TIME", "ORDER",
                                                  "num_possible", "num_sampled", "weight", "offset"))
cat("Using the following predictors:\n")
print(predictor_cols)

library(survival)

# Print summary + AIC + McFadden pseudo-R2 (1 - LL_full/LL_null)
print_fit <- function(fit) {
  print(summary(fit))
  cat("AIC:", round(stats::AIC(fit), 1),
      "  McFadden pseudo-R2:", round(1 - fit$loglik[2] / fit$loglik[1], 4), "\n")
}

model_data <- all_events[, c("TIME", "IS_OBSERVED", predictor_cols)]

# MAIN HR measure = CIRI: keep Fariss under ego_hr_fariss, let ego_hr
# carry the CIRI score (sq/dummy/interaction terms then all use CIRI).
model_data$ego_hr_fariss <- model_data$ego_hr
model_data$ego_hr <- model_data$ego_hr_ciri

# Diagnostic: distribution of OBSERVED events across the four
# (ego_condemned, alter_topic_hr) cells. The condemned x HR cell is
# (near-)empty, which separated the ego_condemned:alter_topic_hr term --
# hence it is NOT included in any spec below. The zero count itself is the
# substantive finding and should be reported descriptively (e.g., "no
# previously-condemned state initiated an HR-condemnation draft").
obs <- model_data[model_data$IS_OBSERVED == 1, ]
cat("\nObserved events by condemned x HR-topic cell:\n")
print(table(`condemned` = obs$ego_condemned, `HR draft` = obs$alter_topic_hr))

# Quadratic HR-score terms -- RAW squares, identical to the sponsoring
# specs so coefficients are comparable across the two models. The
# "diversion" hypothesis predicts the initiation hazard rises steeply
# at the LOW end of the HR scale -- a convex decreasing function of ego_hr
# on the non-HR margin (ego_hr_sq > 0), with ego_hr_sq:alter_topic_hr < 0
# (HR-topic advantage concentrated at the top of the scale).
model_data$ego_hr_sq <- model_data$ego_hr^2
model_data$ego_hr_fariss_sq <- model_data$ego_hr_fariss^2

# Polynomial terms for the nonlinear robustness spec (mirror of the
# sponsoring f_nl: cubics on GDP and the endogenous sponsor count)
model_data$ego_log_gdp_sq      <- model_data$ego_log_gdp^2
model_data$ego_log_gdp_cb      <- model_data$ego_log_gdp^3
model_data$ego_sponsor1year_sq <- model_data$ego_sponsor1year^2
model_data$ego_sponsor1year_cb <- model_data$ego_sponsor1year^3

# Binary alternatives: bottom- vs top-quartile HR scores (middle = baseline).
# Diversion hypothesis: ego_hr_low > 0 AND ego_hr_high > 0 on the non-HR
# margin (U-shape); ego_hr_low:alter_topic_hr <= 0 and
# ego_hr_high:alter_topic_hr > 0 (only high-HR states favor HR topics).
q <- quantile(model_data$ego_hr, c(0.25, 0.75), na.rm = TRUE)
model_data$ego_hr_low  <- as.numeric(model_data$ego_hr <= q[1])
model_data$ego_hr_high <- as.numeric(model_data$ego_hr >= q[2])
cat("\nCIRI cutoffs -- low <=", round(q[1], 2), "| high >=", round(q[2], 2), "\n")
q <- quantile(model_data$ego_hr_fariss, c(0.25, 0.75), na.rm = TRUE)
model_data$ego_hr_fariss_low  <- as.numeric(model_data$ego_hr_fariss <= q[1])
model_data$ego_hr_fariss_high <- as.numeric(model_data$ego_hr_fariss >= q[2])
cat("Fariss cutoffs -- low <=", round(q[1], 2), "| high >=", round(q[2], 2), "\n")

# Risk-set cells are (state, topic): alter_topic_hr distinguishes creating
# an HR-condemnation draft vs another draft; ego_hr:alter_topic_hr tests
# whether a country's HR record shifts it toward creating HR drafts
# specifically (ego main effects = propensity to initiate any draft).
#
# Diversion hypothesis: violators (low HR) create drafts at high rates but
# on non-HR topics -> convex decreasing curve on the non-HR margin
# (ego_hr_sq > 0) with ego_hr_sq:alter_topic_hr < 0 (HR-topic advantage
# concentrated at the top of the scale). The dummy variant predicts
# ego_hr_low > 0 & ego_hr_high > 0 (U-shape), ego_hr_low:topic <= 0 and
# ego_hr_high:topic > 0 (only high-HR states favor HR topics).

# === Specifications (mirror the sponsoring robustness suite) =================
# Main: quadratic HR-score shape + topic interaction.
f_main <- Surv(TIME, IS_OBSERVED) ~ ego_log_gdp + ego_conflict +
  ego_hr + ego_hr_sq + ego_staff + ego_group + #ego_condemned +
  ego_sponsor1year + alter_topic_hr +
  ego_hr:alter_topic_hr + ego_hr_sq:alter_topic_hr

# No endogenous effect: main spec minus ego_sponsor1year (mirror of the
# no-structural-effects sponsoring spec).
f_no_endo <- update(f_main, . ~ . - ego_sponsor1year)

# Robustness: alternative HR measure (Fariss instead of CIRI) in the main
# effect and the interaction -- mirror of the sponsoring *_hr_fariss spec.
f_hr_fariss <- Surv(TIME, IS_OBSERVED) ~ ego_log_gdp + ego_conflict +
  ego_hr_fariss + ego_hr_fariss_sq + ego_staff + ego_group +
  ego_sponsor1year + alter_topic_hr +
  ego_hr_fariss:alter_topic_hr + ego_hr_fariss_sq:alter_topic_hr

# Robustness: bottom/top-quartile HR dummies (mid-range reference)
# replacing the continuous linear+quadratic terms -- discrete check on
# the shape, for each measure (mirror of the sponsoring *_hr_dummies).
f_hr_dummies <- Surv(TIME, IS_OBSERVED) ~ ego_log_gdp + ego_conflict +
  ego_hr_low + ego_hr_high + ego_staff + ego_group +
  ego_sponsor1year + alter_topic_hr +
  ego_hr_low:alter_topic_hr + ego_hr_high:alter_topic_hr

f_hr_fariss_dummies <- Surv(TIME, IS_OBSERVED) ~ ego_log_gdp + ego_conflict +
  ego_hr_fariss_low + ego_hr_fariss_high + ego_staff + ego_group +
  ego_sponsor1year + alter_topic_hr +
  ego_hr_fariss_low:alter_topic_hr + ego_hr_fariss_high:alter_topic_hr

# Robustness: nonlinear extension (mirror of the sponsoring f_nl: cubics
# on GDP and the endogenous sponsor count, quadratic HR shape).
f_nl <- Surv(TIME, IS_OBSERVED) ~
  ego_log_gdp + ego_log_gdp_sq + ego_log_gdp_cb + ego_conflict +
  ego_hr + ego_hr_sq + ego_staff + ego_group +
  ego_sponsor1year + ego_sponsor1year_sq + ego_sponsor1year_cb +
  alter_topic_hr + ego_hr:alter_topic_hr + ego_hr_sq:alter_topic_hr

f_nl_hr_fariss <- Surv(TIME, IS_OBSERVED) ~
  ego_log_gdp + ego_log_gdp_sq + ego_log_gdp_cb + ego_conflict +
  ego_hr_fariss + ego_hr_fariss_sq + ego_staff + ego_group +
  ego_sponsor1year + ego_sponsor1year_sq + ego_sponsor1year_cb +
  alter_topic_hr + ego_hr_fariss:alter_topic_hr +
  ego_hr_fariss_sq:alter_topic_hr

# === Minimal hypothesis specs ("univariate" tests) =====================
# Step-1 descriptive evidence, mirror of the sponsoring
# estimate_univariate_effects.R suite: the unconditional association of
# each hypothesized effect before attribute / endogeneity controls.
# Only H1/H2 blocks exist here -- initiation riskStuy,ehare nraw,rssmO
#:asr_sq ego_hr_fariss_sq are MEAN-CENTERED (defined above), unlike the
# sponsoring specs.
f_univ <- list(
  # H1 (HR shape): quadratic block only
  univariate_h1_hr_quadratic = Surv(TIME, IS_OBSERVED) ~
    ego_hr + ego_hr_sq,
  # H2 (HR x topic): HR block + topic + both interactions (hierarchy:
  # the squared interaction needs the linear one present)
  univariate_h2_hr_x_topic = Surv(TIME, IS_OBSERVED) ~
    ego_hr + ego_hr_sq + alter_topic_hr +
    ego_hr:alter_topic_hr + ego_hr_sq:alter_topic_hr,
  # Robustness: HR as low/high quartile dummies (discrete shape check)
  univariate_h1_hr_dummies = Surv(TIME, IS_OBSERVED) ~
    ego_hr_low + ego_hr_high,
  univariate_h2_hr_dum_x_topic = Surv(TIME, IS_OBSERVED) ~
    ego_hr_low + ego_hr_high + alter_topic_hr +
    ego_hr_low:alter_topic_hr + ego_hr_high:alter_topic_hr,
  # Robustness: Fariss human-rights measure in place of CIRI
  univariate_h1_hr_fariss = Surv(TIME, IS_OBSERVED) ~
    ego_hr_fariss + ego_hr_fariss_sq,
  univariate_h2_hr_fariss_x_topic = Surv(TIME, IS_OBSERVED) ~
    ego_hr_fariss + ego_hr_fariss_sq + alter_topic_hr +
    ego_hr_fariss:alter_topic_hr + ego_hr_fariss_sq:alter_topic_hr)

run_fit <- function(fml, ties = "efron")
  coxph(fml, data = model_data, ties = ties)

fits <- list(
  model_main              = run_fit(f_main),              # quadratic CIRI (main)
  model_no_endo           = run_fit(f_no_endo),           # minus endogenous term
  model_hr_fariss         = run_fit(f_hr_fariss),         # alt measure (Fariss)
  model_hr_dummies        = run_fit(f_hr_dummies),        # CIRI quartile dummies
  model_hr_fariss_dummies = run_fit(f_hr_fariss_dummies), # Fariss quartile dummies
  model_nl_effects        = run_fit(f_nl),                # nonlinear extension
  model_nl_hr_fariss      = run_fit(f_nl_hr_fariss),      # nonlinear + Fariss
  model_main_breslow      = run_fit(f_main, "breslow")    # tie-handling check
)

# Minimal-spec fits (no controls): unconditional H1/H2 associations on
# the same risk sets. Names match the sponsoring univariate files.
univ_fits <- lapply(names(f_univ), function(nm) {
  fit <- tryCatch(run_fit(f_univ[[nm]]), error = function(e) e)
  if (inherits(fit, "error")) {
    warning("univariate fit failed: ", nm, " -- ", conditionMessage(fit))
    NULL
  } else fit
})
names(univ_fits) <- names(f_univ)
fits <- c(fits, univ_fits)

for (nm in names(fits)) {
  cat("\n===", nm, "===\n")
  if (is.null(fits[[nm]])) { cat("(fit failed -- skipped)\n"); next }
  print_fit(fits[[nm]])
}

saveRDS(fits, file.path(CREATION_OUT, "initiating_model_fits.RDS"))
cat("\nFits saved to initiating_model_fits.RDS\n")

# Also save each spec as list(coef, vcov, m = 1) -- the pooled-object
# format used by the delta tests (test_h1_marginal_effects.R) and by
# get_cv() in plot_relative_propensity.R. m = 1: single eventlist, no
# replication pooling on this side.
for (nm in names(fits)) {
  fit <- fits[[nm]]
  if (is.null(fit)) next
  saveRDS(list(coef = stats::coef(fit), vcov = stats::vcov(fit), m = 1L),
          file.path(CREATION_OUT, sprintf("initiating_%s.RDS", nm)))
}
cat("Per-spec coef/vcov objects saved to initiating_*.RDS\n")


