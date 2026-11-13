estimate_models <- function(
  all_events_csv,
  out_rds_path
) {
  # Load data
  all_events <- read.csv(all_events_csv, stringsAsFactors = FALSE)
  # Build model dataset
  predictor_cols <- setdiff(colnames(all_events), c("IS_OBSERVED", "SOURCE", "TARGET", "TIME", "ORDER",
                                                    "num_possible", "num_sampled", "weight", "offset"))
  # NAs are handled by coxph via na.action; here keep as-is, except alter_cons often set to 0
  if ("alter_cons" %in% names(all_events)) {
    all_events$alter_cons[is.na(all_events$alter_cons)] <- 0
  }
  model_data <- all_events[, c("TIME", "IS_OBSERVED", "ORDER", predictor_cols, "offset")]

  # MAIN HR measure = CIRI. Swap: keep the Fariss values under
  # ego_hr_fariss for the alternative-measure robustness spec, then let
  # ego_hr (and every downstream _sq / dummy / interaction term) carry
  # the CIRI score -- all formulas below are unchanged.
  if (all(c("ego_hr", "ego_hr_ciri") %in% names(model_data))) {
    model_data$ego_hr_fariss <- model_data$ego_hr
    model_data$ego_hr <- model_data$ego_hr_ciri
  }

  # Polynomial terms for the non-linear robustness specs (mirrors the
  # fml_alt block in sponsoring-analysis/run_diagnostics_rem.R).
  for (v in c("ego_log_gdp", "ego_hr", "ego_hr_fariss", "alter_sponsor1year",
              "dyad_group_prop", "dyad_regime_prop", "dyad_diplo_av",
              "dyad_trade_gtmean", "dyad_av_diff_hr",
              "dyad_av_diff_hr_fariss", "dyad_av_diff_hr_strict",
              "dyad_av_diff_hr_fariss_strict", "alter_hr_av",
              "alter_hr_av_fariss", "alter_hr_av_strict",
              "alter_hr_av_fariss_strict"))
    if (!is.null(model_data[[v]]))
      model_data[[paste0(v, "_sq")]] <- model_data[[v]]^2
  for (v in c("ego_log_gdp", "alter_sponsor1year"))
    model_data[[paste0(v, "_cb")]] <- model_data[[v]]^3

  # HR quantile dummies for the discrete-shape robustness spec: bottom
  # and top quartile of ego_hr, mid-range is the reference category.
  q_hr <- stats::quantile(model_data$ego_hr, c(.25, .75),
                          na.rm = TRUE, names = FALSE)
  model_data$ego_hr_low  <- as.integer(model_data$ego_hr <= q_hr[1])
  model_data$ego_hr_high <- as.integer(model_data$ego_hr >= q_hr[2])

  # Same dummies on the Fariss score (alternative HR measure).
  if ("ego_hr_fariss" %in% names(model_data)) {
    q_hr_f <- stats::quantile(model_data$ego_hr_fariss, c(.25, .75),
                              na.rm = TRUE, names = FALSE)
    model_data$ego_hr_fariss_low <-
      as.integer(model_data$ego_hr_fariss <= q_hr_f[1])
    model_data$ego_hr_fariss_high <-
      as.integer(model_data$ego_hr_fariss >= q_hr_f[2])
  }

  # Fit models
  requireNamespace("survival", quietly = TRUE)
  Surv <- survival::Surv
  coxph <- survival::coxph
  clogit <- survival::clogit
  strata <- survival::strata

  # Define model formulas
  f_all <- IS_OBSERVED ~ ego_log_gdp + ego_conflict + ego_hr + ego_hr_sq +
    ego_staff + ego_group + ego_sponsor1year +
    alter_topic_hr + alter_topic_coop + alter_topic_ispa + alter_nrefs +
    ego_hr:alter_topic_hr + ego_hr_sq:alter_topic_hr +
    alter_recur + alter_cons + alter_sponsor1year + alter_hr_av +
    dyad_group_prop + dyad_regime_prop +
    dyad_diplo_av + dyad_trade_gtmean + dyad_av_diff_hr +
    dyad_recur + dyad_4cycle_sponsor1year +
    strata(ORDER)

  # Attributes + interpretable DRAFT-HISTORY terms only (recurring
  # draft, lagged sponsor consensus, ego repeat of a recurring draft).
  # The structural endogenous terms -- draft sponsor momentum
  # (alter_sponsor1year) and dyad 4-cycle transitivity -- enter only in
  # f_all, so the contrast isolates the hard-to-interpret structural
  # effects.
  f_no_past <- IS_OBSERVED ~ ego_log_gdp + ego_conflict + ego_hr + ego_hr_sq +
    ego_staff + ego_group +
    alter_topic_hr + alter_topic_coop + alter_topic_ispa + alter_nrefs +
    ego_hr:alter_topic_hr + ego_hr_sq:alter_topic_hr +
    alter_recur + alter_cons + alter_hr_av +
    dyad_group_prop + dyad_regime_prop +
    dyad_diplo_av + dyad_trade_gtmean + dyad_av_diff_hr +
    dyad_recur + strata(ORDER)

  f_ego_only <- IS_OBSERVED ~ ego_log_gdp + ego_conflict + ego_hr + ego_hr_sq +
    ego_staff + ego_group + strata(ORDER)

  f_ego_only_plus_endo <- IS_OBSERVED ~ ego_log_gdp + ego_conflict + ego_hr + ego_hr_sq +
    ego_staff + ego_group + ego_sponsor1year + strata(ORDER)

  f_alter_only <- IS_OBSERVED ~ alter_topic_hr + alter_topic_coop + alter_topic_ispa + alter_nrefs +
    strata(ORDER)

  f_alter_only_plus_endo <- IS_OBSERVED ~ alter_topic_hr + alter_topic_coop + alter_topic_ispa + alter_nrefs +
    alter_recur + alter_cons + alter_sponsor1year + strata(ORDER)

  f_ego_alter <- IS_OBSERVED ~ ego_log_gdp + ego_conflict + ego_hr + ego_hr_sq +
    ego_staff + ego_group + alter_topic_hr +
    alter_topic_coop + alter_topic_ispa + alter_nrefs + 
    ego_hr:alter_topic_hr + ego_hr_sq:alter_topic_hr + strata(ORDER)

  f_ego_alter_plus_endo <- IS_OBSERVED ~ ego_log_gdp + ego_conflict + ego_hr + ego_hr_sq +
    ego_staff + ego_group + ego_sponsor1year +
    alter_topic_hr + alter_topic_coop + alter_topic_ispa + alter_nrefs +
    alter_recur + alter_cons + alter_sponsor1year +
    ego_hr:alter_topic_hr + ego_hr_sq:alter_topic_hr + strata(ORDER)

  # Robustness: full model with alternative human-rights measure
  # (Fariss instead of CIRI) in the main effect and the interaction.
  f_all_hr_fariss <- IS_OBSERVED ~ ego_log_gdp + ego_conflict +
    ego_hr_fariss + ego_hr_fariss_sq + ego_staff + ego_group + ego_sponsor1year +
    alter_topic_hr + alter_topic_coop + alter_topic_ispa + alter_nrefs +
    ego_hr_fariss:alter_topic_hr + ego_hr_fariss_sq:alter_topic_hr +
    alter_recur + alter_cons + alter_sponsor1year + alter_hr_av_fariss +
    dyad_group_prop + dyad_regime_prop +
    dyad_diplo_av + dyad_trade_gtmean + dyad_av_diff_hr_fariss +
    dyad_recur + dyad_4cycle_sponsor1year +
    strata(ORDER)

  # Robustness: ego human rights as low/high DUMMIES (bottom and top
  # quartile of ego_hr; mid-range reference) replacing the continuous
  # linear+quadratic terms -- a discrete check on the HR shape and on
  # its conditioning of alter_topic_hr.
  f_hr_dummies <- IS_OBSERVED ~ ego_log_gdp + ego_conflict +
    ego_hr_low + ego_hr_high +
    ego_staff + ego_group + ego_sponsor1year +
    alter_topic_hr + alter_topic_coop + alter_topic_ispa + alter_nrefs +
    ego_hr_low:alter_topic_hr + ego_hr_high:alter_topic_hr +
    alter_recur + alter_cons + alter_sponsor1year + alter_hr_av +
    dyad_group_prop + dyad_regime_prop +
    dyad_diplo_av + dyad_trade_gtmean + dyad_av_diff_hr +
    dyad_recur + dyad_4cycle_sponsor1year +
    strata(ORDER)

  # Same dummy spec with the Fariss measure: completes the
  # CIRI/Fariss x continuous/dummies robustness grid.
  f_hr_fariss_dummies <- IS_OBSERVED ~ ego_log_gdp + ego_conflict +
    ego_hr_fariss_low + ego_hr_fariss_high +
    ego_staff + ego_group + ego_sponsor1year +
    alter_topic_hr + alter_topic_coop + alter_topic_ispa + alter_nrefs +
    ego_hr_fariss_low:alter_topic_hr + ego_hr_fariss_high:alter_topic_hr +
    alter_recur + alter_cons + alter_sponsor1year + alter_hr_av_fariss +
    dyad_group_prop + dyad_regime_prop +
    dyad_diplo_av + dyad_trade_gtmean + dyad_av_diff_hr_fariss +
    dyad_recur + dyad_4cycle_sponsor1year +
    strata(ORDER)

  # Robustness: full model treating concurrent events as SIMULTANEOUS
  # (all events at a shared TIME compete in one tied risk set), via coxph
  # with Efron and Breslow tie handling. Strict-precedence effects
  # (TIME < time): concurrent events cannot precede each other, matching
  # the simultaneity assumption. ego_sponsor1year and dyad_recur have no
  # strict versions -- kept as is.
  f_all_tied <- Surv(TIME, IS_OBSERVED) ~ ego_log_gdp + ego_conflict + ego_hr + ego_hr_sq +
    ego_staff + ego_group + ego_sponsor1year +
    alter_topic_hr + alter_topic_coop + alter_topic_ispa + alter_nrefs +
    ego_hr:alter_topic_hr + ego_hr_sq:alter_topic_hr +
    alter_recur + alter_cons + alter_sponsor1year_strict +
    alter_hr_av_strict +
    dyad_group_prop_strict + dyad_regime_prop_strict +
    dyad_diplo_av_strict + dyad_trade_gtmean_strict +
    dyad_av_diff_hr_strict + dyad_recur +
    dyad_4cycle_sponsor1year_strict

  # Robustness: trade as raw mean importance over all previous sponsors
  # (zeros included) instead of the count above the matrix mean.
  f_trade_raw <- IS_OBSERVED ~ ego_log_gdp + ego_conflict + ego_hr + ego_hr_sq +
    ego_staff + ego_group + ego_sponsor1year +
    alter_topic_hr + alter_topic_coop + alter_topic_ispa + alter_nrefs +
    ego_hr:alter_topic_hr + ego_hr_sq:alter_topic_hr +
    alter_recur + alter_cons + alter_sponsor1year + alter_hr_av +
    dyad_group_prop + dyad_regime_prop +
    dyad_diplo_av + dyad_trade_av + dyad_av_diff_hr +
    dyad_recur + dyad_4cycle_sponsor1year +
    strata(ORDER)

  # Robustness: trade as the COUNT of previous sponsors with strictly
  # positive trade dependence on the ego.
  f_trade_pos <- IS_OBSERVED ~ ego_log_gdp + ego_conflict + ego_hr + ego_hr_sq +
    ego_staff + ego_group + ego_sponsor1year +
    alter_topic_hr + alter_topic_coop + alter_topic_ispa + alter_nrefs +
    ego_hr:alter_topic_hr + ego_hr_sq:alter_topic_hr +
    alter_recur + alter_cons + alter_sponsor1year + alter_hr_av +
    dyad_group_prop + dyad_regime_prop +
    dyad_diplo_av + dyad_trade_pos + dyad_av_diff_hr +
    dyad_recur + dyad_4cycle_sponsor1year +
    strata(ORDER)

  # Robustness: non-linear extension selected by the spline diagnostics
  # (fml_alt in run_diagnostics_rem.R): cubics for ego_log_gdp and
  # alter_sponsor1year, quadratics for ego_hr and the dyad
  # proportions/averages, plus the ego_hr_sq:alter_topic_hr interaction.
  # 'all' unanimity dummies omitted, as in the diagnostics spec.
  f_nl <- IS_OBSERVED ~ ego_log_gdp + ego_log_gdp_sq + ego_log_gdp_cb +
    ego_conflict + ego_hr + ego_hr_sq +
    ego_staff + ego_group + ego_sponsor1year +
    alter_topic_hr + alter_topic_coop + alter_topic_ispa + alter_nrefs +
    ego_hr:alter_topic_hr + ego_hr_sq:alter_topic_hr +
    alter_recur + alter_cons +
    alter_sponsor1year + alter_sponsor1year_sq + alter_sponsor1year_cb +
    alter_hr_av + alter_hr_av_sq +
    dyad_group_prop + dyad_group_prop_sq +
    dyad_regime_prop + dyad_regime_prop_sq +
    dyad_diplo_av + dyad_diplo_av_sq +
    dyad_trade_gtmean + dyad_trade_gtmean_sq +
    dyad_av_diff_hr + dyad_av_diff_hr_sq +
    dyad_recur + dyad_4cycle_sponsor1year +
    strata(ORDER)

  # Same non-linear spec with the alternative human-rights measure
  # (Fariss instead of CIRI): is the HR pattern robust to the indicator?
  f_nl_hr_fariss <- IS_OBSERVED ~ ego_log_gdp + ego_log_gdp_sq + ego_log_gdp_cb +
    ego_conflict + ego_hr_fariss + ego_hr_fariss_sq +
    ego_staff + ego_group + ego_sponsor1year +
    alter_topic_hr + alter_topic_coop + alter_topic_ispa + alter_nrefs +
    ego_hr_fariss:alter_topic_hr + ego_hr_fariss_sq:alter_topic_hr +
    alter_recur + alter_cons +
    alter_sponsor1year + alter_sponsor1year_sq + alter_sponsor1year_cb +
    alter_hr_av_fariss + alter_hr_av_fariss_sq +
    dyad_group_prop + dyad_group_prop_sq +
    dyad_regime_prop + dyad_regime_prop_sq +
    dyad_diplo_av + dyad_diplo_av_sq +
    dyad_trade_gtmean + dyad_trade_gtmean_sq +
    dyad_av_diff_hr_fariss + dyad_av_diff_hr_fariss_sq +
    dyad_recur + dyad_4cycle_sponsor1year +
    strata(ORDER)

  # Robustness: full model with the case-control OFFSET
  # (offset = log(num_possible/num_sampled) -- a Horvitz-Thompson / Monte
  # Carlo approximation of the full-risk-set conditional likelihood). The
  # main specs above omit the offset: under standard nested case-control
  # sampling the conditional likelihood over the sampled choice set needs
  # no correction for varying original risk-set sizes. Near-identical
  # estimates would confirm the control sampling is not driving results.
  f_all_htoffset <- IS_OBSERVED ~ ego_log_gdp + ego_conflict + ego_hr + ego_hr_sq +
    ego_staff + ego_group + ego_sponsor1year +
    alter_topic_hr + alter_topic_coop + alter_topic_ispa + alter_nrefs +
    ego_hr:alter_topic_hr + ego_hr_sq:alter_topic_hr +
    alter_recur + alter_cons + alter_sponsor1year + alter_hr_av +
    dyad_group_prop + dyad_regime_prop +
    dyad_diplo_av + dyad_trade_gtmean + dyad_av_diff_hr +
    dyad_recur + dyad_4cycle_sponsor1year +
    offset(offset) + strata(ORDER)

  # Robustness of the homophily terms (alter_hr_av + dyad_av_diff_hr,
  # part of f_all): a quadratic extension of the diff term, on CIRI
  # and on the fully-Fariss spec.
  f_all_hrdiff_sq <-
    update(f_all, . ~ . + dyad_av_diff_hr_sq)
  f_all_hrdiff_sq_fariss <-
    update(f_all_hr_fariss, . ~ . + dyad_av_diff_hr_fariss_sq)



  run_model <- function(fml, use = c("clogit", "coxph"), ties = "efron") {
    use <- match.arg(use)
    # robust = TRUE: sandwich (dfbeta) SEs -- the case-control likelihood is a
    # pseudo-likelihood on sampled risk sets, so model-based SEs are invalid.
    # vcov(fit) then returns the sandwich variance used by the pooling step.
    # clogit's exact partial likelihood cannot return dfbeta residuals, so
    # robust=TRUE is impossible there ("dfbeta residuals are not available
    # for the exact method"). Each ORDER stratum has exactly one observed
    # event -- no within-stratum ties -- so a stratified coxph with Efron
    # gives numerically identical coefficients AND supports sandwich SEs.
    model <- if (use == "clogit") {
      coxph(stats::update(fml, Surv(TIME, .) ~ .),
            data = model_data, ties = "efron", robust = TRUE)
    } else {
      coxph(fml, data = model_data, ties = ties, robust = TRUE)
    }
    fit <- tryCatch(
      model,
      error = function(e) e)
    if (inherits(fit, "error")) {
      attr(fit, "fit_error") <- conditionMessage(fit)
      return(NULL)
    }
    # McFadden pseudo-R2: 1 - LL_full/LL_null (null = all coefs zero, same
    # strata/offset). Comparable across specs fit on the same data.
    attr(fit, "pseudo_r2") <- 1 - fit$loglik[2] / fit$loglik[1]
    fit
  }

  results <- list(
    model_all_effects = run_model(f_all),
    model_all_effects_htoffset = run_model(f_all_htoffset),
    model_no_past = run_model(f_no_past),
    model_ego_only = run_model(f_ego_only),
    model_ego_only_plus_endo = run_model(f_ego_only_plus_endo),
    model_alter_only = run_model(f_alter_only),
    model_alter_only_plus_endo = run_model(f_alter_only_plus_endo),
    model_ego_alter = run_model(f_ego_alter),
    model_ego_alter_plus_endo = run_model(f_ego_alter_plus_endo),
    model_all_effects_hr_fariss = run_model(f_all_hr_fariss),
    model_all_effects_hr_dummies = run_model(f_hr_dummies),
    model_all_effects_hr_fariss_dummies = run_model(f_hr_fariss_dummies),
    model_all_effects_trade_raw = run_model(f_trade_raw),
    model_all_effects_trade_pos = run_model(f_trade_pos),
    model_nl_effects = run_model(f_nl),
    model_nl_hr_fariss = run_model(f_nl_hr_fariss),

    model_all_effects_hrdiff_sq = run_model(f_all_hrdiff_sq),
    model_all_effects_hrdiff_sq_fariss =
      run_model(f_all_hrdiff_sq_fariss),
    model_all_effects_tied_efron = run_model(f_all_tied, "coxph", "efron"),
    model_all_effects_tied_breslow = run_model(f_all_tied, "coxph", "breslow")
  )

  dir.create(dirname(out_rds_path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(results, out_rds_path)
  invisible(out_rds_path)
}
