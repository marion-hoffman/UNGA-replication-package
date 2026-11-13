# Master script to run 16 parallel pipelines and pool results with Rubin's rules
# Windows-friendly (PSOCK clusters)

suppressPackageStartupMessages({
  library(parallel)
  library(data.table)
})

source(file.path(Find(function(d) file.exists(file.path(d, "paths.R")),
                      c(".", "..", "../..")), "paths.R"), chdir = TRUE)
self_dir <- SPONSORING_FUNCTIONS_DIR
results_dir <- SPONSORING_OUT
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# Source function files
source(file.path(self_dir, "functions_create_eventlist.R"))
source(file.path(self_dir, "functions_precompute_prev_senders.R"))
source(file.path(self_dir, "functions_add_ego_attribute_effects.R"))
source(file.path(self_dir, "functions_add_alter_attribute_effects.R"))
source(file.path(self_dir, "functions_add_dyad_attribute_effects.R"))
source(file.path(self_dir, "functions_add_ego_endo_effects.R"))
source(file.path(self_dir, "functions_add_alter_endo_effects.R"))
source(file.path(self_dir, "functions_add_dyad_endo_effects.R"))
source(file.path(self_dir, "functions_add_interaction_effects.R"))
source(file.path(self_dir, "functions_estimate_models.R"))
source(file.path(self_dir, "functions_pool_rubin.R"))

# Data paths
drafts_list_path <- DRAFTS_KEEP
country_draft_abs_path <- COUNTRY_DRAFT_ABS

# Draft-level data
recur_title_rds <- RECUR_TITLE_RDS
consensus_rds <- CONSENSUS_RDS
topic_rds <- TOPIC_RDS
nrefs_rds <- NREFS_RDS
recur_dyad_rds <- RECUR_DYAD_RDS

# Country-level data
gdp_csv <- GDP_CSV
conflict_csv <- CONFLICT_CSV
hr_csv <- HR_FARISS_CSV
hr_ciri_csv <- HR_CIRI_CSV
staff_csv <- STAFF_CSV
group_csv <- GROUP_CSV
regime_csv <- REGIME_CSV

# Dyad-level data
same_group_csv <- SAME_GROUP_CSV
same_regime_csv <- SAME_REGIME_CSV
diplo_csv <- DIPLO_CSV
trade_csv <- TRADE_CSV

# Fail fast if inputs are unreachable -- inside workers a missing file only
# surfaces as 'cannot open the connection'.
input_files <- c(drafts_list_path, country_draft_abs_path,
                 recur_title_rds, consensus_rds, topic_rds, nrefs_rds,
                 recur_dyad_rds, gdp_csv, conflict_csv, hr_csv, hr_ciri_csv,
                 staff_csv, group_csv, regime_csv,
                 same_group_csv, same_regime_csv, diplo_csv, trade_csv)
missing_inputs <- input_files[!file.exists(input_files)]
cat("Working directory:", getwd(), "\n")
cat("package root:", ROOT, "\n")
if (length(missing_inputs)) {
  stop("Input files not found (check wd and synced data):\n  ",
       paste(missing_inputs, collapse = "\n  "))
}

# One pipeline run
run_one <- function(rep_id, seed_base = 1000) {
  set.seed(seed_base + rep_id)
  rep_dir <- file.path(results_dir, sprintf("rep_%02d", rep_id))
  if (!dir.exists(rep_dir)) dir.create(rep_dir, recursive = TRUE)

  eventlist_csv <- file.path(rep_dir, "eventlist_for_rem.csv")
  prev_senders_rds <- file.path(rep_dir, "prev_senders_list.RDS")
  model_rds <- file.path(rep_dir, "model_result.RDS")

  # 1) Create eventlist
  create_eventlist(
    N_NON_EVENTS = 20,
    drafts_list_path = drafts_list_path,
    country_draft_abs_path = country_draft_abs_path,
    out_csv_path = eventlist_csv,
    seed = seed_base + rep_id
  )
  cat(sprintf("[rep %02d] done: eventlist\n", rep_id))

  # 2) Precompute previous senders (also emits prev_senders_strict_list.RDS)
  precompute_prev_senders(
    all_events_csv = eventlist_csv,
    country_draft_abs_path = country_draft_abs_path,
    out_rds_path = prev_senders_rds
  )
  cat(sprintf("[rep %02d] done: previous senders\n", rep_id))

  # 3) Add effects (attribute + endogeneity + interaction)
  add_ego_attribute_effects(
    all_events_csv = eventlist_csv,
    gdp_csv = gdp_csv,
    conflict_csv = conflict_csv,
    hr_csv = hr_csv,
    hr_ciri_csv = hr_ciri_csv,
    staff_csv = staff_csv,
    group_csv = group_csv,
    regime_csv = regime_csv,
    out_csv_path = eventlist_csv
  )
  cat(sprintf("[CPU %d, rep %02d] done: ego attributes\n", Sys.getpid(), rep_id))

  add_alter_attribute_effects(
    all_events_csv = eventlist_csv,
    country_draft_abs_path = country_draft_abs_path,
    recur_title_rds = recur_title_rds,
    consensus_rds = consensus_rds,
    topic_rds = topic_rds,
    nrefs_rds = nrefs_rds,
    out_csv_path = eventlist_csv
  )
  cat(sprintf("[rep %02d] done: alter attributes\n", rep_id))

  add_dyad_attribute_effects(
    all_events_csv = eventlist_csv,
    group_csv = same_group_csv,
    regime_csv = same_regime_csv,
    diplo_csv = diplo_csv,
    trade_csv = trade_csv,
    hr_ciri_csv = hr_ciri_csv,
    hr_fariss_csv = hr_csv,
    recur_dyad_rds = recur_dyad_rds,
    country_draft_abs_path = country_draft_abs_path,
    out_csv_path = eventlist_csv
  )
  cat(sprintf("[rep %02d] done: dyad attributes\n", rep_id))

  add_ego_endo_effects(
    all_events_csv = eventlist_csv,
    country_draft_abs_path = country_draft_abs_path,
    out_csv_path = eventlist_csv
  )
  cat(sprintf("[rep %02d] done: ego endo\n", rep_id))

  add_alter_endo_effects(
    all_events_csv = eventlist_csv,
    country_draft_abs_path = country_draft_abs_path,
    prev_senders_rds = prev_senders_rds,
    out_csv_path = eventlist_csv
  )
  cat(sprintf("[rep %02d] done: alter endo\n", rep_id))

  add_dyad_endo_effects(
    all_events_csv = eventlist_csv,
    country_draft_abs_path = country_draft_abs_path,
    prev_senders_rds = prev_senders_rds,
    out_csv_path = eventlist_csv
  )
  cat(sprintf("[rep %02d] done: dyad endo\n", rep_id))

  add_interaction_effects(
    all_events_csv = eventlist_csv,
    out_csv_path = eventlist_csv
  )
  cat(sprintf("[rep %02d] done: interactions\n", rep_id))

  # 4) Estimate models (placeholder)
  estimate_models(
    all_events_csv = eventlist_csv,
    out_rds_path = model_rds
  )
  cat(sprintf("[rep %02d] done: estimation\n", rep_id))

  return(list(rep_id = rep_id, model_rds = model_rds))
}

# Parallel execution
n_reps <- 16#8
rep_ids <- seq_len(n_reps)

n_cores <- min(detectCores(), n_reps)
cl <- makeCluster(n_cores, type = "PSOCK")
# Export needed symbols and source functions on workers
clusterExport(cl, varlist = c(
  "results_dir", "drafts_list_path", "country_draft_abs_path",
  "recur_title_rds", "consensus_rds", "topic_rds", "nrefs_rds", "recur_dyad_rds",
  "gdp_csv", "conflict_csv", "hr_csv", "hr_ciri_csv", "staff_csv", "group_csv", "regime_csv",
  "same_group_csv", "same_regime_csv", "diplo_csv", "trade_csv",
  "run_one", "self_dir"
), envir = environment())
clusterEvalQ(cl, {
  source(file.path(self_dir, "functions_create_eventlist.R"))
  source(file.path(self_dir, "functions_precompute_prev_senders.R"))
  source(file.path(self_dir, "functions_add_ego_attribute_effects.R"))
  source(file.path(self_dir, "functions_add_alter_attribute_effects.R"))
  source(file.path(self_dir, "functions_add_dyad_attribute_effects.R"))
  source(file.path(self_dir, "functions_add_ego_endo_effects.R"))
  source(file.path(self_dir, "functions_add_alter_endo_effects.R"))
  source(file.path(self_dir, "functions_add_dyad_endo_effects.R"))
  source(file.path(self_dir, "functions_add_interaction_effects.R"))
  source(file.path(self_dir, "functions_estimate_models.R"))
  NULL
})

results <- parLapply(cl, rep_ids, function(i) {
  tryCatch(run_one(i), error = function(e)
    sprintf("[rep %02d] %s | in: %s", i, conditionMessage(e),
            paste(deparse(conditionCall(e)), collapse = " ")))
})
rep_errors <- vapply(results, is.character, logical(1))
if (any(rep_errors)) {
  cat(sprintf("\n%d reps failed; errors:\n", sum(rep_errors)))
  cat(paste(unique(unlist(results[rep_errors])), collapse = "\n"), "\n")
}
stopCluster(cl)

# Load model results and pool per model
results_by_model <- list()
for (i in rep_ids) {
  rpath <- file.path(results_dir, sprintf("rep_%02d", i), "model_result.RDS")
  if (!file.exists(rpath)) next
  res_i <- readRDS(rpath)
  if (!is.list(res_i)) next
  for (mname in names(res_i)) {
    results_by_model[[mname]] <- c(results_by_model[[mname]], list(res_i[[mname]]))
  }
}

pooled_dir <- file.path(results_dir, "pooled")
if (!dir.exists(pooled_dir)) dir.create(pooled_dir, recursive = TRUE)

pooled_index <- list()
pooled_models <- list()
aic_by_model <- list()
r2_by_model <- list()
ll_by_model <- list()
if (length(results_by_model) == 0) {
  cat("No model results found to pool.\n")
} else {
  for (mname in names(results_by_model)) {
    # Transform stored fits (coxph) into list of {coef, vcov}
    raw_list <- results_by_model[[mname]]
    rlist <- list()
    for (obj in raw_list) {
      if (is.null(obj)) next
      if (inherits(obj, c("coxph", "clogit"))) {
        cf <- tryCatch(stats::coef(obj), error = function(e) NULL)
        vc <- tryCatch(stats::vcov(obj), error = function(e) NULL)
        if (!is.null(cf) && !is.null(vc)) {
          rlist[[length(rlist)+1]] <- list(coef = cf, vcov = vc)
        }
      } else if (is.list(obj) && !is.null(obj$coef) && !is.null(obj$vcov)) {
        # Backward compatibility
        rlist[[length(rlist)+1]] <- list(coef = obj$coef, vcov = obj$vcov)
      }
    }
    if (length(rlist) == 0) next
    pooled <- tryCatch(pool_rubin(rlist), error = function(e) e)
    # Pseudo-R2 is only meaningful within a rep (risk sets differ across
    # reps), so carry the per-rep values and their mean as attributes.
    pr2 <- vapply(raw_list, function(o) {
      v <- attr(o, "pseudo_r2"); if (is.null(v)) NA_real_ else v
    }, numeric(1))
    attr(pooled, "pseudo_r2_mean") <- mean(pr2, na.rm = TRUE)
    attr(pooled, "pseudo_r2_reps") <- pr2
    aic_i <- vapply(raw_list, function(o)
      tryCatch(stats::AIC(o), error = function(e) NA_real_), numeric(1))
    aic_by_model[[mname]] <- aic_i
    r2_by_model[[mname]] <- pr2
    ll_i <- vapply(raw_list, function(o)
      tryCatch(o$loglik[2], error = function(e) NA_real_), numeric(1))
    ll_by_model[[mname]] <- ll_i
    cat(sprintf("  %-36s pseudo-R2 mean %.4f | mean AIC %.0f | mean logLik %.1f (%d reps)\n",
                mname, mean(pr2, na.rm = TRUE), mean(aic_i, na.rm = TRUE),
                mean(ll_i, na.rm = TRUE), sum(!is.na(pr2))))
    out_path <- file.path(pooled_dir, paste0("pooled_", mname, ".RDS"))
    saveRDS(pooled, out_path)
    pooled_index[[mname]] <- out_path
    pooled_models[[mname]] <- pooled
  }
  saveRDS(pooled_index, file.path(pooled_dir, "pooled_index.RDS"))
  # Save pooled models as a single list, analogous to per-rep model_result.RDS
  saveRDS(pooled_models, file.path(pooled_dir, "pooled_models.RDS"))
  cat("Saved pooled results to ", pooled_dir, "\n", sep = "")

  # Cross-model fit summary over reps. Absolute AIC drifts with each rep's
  # sampled controls, so the informative column is the mean WITHIN-REP
  # difference vs the full model -- clogit specs only, since the tied coxph
  # models use merged risk sets and are on a different likelihood scale.
  mn <- names(aic_by_model)
  if (length(mn) > 0) {
    ref_aic <- aic_by_model[["model_all_effects"]]
    # Case-control specs = stratified fits (formerly clogit, now
    # stratified coxph for robust SEs) -- detect via the strata term.
    has_strata <- function(o) {
      sp <- tryCatch(attr(stats::terms(o), "specials"),
                     error = function(e) NULL)
      !is.null(sp$strata)
    }
    is_clogit <- vapply(mn, function(m) any(vapply(results_by_model[[m]],
      function(o) !is.null(o) && inherits(o, c("clogit", "coxph")) &&
        has_strata(o), logical(1))),
      logical(1))
    fit_tab <- data.frame(
      mean_AIC = vapply(mn, function(m) mean(aic_by_model[[m]], na.rm = TRUE), numeric(1)),
      mean_dAIC_vs_full = vapply(mn, function(m)
        if (is_clogit[m] && !is.null(ref_aic) &&
            length(aic_by_model[[m]]) == length(ref_aic))
          mean(aic_by_model[[m]] - ref_aic, na.rm = TRUE) else NA_real_,
        numeric(1)),
      mean_logLik = vapply(mn, function(m) mean(ll_by_model[[m]], na.rm = TRUE), numeric(1)),
      mean_pseudoR2 = vapply(mn, function(m) mean(r2_by_model[[m]], na.rm = TRUE), numeric(1)),
      sd_pseudoR2 = vapply(mn, function(m) stats::sd(r2_by_model[[m]], na.rm = TRUE), numeric(1)),
      row.names = mn)
    cat("\nFit statistics over reps (dAIC: clogit models vs model_all_effects):\n")
    print(round(fit_tab, 4))
  }
}

cat("Parallel runs finished. Results in:", results_dir, "\n")
