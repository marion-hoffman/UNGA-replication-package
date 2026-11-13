# 02_sensitivity_imputations_controls.R
# Sensitivity of the main sponsoring model (model_all_effects) to
#   1) the number of imputed event sequences: re-pools subsets of the
#      per-rep fits with m = 4, 8, 12, 16 (plus the imputation burden);
#   2) the number of sampled non-events per event: refits the model on
#      each rep with k = 5, 10, 15, 20 controls and pools.
#
# Reads:  output/sponsoring/rep_XX/{model_result.RDS, eventlist_for_rem.csv}
# Writes: output/post_estimation_checks/sensitivity_m_and_controls/

## ---------------------------------------------------------------------------
## 0. Setup
## ---------------------------------------------------------------------------
source(file.path(Find(function(d) file.exists(file.path(d, "paths.R")),
                      c(".", "..", "../..")), "paths.R"), chdir = TRUE)
results_dir <- SPONSORING_OUT
out_dir     <- file.path(POST_CHECKS_OUT, "sensitivity_m_and_controls")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

source(file.path(SPONSORING_FUNCTIONS_DIR, "functions_pool_rubin.R"))  # pool_rubin, summary_pooled

FOCAL_MODEL <- "model_all_effects"
M_GRID      <- c(4, 8, 12, 16)

# Terms highlighted in the selection table (hypothesis-related effects)
HYP_TERMS <- c(
  "ego_sponsor1year", "alter_recur", "alter_cons", "alter_sponsor1year",
  "dyad_group_prop", "dyad_group_all", "dyad_regime_prop", "dyad_regime_all",
  "dyad_diplo_av", "dyad_trade_av", "dyad_recur", "dyad_4cycle_sponsor1year",
  "ego_hr", "alter_topic_hr", "ego_hr:alter_topic_hr"
)

# Pretty labels for LaTeX output
term_labels <- c(
  ego_log_gdp = "GDP (log)", ego_conflict = "Civil conflict",
  ego_hr = "HR score (CIRI)", ego_hr_sq = "HR score^2 (CIRI)",
  ego_hr_fariss = "HR score (Fariss)", ego_hr_fariss_sq = "HR score^2 (Fariss)",
  ego_staff = "Mission size", ego_group = "Group membership",
  ego_sponsor1year = "N prior sponsorings (ego)",
  alter_topic_hr = "Topic (HR)", alter_topic_coop = "Topic (coop)",
  alter_topic_ispa = "Topic (ISPA)", alter_nrefs = "N references",
  alter_recur = "Recurring (alter)", alter_cons = "Recurring & consensus (alter)",
  alter_sponsor1year = "N prior sponsors (alter)",
  alter_sponsor1year_strict = "N prior sponsors (strict version)",
  alter_hr_av = "Av HR score of prior sponsors",
  alter_hr_av_strict = "Av HR score of prior sponsors (strict version)",
  alter_hr_av_fariss = "Av HR score of prior sponsors (Fariss)",
  dyad_av_diff_hr = "Av abs HR difference",
  dyad_av_diff_hr_strict = "Av abs HR difference (strict version)",
  dyad_av_diff_hr_fariss = "Av abs HR difference (Fariss)",
  dyad_group_prop = "% shared groups", 
  dyad_regime_prop = "% same regime", 
  dyad_diplo_av = "Av diplomatic ties", 
  dyad_trade_gtmean = "% trade ties (> mean)", 
  dyad_recur = "Recurring sponsoring (dyad)",
  dyad_group_prop_strict = "% shared groups (strict version)",
  dyad_regime_prop_strict = "% same regime (strict version)",
  dyad_diplo_av_strict = "Av diplomatic ties (strict version)", 
  dyad_trade_gtmean_strict = "% trade ties (strict version)",
  dyad_4cycle_sponsor1year = "Repetition co-sponsoring",
  dyad_4cycle_sponsor1year_strict = "Repetition co-sponsoring (strict version)",
  `ego_hr:alter_topic_hr` = "HR (CIRI) x topic (HR)", 
  `ego_hr_sq:alter_topic_hr` = "HR^2 (CIRI) x topic (HR)", 
  `ego_hr_fariss:alter_topic_hr` = "HR (Fariss) x topic (HR)",
  `ego_hr_fariss_sq:alter_topic_hr` = "HR^2 (Fariss) x topic (HR)"
)
lab <- function(t) ifelse(t %in% names(term_labels), term_labels[t],
                          gsub("_", "\\_", t, fixed = TRUE))

## ---------------------------------------------------------------------------
## 1. Load per-rep fits and convert to {coef, vcov} lists
## ---------------------------------------------------------------------------
rep_files <- sort(Sys.glob(file.path(results_dir, "rep_*", "model_result.RDS")))
if (length(rep_files) == 0) stop("No per-rep model_result.RDS found in ", results_dir)
cat("Found", length(rep_files), "rep result files\n")

per_rep <- lapply(rep_files, readRDS)          # list over reps; each a named list of fits
model_names <- unique(unlist(lapply(per_rep, names)))

extract_cv <- function(fit) {
  if (is.null(fit) || !inherits(fit, "coxph")) return(NULL)
  cf <- tryCatch(stats::coef(fit), error = function(e) NULL)
  vc <- tryCatch(stats::vcov(fit), error = function(e) NULL)
  if (is.null(cf) || is.null(vc)) return(NULL)
  list(coef = cf, vcov = vc)
}
# fits_by_model[[mname]][[rep]] <- {coef, vcov}
fits_by_model <- lapply(model_names, function(mn)
  lapply(per_rep, function(res) extract_cv(res[[mn]])))
names(fits_by_model) <- model_names

## ---------------------------------------------------------------------------
## 2. Imputation burden: how much of the sequence is actually uncertain?
## ---------------------------------------------------------------------------
# Wave order inside a draft is exact; only times are drawn (uniform in each
# draft's (date_or, date_final] window). Order ambiguity therefore exists only
# between events on DIFFERENT drafts whose windows overlap.
#   7a) structural: share of event pairs whose relative order is not fixed
#   7b) realized:   agreement between the m imputed sequences (Kendall tau on
#                   event-time ranks; wave-mates count as tied)

cda_all <- readRDS(COUNTRY_DRAFT_ABS)
dk <- as.character(readRDS(DRAFTS_KEEP))
cda_imp <- cda_all[as.character(cda_all$ldoc_final) %in% dk, , drop = FALSE]
sc <- names(cda_imp)[30:225]
win_lo <- suppressWarnings(as.numeric(as.POSIXct(cda_imp$date_or, tz = "UTC")))
win_hi <- suppressWarnings(as.numeric(as.POSIXct(cda_imp$date_final, tz = "UTC")))

n_join <- integer(nrow(cda_imp)); n_wave <- integer(nrow(cda_imp))
for (r in seq_len(nrow(cda_imp))) {
  v <- suppressWarnings(as.numeric(as.matrix(cda_imp[r, sc, drop = FALSE])))
  w <- v[!is.na(v) & v >= 2]
  n_join[r] <- length(w); n_wave[r] <- if (length(w)) max(w) else 0L
}
n_ev <- sum(n_join); n_pairs <- choose(n_ev, 2)

same_draft <- sum(choose(n_join, 2))          # order fixed or simultaneous
ov <- outer(win_lo, win_hi, `<`) & outer(win_hi, win_lo, `>`)  # windows overlap
ov[is.na(ov)] <- FALSE
ov[upper.tri(ov, diag = TRUE)] <- FALSE       # count each draft pair once
ambiguous <- sum(ov * outer(n_join, n_join))  # cross-draft reorderable pairs

cat("\n=== 7. Imputation burden ===\n")
cat("We know WHICH state joined WHICH draft and in which wave, but not the\n")
cat("exact day. Each event time is drawn uniformly inside its draft's open\n")
cat("window (date_or, date_final], keeping wave order. So 'imputation' here\n")
cat("concerns event TIMING only -- the events themselves are observed.\n\n")

win_days <- (win_hi - win_lo) / 86400
cat(sprintf("Data: %d join events on %d drafts (the other %d drafts in the\n",
            n_ev, sum(n_join > 0), nrow(cda_imp) - sum(n_join > 0)))
cat(sprintf("      subset were only created, never joined). Mean %.1f waves per\n",
            mean(n_wave[n_wave > 0])))
cat(sprintf("      draft (max %d); drafts stay open a median of %.0f days.\n\n",
            max(n_wave), median(win_days[is.finite(win_days)])))

cat("For each pair of events, do we know which happened first?\n")
cat(sprintf("  Total event pairs: %s\n", format(n_pairs, big.mark = ",")))
cat(sprintf("    %5.1f%%  same draft -- order fixed by wave number (or simultaneous)\n",
            100 * same_draft / n_pairs))
cat(sprintf("    %5.1f%%  different drafts, windows never overlap -- order fixed\n",
            100 * (n_pairs - same_draft - ambiguous) / n_pairs))
cat(sprintf("    %5.1f%%  different drafts, windows DO overlap -- order imputed\n\n",
            100 * ambiguous / n_pairs))
cat(sprintf("So at most %.1f%% of event pairs can flip order across imputation\n",
            100 * ambiguous / n_pairs))
cat(sprintf("draws (%.0f%% of draft pairs have overlapping windows).\n",
            100 * sum(ov) / choose(nrow(cda_imp), 2)))

# 7b) realized distance between the m imputed sequences
rep_ev_files <- Sys.glob(file.path(results_dir, "rep_*", "eventlist_for_rem.csv"))
obs_list <- lapply(rep_ev_files, function(f) {
  e <- read.csv(f, stringsAsFactors = FALSE)
  e <- e[e$IS_OBSERVED == 1, ]
  data.frame(key = paste(e$SOURCE, e$TARGET, sep = "|"), TIME = e$TIME)
})
keys <- Reduce(intersect, lapply(obs_list, `[[`, "key"))
if (length(obs_list) > 1 && length(keys) > 0) {
  rank_mat <- sapply(obs_list, function(e)
    rank(e$TIME[match(keys, e$key)], ties.method = "average"))
  pos_sd <- apply(rank_mat, 1, sd) / length(keys)
  set.seed(1)
  rp <- combn(ncol(rank_mat), 2)
  if (ncol(rp) > 60) rp <- rp[, sample(ncol(rp), 60)]
  taus <- apply(rp, 2, function(j)
    cor(rank_mat[, j[1]], rank_mat[, j[2]], method = "kendall"))
  cat(sprintf("\nHow different are the %d imputed sequences from each other?\n",
              ncol(rank_mat)))
  cat(sprintf("  A given event's position in the sequence varies with sd %.3f\n",
              median(pos_sd)))
  cat(sprintf("  of the sequence length (median over events; 0 = never moves).\n"))
  cat(sprintf("  Kendall's tau between two sequences: mean %.3f, worst %.3f\n",
              mean(taus), min(taus)))
  cat(sprintf("  (1.0 = identical ordering). That corresponds to ~%.1f%% of\n",
              mean(taus)))
  cat(sprintf("  event pairs ordered differently between two draws -- below the\n"))
  cat(sprintf("  %.1f%% bound above, as it must be (only overlapping-window\n",
              100 * ambiguous / n_pairs))
  cat("  cross-draft pairs can ever flip; wave-mates are tied in both).\n")
}

cat("\nAll post-analysis outputs in:", out_dir, "\n")

## ---------------------------------------------------------------------------
## 3a. Stability across number of sequences m
## ---------------------------------------------------------------------------
cat("\n=== Stability across m sequences (", FOCAL_MODEL, ") ===\n", sep = "")

fits_focal <- Filter(Negate(is.null), fits_by_model[[FOCAL_MODEL]])
pooled_by_m <- lapply(M_GRID, function(m) {
  pool_rubin(fits_focal[seq_len(min(m, length(fits_focal)))])
})
names(pooled_by_m) <- paste0("m", M_GRID)

ref <- pooled_by_m[[length(pooled_by_m)]]              # largest m as reference
terms_all <- names(ref$coef)
stab <- data.frame(term = terms_all, stringsAsFactors = FALSE)
for (mn in names(pooled_by_m)) {
  p <- pooled_by_m[[mn]]
  stab[[paste0("b_", mn)]]  <- p$coef[terms_all]
  stab[[paste0("se_", mn)]] <- sqrt(diag(p$vcov))[terms_all]
}
stab$abs_diff_vs_ref <- abs(stab$b_m16 - ref$coef[terms_all])
stab$diff_over_se    <- stab$abs_diff_vs_ref / sqrt(diag(ref$vcov))[terms_all]

print(stab[order(-stab$diff_over_se),
           c("term", "b_m4", "b_m8", "b_m12", "b_m16", "diff_over_se")],
      digits = 3)
write.csv(stab, file.path(out_dir, "stability_by_m.csv"), row.names = FALSE)

# Global summary: how much do coefficients move as m grows? Same compact
# format as the control-set check (3b) below: max/median deviation in SE
# units, share of terms with same sign / same 5% significance, mean B/W.
se_ref_m <- sqrt(diag(ref$vcov))[terms_all]
b_ref_m  <- ref$coef[terms_all]
for (j in seq_along(M_GRID)[-length(M_GRID)]) {
  m   <- M_GRID[j]
  b_m <- stab[[paste0("b_m", m)]]
  s_m <- stab[[paste0("se_m", m)]]
  # B/W for the m-subset (SD across sequence fits / avg within-seq SE)
  fm <- fits_focal[seq_len(m)]
  bm <- do.call(rbind, lapply(fm, function(x) x$coef[terms_all]))
  Bm <- apply(bm, 2, stats::sd)
  Wm <- sqrt(diag(Reduce(`+`, lapply(fm, function(x)
    x$vcov[terms_all, terms_all, drop = FALSE])) / m))
  same_sign <- mean(sign(b_m) == sign(b_ref_m))
  sig_ref   <- abs(b_ref_m / se_ref_m) > 1.96
  same_sig  <- mean((abs(b_m / s_m) > 1.96) == sig_ref)
  kept_sig  <- if (any(sig_ref))
    mean((abs(b_m / s_m) > 1.96)[sig_ref]) else NA_real_
  d_m <- abs(b_m - b_ref_m) / se_ref_m
  cat(sprintf(paste0("m=%2d vs m=16: max D* = %.3f | median D* = %.3f | ",
                     "sign match %3.0f%% | sig match %3.0f%% | m16-sig kept %3.0f%% ",
                     "| mean B/W = %.2f\n"),
              m, max(d_m, na.rm = TRUE), median(d_m, na.rm = TRUE),
              100 * same_sign, 100 * same_sig, 100 * kept_sig,
              mean(Bm / Wm)))
}

# Between-sequence variance component B at full m (order-imputation +
# control-sampling variation); should shrink relative to W as m is adequate
m <- length(fits_focal)
betas <- do.call(rbind, lapply(fits_focal, function(r) r$coef[terms_all]))
B_sd <- apply(betas, 2, stats::sd)
W_sd <- sqrt(diag(Reduce(`+`, lapply(fits_focal, function(r)
  r$vcov[terms_all, terms_all, drop = FALSE])) / m))
stab$B_over_W <- B_sd[terms_all] / W_sd[terms_all]
cat("\nB/W (between- over within-sequence SD), top terms:\n")
print(stab[order(-stab$B_over_W), c("term", "B_over_W")], digits = 3)

# ---- Stability plots --------------------------------------------------
# Per-term small multiples of raw pooled coefficient + SE vs m (one page),
# plus a B/W ratio chart. Nested subsets -> successive points correlated;
# these are convergence diagnostics, not independent redraws.
ms     <- M_GRID
b_mat  <- as.matrix(stab[, paste0("b_",  "m", ms)])
se_mat <- as.matrix(stab[, paste0("se_", "m", ms)])

# Per-term detail as small multiples on a single page: raw pooled
# coefficient + SE at each m. Points = pooled estimate, whiskers = +-1 SE,
# dotted line = m16 reference. Convergence = estimate settles and whiskers
# stop widening.
nc  <- min(5L, nrow(stab))
nr  <- ceiling(nrow(stab) / nc)
pdf(file.path(out_dir, "stability_per_term.pdf"),
    width = 4.4 * nc, height = 3.6 * nr)
par(mfrow = c(nr, nc), mar = c(3.2, 3.6, 2.6, .6),
    oma = c(4.5, 4.5, 1, .5), mgp = c(2.4, .9, 0))
for (i in seq_len(nrow(stab))) {
  b  <- as.numeric(b_mat[i, ])
  se <- as.numeric(se_mat[i, ])
  ylim <- range(c(b - 1.96 * se, b + 1.96 * se))
  plot(ms, b, pch = 16, cex = 1.4, ylim = ylim, xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", main = lab(stab$term[i]), cex.main = 2.5)
  axis(1, at = ms, cex.axis = 1.3)
  axis(2, cex.axis = 2)
  arrows(ms, b - 1.96 * se, ms, b + 1.96 * se, angle = 90, code = 3,
         length = .06, lwd = 1.4)
  abline(h = 0, lty = 3, col = "grey50", lwd = 1.4)
  abline(h = b[length(b)], col = "red", lty = 1, lwd = 1.4)  # m16 ref
}
mtext("Number of imputed sequences m", 1, outer = TRUE, line = 2.4,
      cex = 2.5)
mtext("Pooled coefficient (95% CI)", 2, outer = TRUE, line = 2,
      cex = 2.5)
dev.off()
cat("Wrote", file.path(out_dir, "stability_per_term.pdf"), "\n")

# B/W ratios: how much of each coefficient's variance comes from sequence
# imputation vs model sampling noise. v = 1 is equal shares.
ord <- order(stab$B_over_W)
pdf(file.path(out_dir, "stability_bw.pdf"), width = 7,
    height = max(4, 0.28 * nrow(stab)))
par(mar = c(5, 14, 3, 1), mgp = c(3, .7, 0))
dotchart(stab$B_over_W[ord], labels = lab(stab$term[ord]),
         xlab = "B/W (between- over within-sequence SD)",
         main = "Imputation-driven share of coefficient noise",
         pch = 19, cex = .8, cex.lab = .9)
abline(v = 1, lty = 2, col = "grey40")
dev.off()
cat("Wrote", file.path(out_dir, "stability_bw.pdf"), "\n")


## ---------------------------------------------------------------------------
## 3b. Sensitivity to the number of sampled controls (pooled over reps)
## ---------------------------------------------------------------------------
# Within each imputed sequence, the sampled controls in every ORDER
# stratum are randomly ordered once; nested sets keep the first
# 5/10/15/20, so cross-k differences reflect control-set SIZE only. The
# focal model is refit per (rep, k) and Rubin-pooled over reps -- the
# k = 20 column should reproduce the published pooled estimates.
K_GRID <- c(5, 10, 15, 20)
fit_ref  <- per_rep[[1]][[FOCAL_MODEL]]
fml_ref  <- stats::formula(fit_ref)

cat("\n=== Sensitivity: number of controls (", FOCAL_MODEL,
    ", pooled over ", length(rep_files), " reps) ===\n", sep = "")

Surv   <- survival::Surv
strata <- survival::strata
set.seed(1L)

# cv_by_k[[k]][[rep]] <- list(coef, vcov)
cv_by_k <- setNames(vector("list", length(K_GRID)),
                    paste0("k", K_GRID))
for (r in seq_along(rep_files)) {
  ev_r <- read.csv(file.path(dirname(rep_files[r]),
                             "eventlist_for_rem.csv"),
                   stringsAsFactors = FALSE)
  # Same data prep as estimate_models(): alter_cons NA->0, CIRI swap,
  # polynomial terms, HR quartile dummies.
  if ("alter_cons" %in% names(ev_r))
    ev_r$alter_cons[is.na(ev_r$alter_cons)] <- 0
  pcols <- setdiff(names(ev_r), c("IS_OBSERVED", "SOURCE", "TARGET",
                                  "TIME", "ORDER", "num_possible",
                                  "num_sampled", "weight", "offset"))
  md <- ev_r[, c("TIME", "IS_OBSERVED", "ORDER", pcols,
                 intersect("offset", names(ev_r)))]
  if (all(c("ego_hr", "ego_hr_ciri") %in% names(md))) {
    md$ego_hr_fariss <- md$ego_hr
    md$ego_hr        <- md$ego_hr_ciri
  }
  for (v in c("ego_log_gdp", "ego_hr", "ego_hr_fariss",
              "alter_sponsor1year", "dyad_group_prop", "dyad_regime_prop",
              "dyad_diplo_av", "dyad_trade_gtmean", "dyad_av_diff_hr",
              "dyad_av_diff_hr_fariss", "dyad_av_diff_hr_strict",
              "dyad_av_diff_hr_fariss_strict", "alter_hr_av",
              "alter_hr_av_fariss", "alter_hr_av_strict",
              "alter_hr_av_fariss_strict"))
    if (!is.null(md[[v]])) md[[paste0(v, "_sq")]] <- md[[v]]^2
  for (v in c("ego_log_gdp", "alter_sponsor1year"))
    md[[paste0(v, "_cb")]] <- md[[v]]^3
  q_hr <- stats::quantile(md$ego_hr, c(.25, .75), na.rm = TRUE,
                          names = FALSE)
  md$ego_hr_low  <- as.integer(md$ego_hr <= q_hr[1])
  md$ego_hr_high <- as.integer(md$ego_hr >= q_hr[2])

  # Random control order within each stratum -> nested subsets.
  is_ctrl <- md$IS_OBSERVED == 0
  md$ctrl_rank <- NA_real_
  md$ctrl_rank[is_ctrl] <- ave(runif(sum(is_ctrl)), md$ORDER[is_ctrl],
                               FUN = rank)

  for (j in seq_along(K_GRID)) {
    f <- survival::coxph(fml_ref,
                         data = md[!is_ctrl | md$ctrl_rank <= K_GRID[j], ],
                         ties = "efron", robust = TRUE)
    cv_by_k[[j]][[r]] <- list(coef = coef(f), vcov = vcov(f))
  }
  cat(sprintf("  rep %d/%d fitted for all k\n", r, length(rep_files)))
}

# Pool over reps for each k
pooled_k <- lapply(cv_by_k, pool_rubin)
terms_k  <- names(pooled_k[[length(pooled_k)]]$coef)
coef_k   <- sapply(pooled_k, function(p) p$coef[terms_k])
se_k     <- sapply(pooled_k, function(p) sqrt(diag(p$vcov))[terms_k])
b_ref    <- coef_k[, ncol(coef_k)]
se_ref   <- se_k[, ncol(se_k)]

# B/W per k: does shrinking the control set amplify ordering noise?
# Report B and W separately too: the ratio can rise simply because the
# within-model SE (W) shrinks as controls are added.
bw_list <- lapply(cv_by_k, function(cvs) {
  bm <- sapply(cvs, function(x) x$coef[terms_k])
  B  <- apply(bm, 1, stats::sd)
  W  <- sqrt(diag(Reduce(`+`, lapply(cvs, function(x)
    x$vcov[terms_k, terms_k, drop = FALSE])) / length(cvs)))
  list(termwise = B / W, meanB = mean(B), meanW = mean(W))
})
bw_k     <- sapply(bw_list, `[[`, "termwise")
bw_parts <- sapply(bw_list, function(x)
  c(bw = mean(x$termwise), meanB = x$meanB, meanW = x$meanW))

# Sanity: pooled k = max controls should equal the m = 16 pooled fit
pooled_ref <- pooled_by_m[[length(pooled_by_m)]]
max_diff   <- max(abs(coef_k[, ncol(coef_k)] - pooled_ref$coef[terms_k]),
                  na.rm = TRUE)
cat(sprintf("Sanity: max |b_pooled(k=%d) - pooled m=16| = %.4g\n",
            K_GRID[length(K_GRID)], max_diff))

# Long-format table: one row per parameter x control-set size.
sens <- do.call(rbind, lapply(seq_along(K_GRID), function(j)
  data.frame(term      = terms_k,
             controls  = K_GRID[j],
             estimate  = coef_k[, j],
             se        = se_k[, j],
             D         = abs(coef_k[, j] - b_ref),
             Dstar     = abs(coef_k[, j] - b_ref) / se_ref,
             dSE       = abs(se_k[, j] - se_ref) / se_ref,
             bw        = bw_k[, j],
             row.names = NULL)))
write.csv(sens, file.path(out_dir, "controls_sensitivity.csv"),
          row.names = FALSE)
cat("Wrote", file.path(out_dir, "controls_sensitivity.csv"), "\n")

# Sign / significance stability vs the 20-control pooled estimates.
for (j in seq_along(K_GRID)[-length(K_GRID)]) {
  same_sign <- mean(sign(coef_k[, j]) == sign(b_ref))
  sig_ref   <- abs(b_ref / se_ref) > 1.96
  same_sig  <- mean((abs(coef_k[, j] / se_k[, j]) > 1.96) == sig_ref)
  kept_sig  <- if (any(sig_ref))
    mean((abs(coef_k[, j] / se_k[, j]) > 1.96)[sig_ref]) else NA_real_
  cat(sprintf(paste0("k=%2d vs k=%d: max D* = %.3f | sign match %3.0f%% | ",
                     "sig match %3.0f%% | k20-sig kept %3.0f%% | ",
                     "B/W = %.2f (mean B = %.3f, mean W = %.3f)\n"),
              K_GRID[j], K_GRID[length(K_GRID)],
              max(abs(coef_k[, j] - b_ref) / se_ref),
              100 * same_sign, 100 * same_sig, 100 * kept_sig,
              bw_parts["bw", j], bw_parts["meanB", j],
              bw_parts["meanW", j]))
}

# Small multiples: pooled coefficient + 95% CI vs number of controls.
nc <- min(5L, length(terms_k))
nr <- ceiling(length(terms_k) / nc)
pdf(file.path(out_dir, "controls_sensitivity.pdf"),
    width = 4.4 * nc, height = 3.6 * nr)
par(mfrow = c(nr, nc), mar = c(3.2, 3.6, 2.6, .6),
    oma = c(4.5, 4.5, 1, .5), mgp = c(2.4, .9, 0))
for (i in seq_along(terms_k)) {
  b  <- coef_k[i, ]
  se <- se_k[i, ]
  ylim <- range(c(b - 1.96 * se, b + 1.96 * se))
  plot(K_GRID, b, pch = 16, cex = 1.4, ylim = ylim, xaxt = "n",
       yaxt = "n", xlab = "", ylab = "", main = lab(terms_k[i]),
       cex.main = 2.5)
  axis(1, at = K_GRID, cex.axis = 1.3)
  axis(2, cex.axis = 2)
  arrows(K_GRID, b - 1.96 * se, K_GRID, b + 1.96 * se, angle = 90,
         code = 3, length = .06, lwd = 1.4)
  abline(h = 0, col = "red", lty = 2, lwd = 1.2)
  abline(h = b[length(b)], lty = 3, col = "grey50", lwd = 1.4)
}
mtext("Number of sampled controls", 1, outer = TRUE, line = 2.4, cex = 2.5)
mtext("Pooled coefficient (95% CI)", 2, outer = TRUE, line = 2, cex = 2.5)
dev.off()
cat("Wrote", file.path(out_dir, "controls_sensitivity.pdf"), "\n")
