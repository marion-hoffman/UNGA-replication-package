# Diagnostics for the pooled parallel-pipeline results.
# Run after 2_sponsoring_model/01_master_run.R has produced
# output/sponsoring/rep_XX/model_result.RDS (and output/sponsoring/pooled/).
# Answers: is the Rubin's-rules pooling stable? i.e. how much of the pooled
# uncertainty comes from the random control draw, are sign/significance
# conclusions consistent across reps, do the pooled estimates converge in m,
# and is any single rep driving the result.

source(file.path(Find(function(d) file.exists(file.path(d, "paths.R")),
                      c(".", "..", "../..")), "paths.R"), chdir = TRUE)
source(file.path(SPONSORING_FUNCTIONS_DIR, "functions_pool_rubin.R"))

results_dir <- SPONSORING_OUT
stopifnot(dir.exists(results_dir))
main_model <- "model_all_effects"

## ===== 0. Load per-rep fits =====
rep_dirs <- list.dirs(results_dir, full.names = FALSE, recursive = FALSE)
rep_dirs <- sort(rep_dirs[grepl("^rep_", rep_dirs)])
stopifnot(length(rep_dirs) > 0)

fits_by_model <- list()   # model -> named list rep -> fit object
for (rd in rep_dirs) {
  rp <- file.path(results_dir, rd, "model_result.RDS")
  if (!file.exists(rp)) next
  res <- readRDS(rp)
  for (m in names(res)) fits_by_model[[m]][[rd]] <- res[[m]]
}
cat(sprintf("Loaded %d rep(s) [%s]; %d models\n",
            length(rep_dirs), paste(rep_dirs, collapse = ","),
            length(fits_by_model)))

# Extract {coef, vcov} lists per model (valid fits only)
as_cv <- function(fits) {
  lapply(fits, function(f) {
    cf <- tryCatch(coef(f), error = function(e) NULL)
    vc <- tryCatch(vcov(f), error = function(e) NULL)
    if (is.null(cf) || is.null(vc)) return(NULL)
    list(coef = cf, vcov = vc)
  })
}
valid_fits <- function(mname) Filter(function(f)
  !is.null(f) && inherits(f, c("clogit", "coxph")), fits_by_model[[mname]])

## ===== 1. Coverage: which models fitted in which reps =====
cat("\n===== 1. Fit coverage per model =====\n")
cov_tab <- sapply(fits_by_model, function(l)
  vapply(rep_dirs, function(rd)
    !is.null(l[[rd]]) && inherits(l[[rd]], c("clogit", "coxph")),
    logical(1)))
print(cov_tab)
cat("Successful fits per model:\n"); print(colSums(cov_tab))
n_main <- sum(cov_tab[, main_model])
if (n_main < 3) stop("Too few valid ", main_model, " fits to diagnose pooling")

## ===== 2. Variance decomposition: within vs between (FMI) =====
cat("\n===== 2. Within/between decomposition:", main_model, "=====\n")
rl <- as_cv(valid_fits(main_model))
rl <- Filter(Negate(is.null), rl)
m <- length(rl)
params <- Reduce(intersect, lapply(rl, function(r) names(r$coef)))
betas <- do.call(cbind, lapply(rl, function(r) r$coef[params]))          # p x m
ses   <- do.call(cbind, lapply(rl, function(r)
  setNames(sqrt(diag(r$vcov)), rownames(r$vcov))[params]))                # p x m
rownames(betas) <- rownames(ses) <- params

beta_bar <- rowMeans(betas)
B  <- apply(betas, 1, var)              # between-rep variance
W  <- rowMeans(ses^2)                   # within-rep variance
Tv <- W + (1 + 1/m) * B                 # Rubin total variance
se_pooled <- sqrt(Tv)
fmi <- (1 + 1/m) * B / Tv               # fraction of missing information

tab2 <- data.frame(
  pooled    = beta_bar,
  se_pooled = se_pooled,
  sd_reps   = sqrt(B),
  se_within = sqrt(W),
  fmi       = round(fmi, 3),
  sign_agree = rowMeans(sign(betas) == sign(beta_bar)),
  frac_sig   = rowMeans(abs(betas / ses) > 1.96))
rownames(tab2) <- params
print(round(tab2[order(-fmi), ], 3))
cat("\nFMI = share of pooled variance due to the random control draw.\n")
cat("FMI > ~0.3 for a reported coefficient = conclusions may be\n")
cat("draw-dependent; sign_agree < 1 means reps disagree on direction.\n")

## ===== 3. m-sensitivity: do pooled estimates converge in m? =====
cat("\n===== 3. m-sensitivity of pooled estimates =====\n")
# Re-pool on R random subsets of size k; MC sd of the pooled estimate
# should fall ~1/sqrt(k). If estimates still move a lot at k ~ m, more
# reps would change the answer.
ks <- sort(unique(pmin(m, c(4, 8, 12, m))))
R  <- 50
mc_sd <- matrix(NA_real_, nrow = length(params), ncol = length(ks),
                dimnames = list(params, paste0("k=", ks)))
set.seed(1)
for (j in seq_along(ks)) {
  k <- ks[j]
  ests <- replicate(R, {
    idx <- sample(m, k)
    pool_rubin(rl[idx])$coef[params]
  })
  if (is.null(dim(ests))) ests <- matrix(ests, ncol = R, dimnames = list(params, NULL))
  mc_sd[, j] <- apply(ests, 1, sd, na.rm = TRUE)
}
cat("MC sd of pooled estimate vs subset size (cols = k, out of m =", m, "):\n")
print(round(mc_sd, 4))
cat("Expected scaling sd(k) ~ sd_reps / sqrt(k); persistently large\n")
cat("values at k = m signal the pooled estimate is not yet converged.\n")

## ===== 4. Leave-one-rep-out influence =====
cat("\n===== 4. Leave-one-rep-out =====\n")
loo <- do.call(cbind, lapply(seq_len(m), function(j)
  pool_rubin(rl[-j])$coef[params]))
delta <- abs(loo - beta_bar) / se_pooled   # shift in units of pooled SE
worst_rep <- apply(delta, 1, which.max)
loo_tab <- data.frame(
  max_shift_in_SE = round(apply(delta, 1, max), 3),
  rep_dropped     = names(rl)[worst_rep])
rownames(loo_tab) <- params
print(loo_tab[order(-loo_tab$max_shift_in_SE), , drop = FALSE])
cat("Shifts are in units of the pooled SE; > ~0.1 means one rep\n")
cat("materially moves that coefficient.\n")

## ===== 5. Robustness specs: consistency across reps =====
cat("\n===== 5. Robustness spec consistency =====\n")
# Per-rep paired dAIC vs the main model (clogit alternates only; tied
# coxph specs use merged risk sets and a different likelihood scale), and
# pooled-coefficient agreement on shared terms.
aic_of <- function(f) tryCatch(AIC(f), error = function(e) NA_real_)
alt_models <- setdiff(names(fits_by_model),
                      c(main_model, grep("tied", names(fits_by_model), value = TRUE)))
fits_main <- valid_fits(main_model)
for (alt in alt_models) {
  fits_alt <- valid_fits(alt)
  common <- intersect(names(fits_main), names(fits_alt))
  if (length(common) == 0) next
  is_clogit <- inherits(fits_alt[[common[1]]], "clogit")
  if (is_clogit) {
    d <- vapply(common, function(rd)
      aic_of(fits_alt[[rd]]) - aic_of(fits_main[[rd]]), numeric(1))
    cat(sprintf("\n-- %s (%d reps)\n", alt, length(common)))
    cat(sprintf("   dAIC vs %s: mean %+.1f, sd %.1f, reps>0: %.0f%%\n",
                main_model, mean(d, na.rm = TRUE), sd(d, na.rm = TRUE),
                100 * mean(d > 0, na.rm = TRUE)))
  } else cat(sprintf("\n-- %s (coxph, AIC not comparable)\n", alt))
  # Pooled coefficient comparison on shared terms
  pl_alt  <- pool_rubin(as_cv(fits_alt))
  shared_p <- intersect(names(beta_bar), names(pl_alt$coef))
  b1 <- beta_bar[shared_p]; b2 <- pl_alt$coef[shared_p]
  rel <- abs(b2 - b1) / pmax(abs(b1), 1e-9)
  flips <- shared_p[sign(b1) != sign(b2)]
  cat(sprintf("   shared terms: median |rel change| %.2f, max %.2f (%s)",
              median(rel), max(rel),
              if (length(flips)) paste("SIGN FLIPS:", paste(flips, collapse = ", "))
              else "no sign flips"))
  cat("\n")
}

## ===== 6. Plots =====
pdf(file.path(POST_CHECKS_OUT, "diagnostics_pooled.pdf"), width = 9, height = 7)

# 6a. Caterpillar: per-rep estimate + 95% CI per term of the main model,
#     sorted by estimate; horizontal line = pooled estimate.
pp <- 12; np <- ceiling(length(params) / pp)
for (pg in seq_len(np)) {
  idx <- ((pg - 1) * pp + 1):min(pg * pp, length(params))
  nc <- 4; nr <- ceiling(length(idx) / nc)
  par(mfrow = c(nr, nc), mar = c(3, 3, 2.2, .5), mgp = c(1.8, .6, 0))
  for (i in idx) {
    o <- order(betas[i, ])
    yl <- range(c(betas[i, o] - 1.96 * ses[i, o],
                  betas[i, o] + 1.96 * ses[i, o]), finite = TRUE)
    plot(seq_len(m), betas[i, o], pch = 19, cex = .6, ylim = yl,
         xlab = "rep (sorted)", ylab = "coef",
         main = paste0(params[i], "  [FMI=", round(fmi[i], 2), "]"),
         cex.main = .75, xaxt = "n")
    arrows(seq_len(m), betas[i, o] - 1.96 * ses[i, o],
           seq_len(m), betas[i, o] + 1.96 * ses[i, o],
           length = .02, angle = 90, code = 3, col = "grey60")
    abline(h = beta_bar[i], col = "red", lwd = 1.5)
    abline(h = beta_bar[i] + c(-1.96, 1.96) * se_pooled[i],
           col = "red", lty = 3)
    abline(h = 0, col = "grey80")
  }
}

# 6b. m-sensitivity: MC sd of pooled estimate vs k, one line per term,
#     normalized by the between-rep sd (expect ~1/sqrt(k) decay).
par(mfrow = c(1, 1), mar = c(4, 4, 2.5, 1))
norm <- sweep(mc_sd, 1, sqrt(B), `/`)
matplot(ks, t(norm), type = "b", pch = 19, cex = .5, lty = 1,
        col = adjustcolor("grey50", .6), xlab = "subset size k",
        ylab = "MC sd / between-rep sd",
        main = "Pooling convergence in m (grey = one line per term)")
lines(ks, 1/sqrt(ks), col = "red", lwd = 2)
legend("topright", legend = "1/sqrt(k) benchmark", col = "red", lwd = 2,
       bty = "n", cex = .8)

# 6c. Per-rep paired dAIC for each clogit robustness spec
alt_cl <- Filter(function(a) {
  f <- valid_fits(a)
  length(f) > 0 && inherits(f[[1]], "clogit")
}, alt_models)
if (length(alt_cl)) {
  par(mfrow = c(1, length(alt_cl)), mar = c(4, 4, 2.5, .5))
  for (alt in alt_cl) {
    fits_alt <- valid_fits(alt)
    common <- intersect(names(fits_main), names(fits_alt))
    d <- vapply(common, function(rd)
      aic_of(fits_alt[[rd]]) - aic_of(fits_main[[rd]]), numeric(1))
    barplot(d, names.arg = gsub("rep_", "", common),
            ylab = paste("dAIC vs", main_model),
            main = gsub("model_all_effects_", "", alt), cex.main = .8,
            las = 2, cex.names = .7)
    abline(h = 0, col = "red")
  }
}
dev.off()
cat("\nDone. Plots:", file.path(POST_CHECKS_OUT, "diagnostics_pooled.pdf"), "\n")
