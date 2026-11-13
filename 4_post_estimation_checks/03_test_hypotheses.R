# 03_test_hypotheses.R
#
# Hypothesis tests across the specification suite, for either model.
# Choose the model in the settings block below; everything else is shared.
#
# For each hypothesis, directional linear contrasts g'b are evaluated on
# the (pooled) coef/vcov objects with the delta method:
#   H1 (high HR more active):  d(eta)/d(HR) > 0 at HR_max (baseline margin)
#                              dummy specs: top-quartile HR dummy > 0
#   H1 (cond.):                topic delta of the marginal HR effect > 0
#                              at HR_max; dummy specs: high x topic > 0
#   H2 (low HR more active):   d(eta)/d(HR) < 0 at HR_min (baseline margin)
#                              dummy specs: bottom-quartile HR dummy > 0
#   H2 (non-cond.):            topic delta > 0 at HR_min;
#                              dummy specs: low x topic < 0
#   H3a-c (sponsoring only):   dyadic-statistic coefficient > 0
#
# Verdicts: + supported (all conditions right direction, p < .05)
#           (+) weak (right direction at p < .10, nothing wrong at .10)
#           0 null | (-) weakly rejected | - rejected (wrong direction, p < .05)
#           -- terms absent from the spec / estimates missing
#
# Reads:  sponsoring -> output/sponsoring/pooled/pooled_<model>.RDS
#         creation   -> output/creation/initiating_<model>.RDS
# Writes: output/post_estimation_checks/hypotheses_<MODEL>.csv
#         (verdict matrix) and hypotheses_<MODEL>_detail.csv (estimates,
#         z and p of every contrast)

## ---- settings -----------------------------------------------------------
MODEL <- "sponsoring"   # "sponsoring" or "creation"
ALPHA <- 0.05
## -------------------------------------------------------------------------

source(file.path(Find(function(d) file.exists(file.path(d, "paths.R")),
                      c(".", "..", "../..")), "paths.R"), chdir = TRUE)
stopifnot(MODEL %in% c("sponsoring", "creation"))

# Specification columns: label -> estimates file
specs <- if (MODEL == "sponsoring") list(
  c("Full",              "pooled_model_all_effects.RDS"),
  c("Dummies",           "pooled_model_all_effects_hr_dummies.RDS"),
  c("+ diff^2",          "pooled_model_all_effects_hrdiff_sq.RDS"),
  c("Tied (Efron)",      "pooled_model_all_effects_tied_efron.RDS"),
  c("Trade raw",         "pooled_model_all_effects_trade_raw.RDS"),
  c("Trade pos",         "pooled_model_all_effects_trade_pos.RDS"),
  c("NL effects",        "pooled_model_nl_effects.RDS"),
  c("Fariss",            "pooled_model_all_effects_hr_fariss.RDS"),
  c("Fariss dummies",    "pooled_model_all_effects_hr_fariss_dummies.RDS"),
  c("+ diff^2 Fariss",   "pooled_model_all_effects_hrdiff_sq_fariss.RDS"),
  c("NL Fariss",         "pooled_model_nl_hr_fariss.RDS")
) else list(
  c("Full",              "initiating_model_main.RDS"),
  c("Breslow",           "initiating_model_main_breslow.RDS"),
  c("Dummies",           "initiating_model_hr_dummies.RDS"),
  c("NL effects",        "initiating_model_nl_effects.RDS"),
  c("Fariss",            "initiating_model_hr_fariss.RDS"),
  c("Fariss dummies",    "initiating_model_hr_fariss_dummies.RDS"),
  c("NL Fariss",         "initiating_model_nl_hr_fariss.RDS")
)
est_dir <- if (MODEL == "sponsoring") POOLED_DIR else CREATION_OUT
spec_names <- vapply(specs, `[[`, "", 1)

# HR-score ranges for the derivative contrasts
hr_rng <- function(f) range(read.csv(f)$hr_score, na.rm = TRUE)
xrng <- list(ciri = hr_rng(HR_CIRI_CSV), fariss = hr_rng(HR_FARISS_CSV))

# ---- helpers ---------------------------------------------------------------
load_cv <- function(fname) {
  f <- file.path(est_dir, fname)
  if (!file.exists(f)) return(NULL)
  obj <- readRDS(f)
  if (inherits(obj, c("coxph", "clogit")))
    return(list(cf = stats::coef(obj), V = stats::vcov(obj)))
  list(cf = obj$coef, V = obj$vcov)
}

# HR family of a spec: "ciri" | "fariss" | "dum" | "dum_fariss"
spec_family <- function(sname) {
  far <- grepl("Fariss", sname, fixed = TRUE)
  dum <- grepl("dummies", sname, ignore.case = TRUE)
  if (far && dum) "dum_fariss" else if (far) "fariss" else if (dum) "dum" else "ciri"
}
is_dummy <- function(fam) fam %in% c("dum", "dum_fariss")
hr_dummies <- function(fam) {
  if (fam == "dum") c("ego_hr_low", "ego_hr_high")
  else c("ego_hr_fariss_low", "ego_hr_fariss_high")
}
hr_terms <- function(fam) {
  sw <- if (fam == "fariss") "ego_hr_fariss" else "ego_hr"
  c(sw, paste0(sw, "_sq"))
}
hr_int <- function(fam) {
  base <- if (is_dummy(fam)) hr_dummies(fam) else hr_terms(fam)
  paste0(base, ":alter_topic_hr")
}
slope_g <- function(x) c(1, 2 * x)   # d/dx of b1 x + b2 x^2

# Dyadic term names: strict-precedence versions in the tied model,
# alternative trade measures in the trade specs.
dyad_term <- function(sname, base) {
  if (base == "dyad_trade_gtmean" && sname == "Trade raw") return("dyad_trade_av")
  if (base == "dyad_trade_gtmean" && sname == "Trade pos") return("dyad_trade_pos")
  if (sname == "Tied (Efron)") return(paste0(base, "_strict"))
  base
}

# one contrast: list(terms, weights, expected sign)
cond_eval <- function(cv, terms, g, dir) {
  if (!all(terms %in% names(cv$cf))) return(list(available = FALSE))
  b <- cv$cf[terms]; V <- cv$V[terms, terms, drop = FALSE]
  est <- sum(g * b); se <- sqrt(drop(t(g) %*% V %*% g))
  p <- 2 * pnorm(-abs(est / se))
  list(available = TRUE, est = est, z = est / se, p = p,
       right = if (dir == ">") est > 0 else est < 0)
}

verdict <- function(conds) {
  if (any(!vapply(conds, `[[`, TRUE, "available"))) return(NA_character_)
  p <- vapply(conds, `[[`, 1, "p"); right <- vapply(conds, `[[`, TRUE, "right")
  if (any(!right & p < ALPHA)) return("-")
  if (all(right & p < ALPHA)) return("+")
  if (any(right & p < 0.10) && !any(!right & p < 0.10)) return("(+)")
  if (any(!right & p < 0.10) && !any(right & p < 0.10)) return("(-)")
  "0"
}

# ---- hypotheses ------------------------------------------------------------
sfx  <- if (MODEL == "sponsoring") "b" else "a"
verb <- if (MODEL == "sponsoring") "more active" else "creates more"

rows <- list(
  list(id = paste0("H1", sfx, " (all)"), lab = paste("High HR", verb, "(all drafts)"),
       check = function(fam, sname)
         if (is_dummy(fam)) list(list(hr_dummies(fam)[2], 1, ">"))
         else list(list(hr_terms(fam), slope_g(xrng[[fam]][2]), ">"))),
  list(id = paste0("H1", sfx, " (cond.)"), lab = "High HR on condemning drafts",
       check = function(fam, sname)
         if (is_dummy(fam)) list(list(hr_int(fam)[2], 1, ">"))
         else list(list(hr_int(fam), slope_g(xrng[[fam]][2]), ">"))),
  list(id = paste0("H2", sfx, " (all)"), lab = paste("Low HR", verb, "(all drafts)"),
       check = function(fam, sname)
         if (is_dummy(fam)) list(list(hr_dummies(fam)[1], 1, ">"))
         else list(list(hr_terms(fam), slope_g(xrng[[fam]][1]), "<"))),
  list(id = paste0("H2", sfx, " (non-cond.)"), lab = "Low HR on non-condemning drafts",
       check = function(fam, sname)
         if (is_dummy(fam)) list(list(hr_int(fam)[1], 1, "<"))
         else list(list(hr_int(fam), slope_g(xrng[[fam]][1]), ">")))
)
if (MODEL == "sponsoring") {
  dyad_row <- function(id, lab, base)
    list(id = id, lab = lab, check = function(fam, sname)
      list(list(dyad_term(sname, base), 1, ">")))
  rows <- c(rows, list(
    dyad_row("H3a", "Shared groups", "dyad_group_prop"),
    dyad_row("H3b (diplo)", "Diplomatic ties", "dyad_diplo_av"),
    dyad_row("H3b (trade)", "Trade ties", "dyad_trade_gtmean"),
    dyad_row("H3c", "Same regime", "dyad_regime_prop")))
}

# ---- run -------------------------------------------------------------------
cvs <- lapply(specs, function(s) load_cv(s[2]))
mat <- matrix("--", length(rows), length(specs),
              dimnames = list(vapply(rows, function(r) paste0(r$id, ": ", r$lab), ""),
                              spec_names))
detail <- data.frame()
for (ri in seq_along(rows)) for (si in seq_along(specs)) {
  if (is.null(cvs[[si]])) next
  sname <- spec_names[si]
  conds <- lapply(rows[[ri]]$check(spec_family(sname), sname),
                  function(c) cond_eval(cvs[[si]], c[[1]], c[[2]], c[[3]]))
  v <- verdict(conds)
  if (is.na(v)) next
  mat[ri, si] <- v
  for (c in conds) detail <- rbind(detail, data.frame(
    hypothesis = rows[[ri]]$id, spec = sname, verdict = v,
    estimate = c$est, z = c$z, p = c$p))
}

cat(sprintf("Hypothesis tests -- %s model\n\n", MODEL))
print(noquote(mat))
cat("\n+ supported (p<.05) | (+) weak (p<.10) | 0 null | (-) weakly rejected",
    "| - rejected | -- not estimable\n")

out <- file.path(POST_CHECKS_OUT, sprintf("hypotheses_%s.csv", MODEL))
write.csv(data.frame(hypothesis = rownames(mat), mat, check.names = FALSE),
          out, row.names = FALSE)
write.csv(detail, sub("\\.csv$", "_detail.csv", out), row.names = FALSE)
cat("Wrote", out, "\n")
