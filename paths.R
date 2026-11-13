# paths.R -- every file location used by the replication scripts, relative
# to the package root (the folder containing this file).
#
# Each script starts with
#   source(file.path(Find(function(d) file.exists(file.path(d, "paths.R")),
#                         c(".", "..", "../..")), "paths.R"), chdir = TRUE)
# so it can be run with the working directory set either to the package
# root or to the script's own folder, e.g.
#   Rscript 1_creation_model/01_create_eventlist.R      (from the root)
#   Rscript 01_create_eventlist.R                       (from 1_creation_model/)
# Nothing in this file depends on where the package is stored.

ROOT <- normalizePath(getwd(), winslash = "/")

# ---- input data ---------------------------------------------------------
DATA_DIR      <- file.path(ROOT, "data")
DRAFTS_DIR    <- file.path(DATA_DIR, "drafts")     # draft-level objects
COUNTRIES_DIR <- file.path(DATA_DIR, "countries")  # country-year attributes
DYADS_DIR     <- file.path(DATA_DIR, "dyads")      # country-pair attributes

# drafts
COUNTRY_DRAFT_ABS  <- file.path(DRAFTS_DIR, "country_draft_abs.RDS")
DRAFTS_KEEP        <- file.path(DRAFTS_DIR, "drafts_2017_2020.RDS")
TOPIC_RDS          <- file.path(DRAFTS_DIR, "drafts_topic_condemn_s70_s74_new.rds")
RECUR_TITLE_RDS    <- file.path(DRAFTS_DIR, "drafts_title_recurrence_s70_s74.rds")
CONSENSUS_RDS      <- file.path(DRAFTS_DIR, "recur_covariate_consensus.RDS")
NREFS_RDS          <- file.path(DRAFTS_DIR, "drafts_nrefs_s70on.rds")
RECUR_DYAD_RDS     <- file.path(DRAFTS_DIR, "recur_covariate_spons_last.RDS")
CONDEMNATION_CSV   <- file.path(DRAFTS_DIR, "country_condemnation_variables_final.csv")

# countries
GDP_CSV      <- file.path(COUNTRIES_DIR, "gdp_2015_2020.csv")
CONFLICT_CSV <- file.path(COUNTRIES_DIR, "civil_confl.csv")
HR_FARISS_CSV <- file.path(COUNTRIES_DIR, "hr_fariss_2015_2020.csv")
HR_CIRI_CSV  <- file.path(COUNTRIES_DIR, "hr_ciri_2015_2020.csv")
STAFF_CSV    <- file.path(COUNTRIES_DIR, "mission_size_2015_2020.csv")
GROUP_CSV    <- file.path(COUNTRIES_DIR, "ego_gr_2015_2020.csv")
REGIME_CSV   <- file.path(COUNTRIES_DIR, "ego_regime_2015_2020.csv")

# dyads
SAME_GROUP_CSV  <- file.path(DYADS_DIR, "same_gr_2015_2020.csv")
SAME_REGIME_CSV <- file.path(DYADS_DIR, "same_regime_2015_2020.csv")
DIPLO_CSV       <- file.path(DYADS_DIR, "diplo_2015_2020.csv")
TRADE_CSV       <- file.path(DYADS_DIR, "red_2015_2020.csv")

# ---- code ---------------------------------------------------------------
SPONSORING_FUNCTIONS_DIR <- file.path(ROOT, "2_sponsoring_model", "functions")

# ---- outputs (created on demand) ----------------------------------------
OUTPUT_DIR <- file.path(ROOT, "output")
# creation model: eventlist + fitted models
CREATION_OUT <- file.path(OUTPUT_DIR, "creation")
CREATION_EVENTLIST <- file.path(CREATION_OUT, "eventlist_for_rem.csv")
# sponsoring model: one folder per imputed sequence (rep_XX/, large and
# safe to delete -- 2_sponsoring_model/01_master_run.R rebuilds them) and
# the Rubin-pooled estimates (pooled/, small)
SPONSORING_OUT <- file.path(OUTPUT_DIR, "sponsoring")
POOLED_DIR     <- file.path(SPONSORING_OUT, "pooled")
REP01_EVENTLIST <- file.path(SPONSORING_OUT, "rep_01", "eventlist_for_rem.csv")
# checks and figures
PRE_CHECKS_OUT  <- file.path(OUTPUT_DIR, "pre_estimation_checks")
POST_CHECKS_OUT <- file.path(OUTPUT_DIR, "post_estimation_checks")
FIGURES_OUT     <- file.path(OUTPUT_DIR, "figures")

for (.d in c(CREATION_OUT, SPONSORING_OUT, PRE_CHECKS_OUT, POST_CHECKS_OUT,
             FIGURES_OUT))
  dir.create(.d, recursive = TRUE, showWarnings = FALSE)
rm(.d)
