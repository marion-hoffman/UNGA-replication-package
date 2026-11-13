# Replication package: draft creation and joint sponsorship at the UN General Assembly

This package contains the data and R code to estimate:

- the **creation model**: a conditional (Cox) choice model of which states initiate draft resolutions;
- the **sponsoring model**: a relational-event-style conditional model of which states join a draft's sponsors. It is estimated on 16 imputed event sequences in parallel and pooled with Rubin's rules.

It also contains the pre- and post-estimation checks and the two relative-propensity figures.

## Requirements

- R (developed with R 4.4.2)
- Packages: `survival`, `data.table`, `ggplot2`, `parallel` (`parallel` ships with base R)

```r
install.packages(c("survival", "data.table", "ggplot2"))
```

## Folder structure

```
UNGA_replication/
  paths.R                      every file location, relative to this folder
  UNGA_replication.Rproj       opening it sets the working directory to this folder
  data/
    drafts/                    draft-level objects (sponsorship matrix, topics, recurrence, ...)
    countries/                 country-year attributes (GDP, conflict, HR scores, mission size, ...)
    dyads/                     country-pair attributes (shared groups, regime, diplomatic ties, trade)
  1_creation_model/            build the creation eventlist and estimate the creation models
  2_sponsoring_model/          parallel pipeline for the sponsoring models
    functions/                 eventlist construction, covariates, estimation, Rubin pooling
  3_pre_estimation_checks/     specification diagnostics on one imputed sequence
  4_post_estimation_checks/    pooling diagnostics, sensitivity to m and k, hypothesis tests
  5_figures/                   relative-propensity plots
  output/                      everything the scripts write (created on demand)
    creation/                  creation eventlist and fitted models
    sponsoring/rep_01..rep_16/ per-sequence eventlists and fits (large, see below)
    sponsoring/pooled/         Rubin-pooled estimates
    pre_estimation_checks/  post_estimation_checks/  figures/
```

## Paths

There are no absolute paths in the code. Every script starts by locating `paths.R`, either in the working directory or one or two folders above it. So a script can be run with the working directory set to the package root or to the script's own folder:

```
Rscript 1_creation_model/01_create_eventlist.R     # from the package root
cd 1_creation_model && Rscript 01_create_eventlist.R
```

In RStudio, open `UNGA_replication.Rproj` and `source()` the scripts.

## Run order

| Step | Script | Reads | Writes |
|---|---|---|---|
| 1 | `1_creation_model/01_create_eventlist.R` | `data/drafts` | `output/creation/eventlist_for_rem.csv` |
| 2 | `1_creation_model/02_add_ego_attribute_effects.R` | creation eventlist, `data/countries` | creation eventlist (adds columns) |
| 3 | `1_creation_model/03_add_ego_endo_effects.R` | creation eventlist, `data/drafts` | creation eventlist (adds columns) |
| 4 | `1_creation_model/04_estimate_models.R` | creation eventlist | `output/creation/initiating_*.RDS` |
| 5 | `2_sponsoring_model/01_master_run.R` | `data/` | `output/sponsoring/rep_XX/`, `output/sponsoring/pooled/` |
| 6 | `2_sponsoring_model/02_estimate_univariate_effects.R` | `rep_XX` eventlists | `pooled/pooled_univariate_*.RDS` |
| 7 | `3_pre_estimation_checks/run_diagnostics_rem.R` | `rep_01` eventlist | `output/pre_estimation_checks/diagnostics_rem.pdf` |
| 8 | `4_post_estimation_checks/01_run_diagnostics_pooled.R` | `rep_XX/model_result.RDS` | `output/post_estimation_checks/diagnostics_pooled.pdf` |
| 9 | `4_post_estimation_checks/02_sensitivity_imputations_controls.R` | `rep_XX` fits and eventlists | `output/post_estimation_checks/sensitivity_m_and_controls/` |
| 10 | `4_post_estimation_checks/03_test_hypotheses.R` | pooled / creation estimates | `output/post_estimation_checks/hypotheses_<model>.csv` |
| 11 | `5_figures/plot_relative_propensity_hr.R` | pooled estimates, creation eventlist, `rep_01` eventlist | `output/figures/relative_propensity_hr.{pdf,png}` |
| 12 | `5_figures/plot_relative_propensity_dyad.R` | pooled estimates, `rep_01` eventlist | `output/figures/relative_propensity_dyad.{pdf,png}` |

Notes:

- **Step 5** runs the full sponsoring pipeline for each of the 16 sequences, one PSOCK worker per sequence (`min(detectCores(), 16)` workers). Each sequence draws its own random ordering of simultaneous joins and its own 20 sampled non-events per event (seed `1000 + rep`). This step is the expensive one: it takes several hours and a lot of memory. It then pools every specification in `functions/functions_estimate_models.R` with Rubin's rules. WARNING: this takes a long time!
- **Step 10**: set `MODEL <- "sponsoring"` or `MODEL <- "creation"` at the top of the script.
- **Step 9** refits the main model 64 times (16 sequences × 4 control-set sizes).

## Shipped results and the per-sequence folders

`output/` already contains the results reported in the paper:

- `output/creation/`: the creation eventlist and fitted models (small);
- `output/sponsoring/pooled/`: the pooled sponsoring estimates (small);
- `output/sponsoring/rep_01` … `rep_16`: the per-sequence eventlists and fitted models (about 9 GB in total).

With these files, steps 7–12 can be run without re-estimating anything.

The `rep_XX` folders are not included in the GitHub repository, because their files exceed GitHub's 100 MB file-size limit. Run step 5 to create them.

The `rep_XX` folders can be deleted to make the package smaller. Steps 10 and 11 only need the pooled and creation estimates: without `rep_01`, the HR plot takes its reference values from the creation eventlist. Steps 6–9 and 12 need the `rep_XX` folders, which are rebuilt by running step 5.

## Data

| File | Content |
|---|---|
| `drafts/country_draft_abs.RDS` | draft × state sponsorship matrix (join order; 1 = creator), with draft dates and identifiers |
| `drafts/drafts_2017_2020.RDS` | drafts included in the analysis |
| `drafts/drafts_topic_condemn_s70_s74_new.rds` | draft topics (HR condemnation, cooperation, Israel–Palestine) |
| `drafts/drafts_title_recurrence_s70_s74.rds` | recurrence of draft titles |
| `drafts/recur_covariate_consensus.RDS` | adoption by consensus in the previous cycle |
| `drafts/recur_covariate_spons_last.RDS` | sponsorship of the draft in the previous cycle |
| `drafts/drafts_nrefs_s70on.rds` | number of references in the draft |
| `drafts/country_condemnation_variables_final.csv` | states condemned by the GA (creation model) |
| `countries/gdp_2015_2020.csv` | GDP |
| `countries/civil_confl.csv` | civil conflict |
| `countries/hr_ciri_2015_2020.csv` | CIRI human-rights score (main measure) |
| `countries/hr_fariss_2015_2020.csv` | Fariss latent human-rights score (robustness) |
| `countries/mission_size_2015_2020.csv` | size of the permanent mission |
| `countries/ego_gr_2015_2020.csv` | number of group memberships |
| `countries/ego_regime_2015_2020.csv` | regime type |
| `dyads/same_gr_2015_2020.csv` | shared group membership |
| `dyads/same_regime_2015_2020.csv` | same regime type |
| `dyads/diplo_2015_2020.csv` | diplomatic ties |
| `dyads/red_2015_2020.csv` | trade dependence |
