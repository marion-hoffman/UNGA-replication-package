add_dyad_attribute_effects <- function(
  all_events_csv,
  group_csv,
  regime_csv,
  diplo_csv,
  trade_csv,
  hr_ciri_csv,
  hr_fariss_csv,
  recur_dyad_rds,
  country_draft_abs_path,
  out_csv_path = all_events_csv
) {
  data.table::setDTthreads(1)
  dt <- data.table::fread(all_events_csv)
  group_file <- read.csv(group_csv, stringsAsFactors = FALSE)
  regime_file <- read.csv(regime_csv, stringsAsFactors = FALSE)
  diplo_file <- read.csv(diplo_csv, stringsAsFactors = FALSE)
  trade_file <- read.csv(trade_csv, stringsAsFactors = FALSE)
  hr_file <- read.csv(hr_ciri_csv, stringsAsFactors = FALSE)
  hr_fariss_file <- read.csv(hr_fariss_csv, stringsAsFactors = FALSE)

  # Trade enters as share-of-previous-sponsors covariates:
  #   dyad_trade_gtmean = proportion of previous sponsors whose dependence
  #     on the ego exceeds this year's matrix mean (main spec);
  #   dyad_trade_prop   = share with nonzero trade;
  #   dyad_trade_av     = raw mean importance (robustness spec);
  #   dyad_trade_pos    = proportion with strictly positive trade
  #     (robustness spec).

  # collect all state ISO codes
  states <- unique(c(group_file$iso1, group_file$iso2,
                     regime_file$iso1, regime_file$iso2,
                     diplo_file$iso_o, diplo_file$iso_d,
                     trade_file$iso3_ego, trade_file$iso3_alter))
  states <- states[!is.na(states) & states != ""]

  years <- sort(unique(c(group_file$year, regime_file$year,
                         diplo_file$Year, trade_file$year)))

  # create yearly adjacency matrices (vectorized fills instead of row loops)
  adj_group   <- list(); adj_regime <- list()
  adj_diplo   <- list(); adj_trade  <- list()
  trade_mean  <- numeric(0)   # mean trade dependence per year (all dyad cells)

  for (yr in years) {
    mat_g <- matrix(0, length(states), length(states), dimnames=list(states, states))
    sub_g <- group_file[group_file$year == yr, ]
    mat_g[cbind(sub_g$iso1, sub_g$iso2)] <- sub_g$shared_groups_any
    mat_g[cbind(sub_g$iso2, sub_g$iso1)] <- sub_g$shared_groups_any
    adj_group[[as.character(yr)]] <- mat_g

    mat_r <- matrix(0, length(states), length(states), dimnames=list(states, states))
    sub_r <- regime_file[regime_file$year == yr, ]
    mat_r[cbind(sub_r$iso1, sub_r$iso2)] <- sub_r$same_regime
    mat_r[cbind(sub_r$iso2, sub_r$iso1)] <- sub_r$same_regime
    adj_regime[[as.character(yr)]] <- mat_r

    mat_d <- matrix(0, length(states), length(states), dimnames=list(states, states))
    sub_d <- diplo_file[diplo_file$Year == yr, ]
    # LOR is directional (iso_o = representation host -> iso_d = receiver);
    # do NOT symmetrize -- an embassy a->b does not imply b->a.
    mat_d[cbind(sub_d$iso_o, sub_d$iso_d)] <- sub_d$LOR
    adj_diplo[[as.character(yr)]] <- mat_d

    mat_t <- matrix(0, length(states), length(states), dimnames=list(states, states))
    sub_t <- trade_file[trade_file$year == yr, ]
    mat_t[cbind(sub_t$iso3_ego, sub_t$iso3_alter)] <- sub_t$RED_export_importance
    adj_trade[[as.character(yr)]] <- mat_t
    trade_mean[as.character(yr)] <- mean(mat_t)
  }

  # HR-score lookups for dyad_av_diff_hr(+_fariss): iso x year matrices
  # on the CIRI (main) and Fariss scales. Missing country-years get the
  # year's median (same imputation as the ego attribute), remaining gaps
  # the global median.
  build_hr_mat <- function(file) {
    hf <- file[!is.na(file$iso) & file$iso %in% states, ]
    yrs <- sort(unique(hf$year))
    m <- matrix(NA_real_, nrow = length(states), ncol = length(yrs),
                dimnames = list(states, as.character(yrs)))
    m[cbind(hf$iso, as.character(hf$year))] <- hf$hr_score
    for (yr in yrs) {
      yc <- as.character(yr)
      m[is.na(m[, yc]), yc] <- median(m[, yc], na.rm = TRUE)
    }
    m[is.na(m)] <- median(m, na.rm = TRUE)
    list(mat = m, years = yrs,
         year = function(y) yrs[which.min(abs(yrs - y))])
  }
  hr  <- build_hr_mat(hr_file)
  hrf <- build_hr_mat(hr_fariss_file)

  event_year <- function(t) {
    y <- as.integer(format(as.POSIXct(t, origin="1970-01-01", tz="UTC"), "%Y"))
    nearest <- years[which.min(abs(years - y))]
    nearest
  }

  # Previous senders are per ROW (built from the row's own TARGET draft);
  # the year matrices and time window are per-EVENT, so loop over unique
  # ORDER values for the shared parts and rows for the prev-sender sets.
  prev_rds <- file.path(dirname(all_events_csv), "prev_senders_list.RDS")
  strict_rds <- file.path(dirname(all_events_csv), "prev_senders_strict_list.RDS")
  prev_senders_list <- readRDS(prev_rds)
  if (file.exists(strict_rds)) {
    prev_senders_strict_list <- readRDS(strict_rds)
  } else {
    stop("prev_senders_strict_list.RDS not found next to prev_senders_list.RDS -- ",
         "rerun precompute_prev_senders() (it now emits both files).")
  }
  n <- nrow(dt)
  if (length(prev_senders_list) != n || length(prev_senders_strict_list) != n)
    stop("prev_senders lists do not match eventlist rows (", n,
         ") -- stale RDS? rerun precompute_prev_senders().")

  # per-row accumulation vectors (all 0 for events with no previous senders)
  v_group_prop <- numeric(n); v_group_all <- integer(n)
  v_regime_prop <- numeric(n); v_regime_all <- integer(n)
  v_diplo_av <- numeric(n)
  v_trade_prop <- numeric(n); v_trade_av <- numeric(n)
  v_trade_gtmean <- numeric(n); v_trade_pos <- numeric(n)
  v_group_prop_s <- numeric(n); v_group_all_s <- integer(n)
  v_regime_prop_s <- numeric(n); v_regime_all_s <- integer(n)
  v_diplo_av_s <- numeric(n)
  v_trade_prop_s <- numeric(n); v_trade_av_s <- numeric(n)
  v_trade_gtmean_s <- numeric(n); v_trade_pos_s <- numeric(n)
  v_av_diff_hr <- numeric(n); v_av_diff_hr_s <- numeric(n)
  v_av_diff_hr_f <- numeric(n); v_av_diff_hr_fs <- numeric(n)
  v_alter_hr_av <- numeric(n); v_alter_hr_av_s <- numeric(n)
  v_alter_hr_av_f <- numeric(n); v_alter_hr_av_fs <- numeric(n)

  # Matrices are fetched once per ORDER (one event year), but the
  # previous-sender sets are per ROW: each row's covariates use the
  # senders of ITS OWN TARGET draft.
  ord_groups <- split(seq_len(n), dt$ORDER)
  ng <- length(ord_groups)
  for (g in seq_along(ord_groups)) {
    if (g %% max(1, floor(ng/20)) == 0) cat(sprintf("%d%% done\n", round(100*g/ng)))
    rows <- ord_groups[[g]]
    year <- event_year(dt$TIME[rows[1]])
    group_mat <- adj_group[[as.character(year)]]
    regime_mat <- adj_regime[[as.character(year)]]
    diplo_mat <- adj_diplo[[as.character(year)]]
    trade_mat <- adj_trade[[as.character(year)]]
    hr_vec  <- hr$mat[,  as.character(hr$year(year))]
    hrf_vec <- hrf$mat[, as.character(hrf$year(year))]
    # year median for scores missing from the matrix (unlisted ISOs)
    hr_med  <- median(hr_vec,  na.rm = TRUE)
    hrf_med <- median(hrf_vec, na.rm = TRUE)

    for (i in rows) {
      sender <- dt$SOURCE[i]
      prev_senders <- prev_senders_list[[i]]
      if (length(prev_senders) == 0) next

      # GROUP
      group_vals <- group_mat[sender, prev_senders]
      v_group_prop[i] <- mean(group_vals)
      v_group_all[i] <- as.integer(all(group_vals == 1))
      # REGIME
      regime_vals <- regime_mat[sender, prev_senders]
      v_regime_prop[i] <- mean(regime_vals)
      v_regime_all[i] <- as.integer(all(regime_vals == 1))
      # DIPLO
      v_diplo_av[i] <- mean(diplo_mat[sender, prev_senders])
      # TRADE: prop = share of previous sponsors with nonzero trade (main
      # spec); av = raw mean importance over all previous sponsors
      # (robustness spec)
      trade_vals <- trade_mat[sender, prev_senders]
      v_trade_prop[i] <- mean(trade_vals > 0)
      v_trade_av[i] <- mean(trade_vals)
      # Alternative trade summaries: proportion above this year's matrix
      # mean (main spec) and proportion strictly positive (robustness).
      v_trade_gtmean[i] <- mean(trade_vals > trade_mean[as.character(year)])
      v_trade_pos[i] <- mean(trade_vals > 0)
      # AV |HR DIFF|: mean absolute HR-score difference between the ego
      # and each previous sponsor (homophily robustness check) -- CIRI
      # main measure, Fariss alternative. Scores missing from the
      # matrix (unlisted ISOs) are imputed with the year median.
      ego_v <- hr_vec[sender]; if (is.na(ego_v)) ego_v <- hr_med
      alt_v <- hr_vec[prev_senders]; alt_v[is.na(alt_v)] <- hr_med
      v_av_diff_hr[i] <- mean(abs(ego_v - alt_v))
      v_alter_hr_av[i] <- mean(alt_v)
      ego_vf <- hrf_vec[sender]; if (is.na(ego_vf)) ego_vf <- hrf_med
      alt_vf <- hrf_vec[prev_senders]; alt_vf[is.na(alt_vf)] <- hrf_med
      v_av_diff_hr_f[i] <- mean(abs(ego_vf - alt_vf))
      v_alter_hr_av_f[i] <- mean(alt_vf)

      # STRICT PRECEDENCE: only strictly preceding sponsors (TIME < time);
      # computed only when the row has any previous senders (as before).
      prev_senders_strict <- prev_senders_strict_list[[i]]
      if (length(prev_senders_strict) > 0) {
        group_vals_s <- group_mat[sender, prev_senders_strict]
        v_group_prop_s[i] <- mean(group_vals_s)
        v_group_all_s[i] <- as.integer(all(group_vals_s == 1))
        regime_vals_s <- regime_mat[sender, prev_senders_strict]
        v_regime_prop_s[i] <- mean(regime_vals_s)
        v_regime_all_s[i] <- as.integer(all(regime_vals_s == 1))
        v_diplo_av_s[i] <- mean(diplo_mat[sender, prev_senders_strict])
        trade_vals_s <- trade_mat[sender, prev_senders_strict]
        v_trade_prop_s[i] <- mean(trade_vals_s > 0)
        v_trade_av_s[i] <- mean(trade_vals_s)
        v_trade_gtmean_s[i] <- mean(trade_vals_s > trade_mean[as.character(year)])
        v_trade_pos_s[i] <- mean(trade_vals_s > 0)
        ego_vs <- hr_vec[sender]; if (is.na(ego_vs)) ego_vs <- hr_med
        alt_vs <- hr_vec[prev_senders_strict]
        alt_vs[is.na(alt_vs)] <- hr_med
        v_av_diff_hr_s[i] <- mean(abs(ego_vs - alt_vs))
        v_alter_hr_av_s[i] <- mean(alt_vs)
        ego_vfs <- hrf_vec[sender]; if (is.na(ego_vfs)) ego_vfs <- hrf_med
        alt_vfs <- hrf_vec[prev_senders_strict]
        alt_vfs[is.na(alt_vfs)] <- hrf_med
        v_av_diff_hr_fs[i] <- mean(abs(ego_vfs - alt_vfs))
        v_alter_hr_av_fs[i] <- mean(alt_vfs)
      }
    }
  }

  dt[, `:=`(
    dyad_group_prop = v_group_prop,
    dyad_group_all = v_group_all,
    dyad_regime_prop = v_regime_prop,
    dyad_regime_all = v_regime_all,
    dyad_diplo_av = v_diplo_av,
    dyad_trade_prop = v_trade_prop,
    dyad_trade_av = v_trade_av,
    dyad_trade_gtmean = v_trade_gtmean,
    dyad_trade_pos = v_trade_pos,
    dyad_group_prop_strict = v_group_prop_s,
    dyad_group_all_strict = v_group_all_s,
    dyad_regime_prop_strict = v_regime_prop_s,
    dyad_regime_all_strict = v_regime_all_s,
    dyad_diplo_av_strict = v_diplo_av_s,
    dyad_trade_prop_strict = v_trade_prop_s,
    dyad_trade_av_strict = v_trade_av_s,
    dyad_trade_gtmean_strict = v_trade_gtmean_s,
    dyad_trade_pos_strict = v_trade_pos_s,
    dyad_av_diff_hr = v_av_diff_hr,
    dyad_av_diff_hr_strict = v_av_diff_hr_s,
    dyad_av_diff_hr_fariss = v_av_diff_hr_f,
    dyad_av_diff_hr_fariss_strict = v_av_diff_hr_fs,
    alter_hr_av = v_alter_hr_av,
    alter_hr_av_strict = v_alter_hr_av_s,
    alter_hr_av_fariss = v_alter_hr_av_f,
    alter_hr_av_fariss_strict = v_alter_hr_av_fs,
    dyad_recur = 0
  )]

  # ---- dyad_recur: vectorized lookup into the draft x sender matrix ----
  df_recur_dyad <- readRDS(recur_dyad_rds)
  cat("Reading", recur_dyad_rds, "\n")
  recur_id <- as.character(df_recur_dyad$ldoc_final)
  rmat <- as.matrix(df_recur_dyad[, setdiff(names(df_recur_dyad), "ldoc_final"),
                                drop = FALSE])
  ri <- match(as.character(dt$TARGET), recur_id)
  ci <- match(as.character(dt$SOURCE), colnames(rmat))
  val <- rmat[cbind(ri, ci)]
  dt[, dyad_recur := as.integer(!is.na(val) & val == 1)]

  data.table::fwrite(dt, out_csv_path)
  invisible(out_csv_path)
}
