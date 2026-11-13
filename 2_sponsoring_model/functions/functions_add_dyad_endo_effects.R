add_dyad_endo_effects <- function(
  all_events_csv,
  country_draft_abs_path,
  prev_senders_rds = file.path(dirname(all_events_csv), "prev_senders_list.RDS"),
  out_csv_path = all_events_csv
) {
  data.table::setDTthreads(1)
  dt <- data.table::fread(all_events_csv)
  country_draft_abs <- readRDS(country_draft_abs_path)
  data.table::setDT(country_draft_abs)
  date_or_posix <- as.POSIXct(country_draft_abs$date_or, tz="UTC")
  country_draft_abs[, date_or_numeric := as.numeric(date_or_posix)]

  prev_senders_list <- readRDS(prev_senders_rds)
  strict_rds <- file.path(dirname(all_events_csv), "prev_senders_strict_list.RDS")
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

  # Draft x state indicator matrix: D[d, s] = 1 iff state s sponsored draft d
  # (value present and > 0). crossprod of the windowed rows then gives, for
  # every sender pair, the count of drafts in the window they co-sponsored --
  # one matrix multiply per event instead of a data.table filter per pair.
  state_cols <- names(country_draft_abs)[30:225]
  D <- as.matrix(country_draft_abs[, ..state_cols])
  D <- (!is.na(D) & D > 0) * 1
  dimnames(D) <- list(as.character(country_draft_abs$ldoc_final), state_cols)
  date_or_num <- country_draft_abs$date_or_numeric

  v_4cycle <- numeric(n); v_4cycle_s <- numeric(n)

  # The 1-year draft window is per EVENT (ORDER), so the co-sponsorship
  # matrix is computed once per event over the union of all states the
  # group's rows need; each row's own prev-sender set (keyed by its
  # TARGET) then picks its mean.
  ord_groups <- split(seq_len(n), dt$ORDER)
  ng <- length(ord_groups)
  for (g in seq_along(ord_groups)) {
    if (g %% max(1, floor(ng/20)) == 0) cat(sprintf("%d%% done\n", round(100*g/ng)))
    rows <- ord_groups[[g]]

    time_num <- as.numeric(dt$TIME[rows[1]])
    time_window_start <- time_num - 365*24*60*60
    in_window <- which(date_or_num >= time_window_start & date_or_num < time_num)

    # union of every state any row in this event needs as a matrix index
    rel <- unique(c(dt$SOURCE[rows],
                    unlist(prev_senders_list[rows], use.names = FALSE),
                    unlist(prev_senders_strict_list[rows], use.names = FALSE)))
    co <- crossprod(D[in_window, rel, drop = FALSE])
    dimnames(co) <- list(rel, rel)

    for (i in rows) {
      all_prev_senders <- prev_senders_list[[i]]
      if (length(all_prev_senders) == 0) next
      v_4cycle[i] <- mean(co[dt$SOURCE[i], all_prev_senders])

      # STRICT PRECEDENCE: only strictly preceding sponsors (TIME < time);
      # computed only when the row has any previous senders (as before).
      prev_senders_strict <- prev_senders_strict_list[[i]]
      if (length(prev_senders_strict) > 0) {
        v_4cycle_s[i] <- mean(co[dt$SOURCE[i], prev_senders_strict])
      }
    }
  }

  dt[, `:=`(dyad_4cycle_sponsor1year = v_4cycle,
            dyad_4cycle_sponsor1year_strict = v_4cycle_s)]

  data.table::fwrite(dt, file = out_csv_path)
  invisible(out_csv_path)
}
