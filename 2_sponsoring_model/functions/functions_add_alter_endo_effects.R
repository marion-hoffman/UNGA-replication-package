add_alter_endo_effects <- function(
  all_events_csv,
  country_draft_abs_path,
  prev_senders_rds = file.path(dirname(all_events_csv), "prev_senders_list.RDS"),
  out_csv_path = all_events_csv
) {
  all_events <- read.csv(all_events_csv, stringsAsFactors = FALSE)
  prev_senders_list <- readRDS(prev_senders_rds)
  strict_rds <- file.path(dirname(all_events_csv), "prev_senders_strict_list.RDS")
  if (file.exists(strict_rds)) {
    prev_senders_strict_list <- readRDS(strict_rds)
  } else {
    stop("prev_senders_strict_list.RDS not found next to prev_senders_list.RDS -- ",
         "rerun precompute_prev_senders() (it now emits both files).")
  }
  n <- nrow(all_events)
  if (length(prev_senders_list) != n || length(prev_senders_strict_list) != n)
    stop("prev_senders lists do not match eventlist rows (", n,
         ") -- stale RDS? rerun precompute_prev_senders().")

  # alter_sponsor1year is just the size of each event's previous-sender set;
  # strict version uses strictly preceding sponsors (TIME < time). Both sets
  # are per-event and already computed by precompute_prev_senders().
  all_events$alter_sponsor1year <- lengths(prev_senders_list)
  all_events$alter_sponsor1year_strict <- lengths(prev_senders_strict_list)
  dir.create(dirname(out_csv_path), showWarnings = FALSE, recursive = TRUE)
  write.csv(all_events, file = out_csv_path, row.names = FALSE)
  invisible(out_csv_path)
}
