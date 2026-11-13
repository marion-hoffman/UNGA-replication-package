add_ego_endo_effects <- function(
  all_events_csv,
  country_draft_abs_path,
  out_csv_path = all_events_csv
) {
  all_events <- read.csv(all_events_csv, stringsAsFactors = FALSE)
  country_draft_abs <- readRDS(country_draft_abs_path)

  all_events$ego_sponsor1year <- NA_integer_
  for (i in seq_len(nrow(all_events))) {
    sender <- all_events$SOURCE[i]
    event_time <- all_events$TIME[i]
    time_window_start <- event_time - 365*24*60*60
    if (!(sender %in% colnames(country_draft_abs))) {
      all_events$ego_sponsor1year[i] <- NA
      next
    }
    sponsored <- sum(
      as.integer(as.POSIXct(country_draft_abs$date_or, origin="1970-01-01", tz="UTC")) >= time_window_start &
        as.integer(as.POSIXct(country_draft_abs$date_final, origin="1970-01-01", tz="UTC")) < event_time &
        !is.na(country_draft_abs[[sender]])
    )
    all_events$ego_sponsor1year[i] <- sponsored
  }

  dir.create(dirname(out_csv_path), showWarnings = FALSE, recursive = TRUE)
  write.csv(all_events, file = out_csv_path, row.names = FALSE)
  invisible(out_csv_path)
}
