add_interaction_effects <- function(
  all_events_csv,
  out_csv_path = all_events_csv
) {
  all_events <- read.csv(all_events_csv, stringsAsFactors = FALSE)
  if (!all(c("ego_hr", "alter_topic_hr") %in% colnames(all_events))) {
    stop("Missing required columns ego_hr or alter_topic_hr in eventlist.")
  }
  all_events$ego_hrXalter_topic_hr <- all_events$ego_hr * all_events$alter_topic_hr
  dir.create(dirname(out_csv_path), showWarnings = FALSE, recursive = TRUE)
  write.csv(all_events, file = out_csv_path, row.names = FALSE)
  invisible(out_csv_path)
}
