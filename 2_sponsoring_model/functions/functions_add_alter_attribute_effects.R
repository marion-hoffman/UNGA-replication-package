add_alter_attribute_effects <- function(
  all_events_csv,
  country_draft_abs_path,
  recur_title_rds,
  consensus_rds,
  topic_rds,
  nrefs_rds,
  out_csv_path = all_events_csv
) {
  all_events <- read.csv(all_events_csv, stringsAsFactors = FALSE)
  country_draft_abs <- readRDS(country_draft_abs_path)
  df_recur_title <- readRDS(recur_title_rds)
  df_cons <- readRDS(consensus_rds)
  df_topic <- readRDS(topic_rds)
  df_nrefs <- readRDS(nrefs_rds)

  all_events$alter_recur <- 0
  all_events$alter_cons <- 0
  all_events$alter_topic_hr <- 0
  all_events$alter_topic_coop <- 0
  all_events$alter_topic_ispa <- 0
  all_events$alter_nrefs <- 0

  # Topic/recurrence tables are keyed by the ORIGINAL title ldoc_or;
  # drafts are identified everywhere else by their final title
  # ldoc_final. Translate via the all-drafts table (holds both);
  # drafts whose ldoc_or has no ldoc_final match stay NA (a few
  # titles never line up).
  stopifnot(all(c("ldoc_or", "ldoc_final") %in% names(country_draft_abs)))
  draft_map <- setNames(as.character(country_draft_abs$ldoc_final),
                        as.character(country_draft_abs$ldoc_or))
  if (anyDuplicated(names(draft_map)))
    warning("country_draft_abs has duplicated ldoc_or values")
  key_final <- function(d) {
    if ("ldoc_or" %in% names(d)) draft_map[as.character(d$ldoc_or)]
    else as.character(d$ldoc_final)
  }

  recur_ldoc_final <- unique(na.omit(key_final(df_recur_title[df_recur_title$reoccur == 1, , drop = FALSE])))
  cons_lookup <- setNames(df_cons$lag_cons, as.character(df_cons$ldoc_final))
  hr_lookup <- setNames(df_topic$topic_hr_violation, key_final(df_topic))
  coop_lookup <- setNames(df_topic$topic_int_coop, key_final(df_topic))
  ispa_lookup <- setNames(df_topic$topic_israel_palestine, key_final(df_topic))
  nrefs_lookup <- setNames(df_nrefs$nrefs, as.character(df_nrefs$ldoc_final))

  for (i in seq_len(nrow(all_events))) {
    draft <- as.character(all_events$TARGET[i])
    all_events$alter_recur[i] <- as.integer(draft %in% recur_ldoc_final)
    all_events$alter_cons[i] <- cons_lookup[draft]
    all_events$alter_topic_hr[i] <- hr_lookup[draft]
    all_events$alter_topic_coop[i] <- coop_lookup[draft]
    all_events$alter_topic_ispa[i] <- ispa_lookup[draft]
    all_events$alter_nrefs[i] <- log1p(nrefs_lookup[draft])
  }
  all_events$alter_cons[is.na(all_events$alter_cons)] <- 0

  dir.create(dirname(out_csv_path), showWarnings = FALSE, recursive = TRUE)
  write.csv(all_events, file = out_csv_path, row.names = FALSE)
  invisible(out_csv_path)
}
