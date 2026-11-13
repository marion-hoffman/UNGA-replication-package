# Add ego endogeneity effect to eventlist for the initiation REM
# ego_sponsor1year = number of drafts the candidate country sponsored in the
# previous year (identical calculation as in the sponsoring analysis:
# drafts with date_or within the last 365 days and date_final before t).

source(file.path(Find(function(d) file.exists(file.path(d, "paths.R")),
                      c(".", "..", "../..")), "paths.R"), chdir = TRUE)

# Load eventlist
all_events <- read.csv(CREATION_EVENTLIST, stringsAsFactors = FALSE)
n <- nrow(all_events)
# Load country_draft_abs (country x draft sponsorship matrix with draft dates)
country_draft_abs <- readRDS(COUNTRY_DRAFT_ABS)

# Add ego_sponsor1year column
all_events$ego_sponsor1year <- NA_integer_

# For each row, count number of times the candidate country sponsored a
# draft in the previous year. Same definition as sponsoring analysis:
# sponsorship = non-NA entry in the country column, counted when the draft
# was created within the last year (date_or >= t-1y) and already finished
# before t (date_final < t).
for (i in seq_len(n)) {
  if (i %% max(1, floor(n/100)) == 0) cat(sprintf("%d%% done\n", round(100*i/n)))
  sender <- all_events$SOURCE[i]
  event_time <- all_events$TIME[i]
  time_window_start <- event_time - 365*24*60*60
  if (!(sender %in% colnames(country_draft_abs))) {
    # State absent from the sponsorship matrix = it sponsored zero drafts
    # in the window (full state universe is at risk here, unlike sponsoring)
    all_events$ego_sponsor1year[i] <- 0L
    next
  }
  sponsored <- sum(
    as.integer(as.POSIXct(country_draft_abs$date_or, origin="1970-01-01", tz="UTC")) >= time_window_start &
    as.integer(as.POSIXct(country_draft_abs$date_final, origin="1970-01-01", tz="UTC")) < event_time &
    !is.na(country_draft_abs[[sender]]), na.rm = TRUE)
  all_events$ego_sponsor1year[i] <- sponsored
}

# Re-export file
cat("Remaining NAs per column:\n")
print(colSums(is.na(all_events)))

write.csv(all_events, file = CREATION_EVENTLIST, row.names = FALSE)
cat("Added ego_sponsor1year and exported to eventlist_for_rem.csv\n")
