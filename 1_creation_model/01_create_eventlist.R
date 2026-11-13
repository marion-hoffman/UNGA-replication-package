# Create event list for the INITIATION analysis
# Observed events = draft creation events (cells == 1 in country_draft_abs);
# event, ALL states are at risk of being the creator -> full enumeration of
# the risk set, no non-event sampling and therefore no sampling offset
# (num_sampled = num_possible, weight = 1, offset = 0 for every row).

set.seed(123)
source(file.path(Find(function(d) file.exists(file.path(d, "paths.R")),
                      c(".", "..", "../..")), "paths.R"), chdir = TRUE)

# Rebuild creation events from the join-order matrix in country_draft_abs:
# one row per draft (ldoc_final), cols 30:225 = states, cell = join order
# (1 = creator, 2 = first joiner wave, ...; NA = never sponsored, -9 = not
# a UN member then). Creation events = cells == 1 at time date_or; join
# events = cells >= 2. Restricted to drafts_2017_2020.RDS.
country_draft_abs <- readRDS(COUNTRY_DRAFT_ABS)
drafts_keep <- as.character(readRDS(DRAFTS_KEEP))
state_cols <- names(country_draft_abs)[30:225]
keep <- as.character(country_draft_abs$ldoc_final) %in% drafts_keep
cda <- country_draft_abs[keep, , drop = FALSE]

spons <- do.call(rbind, lapply(seq_len(nrow(cda)), function(r) {
  v <- suppressWarnings(as.numeric(as.matrix(cda[r, state_cols, drop = FALSE])))
  ok <- which(!is.na(v) & v >= 1)
  if (length(ok) == 0) return(NULL)
  data.frame(
    sender   = state_cols[ok],
    receiver = as.character(cda$ldoc_final[r]),
    wave     = v[ok],
    time     = suppressWarnings(as.numeric(as.POSIXct(cda$date_or[r], tz = "UTC"))),
    stringsAsFactors = FALSE
  )
}))
if (is.null(spons) || nrow(spons) == 0) {
  stop("No sponsorship records found in country_draft_abs for the selected drafts")
}
# Earliest wave per (draft, state) in case of duplicate ldoc_final rows
spons <- spons[order(spons$receiver, spons$wave), ]
spons <- spons[!duplicated(paste(spons$receiver, spons$sender)), ]

# Creation events = wave-1 cells; one draft can have several creators ->
# simultaneous events at date_or.
events <- spons[spons$wave == 1, c("sender", "receiver", "time")]
if (any(!is.finite(events$time))) {
  warning(sprintf("%d creation events dropped: missing date_or", sum(!is.finite(events$time))))
  events <- events[is.finite(events$time), ]
}

# State universe = senders of join events (cells >= 2) -- the same ~194
# states used as the sponsoring analysis risk set. States that never
# initiated a draft stay in the risk set.
states <- sort(unique(as.character(spons$sender[spons$wave >= 2])))
# Safety: every observed creator must be in the risk set -- if a creator is
# missing its observed cell would never exist in all_events, so warn loudly.
missing_creators <- setdiff(as.character(events$sender), states)
if (length(missing_creators) > 0) {
  warning(sprintf(
    "%d creation events involve senders outside the %d-state universe (%s) -- these are added to the risk set",
    length(missing_creators), length(states), paste(missing_creators, collapse = ", ")))
  states <- sort(union(states, missing_creators))
}
n_states <- length(states)

# Prepare observed events. ORDER is a plain event id sorted by time -- no
# randomization needed: all covariates are exogenous to the ordering of
# concurrent events, and estimation uses coxph with tie handling instead of
# per-event strata.
n_events <- nrow(events)
events <- events[order(events$time), ]
events$IS_OBSERVED <- 1L
events$SOURCE <- as.character(events$sender)
events$TIME <- as.integer(events$time)
events$ORDER <- seq_len(n_events)
# No TARGET column: the receiver is not modelled directly, but it IS needed
# to look up the created draft's topic (HR condemnation vs other).

# Look up whether each created draft is an HR-condemnation draft. Event
# receivers are ldoc_final; the topic table is keyed by the ORIGINAL
# title (ldoc_or), so its key is translated to ldoc_final via the
# all-drafts table country_draft_abs (which holds both). Drafts whose
# ldoc_or has no ldoc_final match get NA and are dropped below.
df_topic <- readRDS(TOPIC_RDS)
stopifnot("ldoc_or" %in% names(df_topic),
          all(c("ldoc_or", "ldoc_final") %in% names(country_draft_abs)))
draft_map <- setNames(as.character(country_draft_abs$ldoc_final),
                      as.character(country_draft_abs$ldoc_or))
if (anyDuplicated(names(draft_map)))
  warning("country_draft_abs has duplicated ldoc_or values")
topic_key <- draft_map[as.character(df_topic$ldoc_or)]
if (anyNA(topic_key))
  cat(sprintf("%d topic-table drafts have no ldoc_final match and are unusable\n",
              sum(is.na(topic_key))))
hr_lookup <- setNames(df_topic$topic_hr_violation, topic_key)
events$event_topic_hr <- hr_lookup[as.character(events$receiver)]
# Drop events whose draft has no topic entry -- their observed cell is
# undefined
n_dropped <- sum(is.na(events$event_topic_hr))
if (n_dropped > 0) {
  cat(sprintf("Dropping %d creation events with unknown draft topic\n", n_dropped))
  events <- events[!is.na(events$event_topic_hr), ]
  events$ORDER <- seq_len(nrow(events))
  n_events <- nrow(events)
}
events$event_topic_hr <- as.integer(events$event_topic_hr)

# Expand: risk set = all states x topic values (alter_topic_hr in {0,1}),
# so each (state, topic) cell is one candidate. The observed cell is the
# creator x the actual topic of the created draft.
all_events <- data.frame(
  IS_OBSERVED = 0L,
  SOURCE = rep(rep(states, each = 2L), times = n_events),
  alter_topic_hr = rep(rep(c(0L, 1L), times = n_states), times = n_events),
  TIME = rep(events$TIME, each = n_states * 2L),
  ORDER = rep(events$ORDER, each = n_states * 2L),
  stringsAsFactors = FALSE
)
key_all <- paste(all_events$ORDER, all_events$SOURCE, all_events$alter_topic_hr)
key_obs <- paste(events$ORDER, events$SOURCE, events$event_topic_hr)
all_events$IS_OBSERVED[match(key_obs, key_all)] <- 1L

# Full enumeration: no sampling weight needed (kept for schema parity with
# the sponsoring eventlist)
all_events$num_possible <- n_states * 2L
all_events$num_sampled <- n_states * 2L
all_events$weight <- 1
all_events$offset <- 0

# Observed row first within each stratum
all_events <- all_events[order(all_events$ORDER, -all_events$IS_OBSERVED), ]

write.csv(all_events, file = CREATION_EVENTLIST, row.names = FALSE)
cat(sprintf("Created eventlist_for_rem.csv: %d creation events x %d states x 2 topics = %d rows\n",
            n_events, n_states, nrow(all_events)))
