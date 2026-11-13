create_eventlist <- function(
  N_NON_EVENTS = 10,
  drafts_list_path,
  country_draft_abs_path,
  out_csv_path,
  seed = 123
) {
  set.seed(seed)
  country_draft_abs <- readRDS(country_draft_abs_path)
  drafts_keep <- as.character(readRDS(drafts_list_path))

  # --- Rebuild join events from the join-order matrix ---------------------
  # country_draft_abs: one row per draft (ldoc_final); cols 30:225 are
  # states and each cell holds the join order: 1 = creator, 2 = first
  # joiner wave, 3 = second wave, ...; NA = never sponsored, -9 = not a UN
  # member at the time. Join events = cells >= 2, restricted to drafts in
  # drafts_2017_2020.RDS (the object covers a longer period than studied).
  state_cols <- names(country_draft_abs)[30:225]
  keep <- as.character(country_draft_abs$ldoc_final) %in% drafts_keep
  cda <- country_draft_abs[keep, , drop = FALSE]

  # Draft-open windows: a draft is available in (date_or, date_final].
  # Same source as the imputed wave times, so the focal draft is always
  # open at its own event time; replaces the draft_presence event log.
  draft_ids <- as.character(cda$ldoc_final)
  draft_lo  <- suppressWarnings(as.numeric(as.POSIXct(cda$date_or,    tz = "UTC")))
  draft_hi  <- suppressWarnings(as.numeric(as.POSIXct(cda$date_final, tz = "UTC")))

  ev <- do.call(rbind, lapply(seq_len(nrow(cda)), function(r) {
    v <- suppressWarnings(as.numeric(as.matrix(cda[r, state_cols, drop = FALSE])))
    ok <- which(!is.na(v) & v >= 2)
    if (length(ok) == 0) return(NULL)
    data.frame(
      sender   = state_cols[ok],
      receiver = as.character(cda$ldoc_final[r]),
      wave     = v[ok],
      t_lo     = suppressWarnings(as.numeric(as.POSIXct(cda$date_or[r],    tz = "UTC"))),
      t_hi     = suppressWarnings(as.numeric(as.POSIXct(cda$date_final[r], tz = "UTC"))),
      stringsAsFactors = FALSE
    )
  }))
  if (is.null(ev) || nrow(ev) == 0) stop("No join events found in country_draft_abs for the selected drafts")
  # A state can appear in several rows of the same ldoc_final: keep earliest wave
  ev <- ev[order(ev$receiver, ev$wave), ]
  ev <- ev[!duplicated(paste(ev$receiver, ev$sender)), ]

  # --- Monte-Carlo timing imputation (approach 1: batched waves) ----------
  # For each draft with J join waves, draw J distinct uniform times in
  # (date_or, date_final]; wave j takes the j-th smallest drawn time, so
  # wave order is preserved and wave-mates stay simultaneous. The tie-break
  # below then randomizes within-wave order. Same-day drafts (date_final
  # <= date_or) collapse all joiners to a single time; drafts with no
  # usable window are dropped with a warning.
  ev <- do.call(rbind, lapply(split(ev, ev$receiver), function(d) {
    u <- sort(unique(d$wave))
    t_lo <- suppressWarnings(min(d$t_lo)); t_hi <- suppressWarnings(max(d$t_hi))
    if (is.finite(t_lo) && is.finite(t_hi)) {
      if (t_hi > t_lo) {
        wt <- sort(stats::runif(length(u), min = t_lo, max = t_hi))
        d$time <- wt[match(d$wave, u)]
      } else {
        d$time <- t_hi
      }
    } else {
      d$time <- NA_real_
    }
    d
  }))
  n_dropped <- sum(!is.finite(ev$time))
  if (n_dropped > 0) {
    warning(sprintf("%d join events dropped: draft window (date_or/date_final) missing", n_dropped))
    ev <- ev[is.finite(ev$time), ]
  }
  events <- ev[, c("sender", "receiver", "time")]

  # Randomize order among events that occur at the exact same time (per-thread reproducible)
  if (!is.null(events$time)) {
    events$.__rand <- stats::runif(nrow(events))
    events <- events[order(events$time, events$.__rand), ]
    events$.__rand <- NULL
  }

  n_events <- nrow(events)
  events$IS_OBSERVED <- 1
  events$SOURCE <- as.character(events$sender)
  events$TARGET <- as.character(events$receiver)
  events$TIME <- as.integer(events$time)
  events$ORDER <- seq_len(n_events)
  observed <- events[, c("IS_OBSERVED", "SOURCE", "TARGET", "TIME", "ORDER")]

  all_states <- unique(events$SOURCE)
  non_events_list <- vector("list", n_events)

  for (i in seq_len(n_events)) {
    obs_source <- observed$SOURCE[i]
    obs_target <- observed$TARGET[i]
    obs_time   <- observed$TIME[i]
    obs_order  <- observed$ORDER[i]

    # Drafts open at obs_time: (date_or, date_final]. The focal draft
    # qualifies by construction (its wave time was drawn in the window),
    # but as.integer() truncation of TIME can push a draw in
    # (date_or, date_or + 1s] exactly onto date_or -- keep it anyway.
    available_drafts <- unique(draft_ids[
      !is.na(draft_ids) & !is.na(draft_lo) & !is.na(draft_hi) &
      draft_lo < obs_time & draft_hi >= obs_time])
    available_drafts <- unique(c(available_drafts, obs_target))

    possible <- expand.grid(SOURCE = all_states, TARGET = available_drafts, stringsAsFactors = FALSE)
    possible <- possible[!(possible$SOURCE == obs_source & possible$TARGET == obs_target), ]

    prev_pairs <- observed[observed$ORDER < obs_order, c("SOURCE","TARGET")]
    if (nrow(prev_pairs) > 0) {
      prev_keys <- paste(prev_pairs$SOURCE, prev_pairs$TARGET, sep="|")
      pair_keys <- paste(possible$SOURCE, possible$TARGET, sep="|")
      possible <- possible[!(pair_keys %in% prev_keys), ]
    }

    abs_df <- country_draft_abs[country_draft_abs$ldoc_final %in% available_drafts, , drop=FALSE]
    if (nrow(abs_df) > 0) {
      sender_cols <- names(abs_df)[30:225]
      mat <- as.matrix(abs_df[, sender_cols, drop=FALSE])
      rownames(mat) <- abs_df$ldoc_final
      idx <- which(mat == 1, arr.ind = TRUE)
      if (nrow(idx) > 0) {
        forbidden <- data.frame(
          TARGET = rownames(mat)[idx[, 1]],
          SOURCE = colnames(mat)[idx[, 2]],
          stringsAsFactors = FALSE
        )
        fk <- paste(forbidden$SOURCE, forbidden$TARGET, sep = "|")
        pk <- paste(possible$SOURCE, possible$TARGET, sep = "|")
        possible <- possible[!(pk %in% fk), ]
      }
    }

    num_possible <- nrow(possible)
    num_sampled <- min(N_NON_EVENTS, num_possible)
    if (num_possible == 0) {
      sampled <- possible
    } else {
      sampled <- possible[sample(num_possible, num_sampled), ]
    }

    non_events <- data.frame(
      IS_OBSERVED = 0,
      SOURCE = sampled$SOURCE,
      TARGET = sampled$TARGET,
      TIME = obs_time,
      ORDER = obs_order,
      num_possible = num_possible,
      num_sampled = num_sampled,
      weight = ifelse(num_sampled == 0, NA, num_possible / num_sampled),
      offset = ifelse(num_sampled == 0, NA, log(num_possible / num_sampled)),
      stringsAsFactors = FALSE
    )
    non_events_list[[i]] <- non_events
  }

  # Observed events carry weight 1 / offset 0
  observed$num_possible <- NA_integer_
  observed$num_sampled <- NA_integer_
  observed$weight <- 1
  observed$offset <- 0

  non_events_df <- do.call(rbind, non_events_list)
  all_events <- rbind(observed, non_events_df)
  all_events <- all_events[order(all_events$ORDER, -all_events$IS_OBSERVED), ]
  dir.create(dirname(out_csv_path), showWarnings = FALSE, recursive = TRUE)
  write.csv(all_events, file = out_csv_path, row.names = FALSE)
  invisible(out_csv_path)
}
