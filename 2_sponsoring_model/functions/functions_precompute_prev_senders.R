precompute_prev_senders <- function(
  all_events_csv,
  country_draft_abs_path,
  out_rds_path
) {
  all_events <- read.csv(all_events_csv, stringsAsFactors = FALSE)
  country_draft_abs <- readRDS(country_draft_abs_path)

  n <- nrow(all_events)

  # Abstract-stage senders (draft creators) per draft, computed once.
  state_cols <- names(country_draft_abs)[30:225]
  draft_mat <- as.matrix(country_draft_abs[, state_cols])
  abs_by_draft <- lapply(seq_len(nrow(draft_mat)),
                         function(r) state_cols[which(draft_mat[r, ] == 1)])
  names(abs_by_draft) <- as.character(country_draft_abs$ldoc_final)

  # Previous senders are defined PER ROW: for row i they are the senders of
  # observed events on the ROW's TARGET draft (the sampled candidate for
  # non-event rows, not the focal draft) with ORDER < the focal event's
  # order -- plus the candidate draft's abstract-stage creators. So the set
  # depends on (TARGET_i, ORDER_i); rows of the same event can differ.
  # Fast version: per TARGET draft, build cumulative sender sets once
  # (indexed by ORDER for non-strict, TIME for strict), then a cheap
  # count-based lookup per row.
  obs <- all_events[all_events$IS_OBSERVED == 1, ]

  prev_map <- new.env(hash = TRUE, parent = emptyenv())   # by ORDER < o
  strict_map <- new.env(hash = TRUE, parent = emptyenv()) # by TIME < t

  for (d in unique(all_events$TARGET)) {
    od <- obs[obs$TARGET == d, ]
    abs_s <- abs_by_draft[[d]]
    if (is.null(abs_s)) abs_s <- character(0)

    # cum[[k+1]] = creators + senders of the first k events (by ORDER)
    oo <- od[order(od$ORDER), ]
    ords <- oo$ORDER
    cum <- vector("list", length(ords) + 1)
    cum[[1]] <- abs_s
    seen <- abs_s
    for (r in seq_len(nrow(oo))) {
      if (!(oo$SOURCE[r] %in% seen)) seen <- c(seen, oo$SOURCE[r])
      cum[[r + 1]] <- seen
    }
    assign(d, list(ords = ords, cum = cum), envir = prev_map)

    # cumt[[j+1]] = creators + senders of events at the first j unique TIMEs
    ut <- sort(unique(od$TIME))
    cumt <- vector("list", length(ut) + 1)
    cumt[[1]] <- abs_s
    seen <- abs_s
    for (j in seq_along(ut)) {
      seen <- union(seen, od$SOURCE[od$TIME == ut[j]])
      cumt[[j + 1]] <- seen
    }
    assign(d, list(times = ut, cum = cumt), envir = strict_map)
  }

  prev_senders_list <- vector("list", n)
  prev_senders_strict_list <- vector("list", n)
  tgt <- as.character(all_events$TARGET)
  for (i in seq_len(n)) {
    e <- prev_map[[tgt[i]]]
    if (is.null(e)) {
      prev_senders_list[[i]] <- character(0)
    } else {
      prev_senders_list[[i]] <- e$cum[[sum(e$ords < all_events$ORDER[i]) + 1]]
    }
    e <- strict_map[[tgt[i]]]
    if (is.null(e)) {
      prev_senders_strict_list[[i]] <- character(0)
    } else {
      prev_senders_strict_list[[i]] <- e$cum[[sum(e$times < all_events$TIME[i]) + 1]]
    }
  }

  dir.create(dirname(out_rds_path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(prev_senders_list, file = out_rds_path)
  saveRDS(prev_senders_strict_list,
          file = file.path(dirname(out_rds_path), "prev_senders_strict_list.RDS"))
  invisible(out_rds_path)
}
