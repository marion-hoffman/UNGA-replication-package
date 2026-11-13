add_ego_attribute_effects <- function(
  all_events_csv,
  gdp_csv,
  conflict_csv,
  hr_csv,
  hr_ciri_csv,
  staff_csv,
  group_csv,
  regime_csv,
  out_csv_path = all_events_csv
) {
  all_events <- read.csv(all_events_csv, stringsAsFactors = FALSE)
  attr_gdp <- read.csv(gdp_csv, stringsAsFactors = FALSE)
  attr_conflict <- read.csv(conflict_csv, stringsAsFactors = FALSE)
  attr_hr <- read.csv(hr_csv, stringsAsFactors = FALSE)
  attr_hr_ciri <- read.csv(hr_ciri_csv, stringsAsFactors = FALSE)
  attr_staff <- read.csv(staff_csv, stringsAsFactors = FALSE)
  attr_group <- read.csv(group_csv, stringsAsFactors = FALSE)
  attr_regime <- read.csv(regime_csv, stringsAsFactors = FALSE)

  all_events$ego_log_gdp <- NA
  all_events$ego_conflict <- NA
  all_events$ego_hr <- NA
  all_events$ego_hr_ciri <- NA
  all_events$ego_staff <- NA
  all_events$ego_group <- NA
  all_events$ego_dem <- NA

  get_year <- function(sec) as.integer(format(as.POSIXct(sec, origin="1970-01-01", tz="UTC"), "%Y"))
  years <- get_year(all_events$TIME)

  med_gdp <- tapply(attr_gdp$gdp_log, attr_gdp$year, median, na.rm=TRUE)
  for (i in seq_len(nrow(all_events))) {
    yr <- years[i]
    iso <- all_events$SOURCE[i]
    val <- attr_gdp$gdp_log[attr_gdp$iso == iso & attr_gdp$year == yr]
    if (length(val) == 0) val <- med_gdp[as.character(yr)]
    all_events$ego_log_gdp[i] <- val[1]
  }

  med_conflict <- tapply(attr_conflict$intensity_level, attr_conflict$year, median, na.rm=TRUE)
  for (i in seq_len(nrow(all_events))) {
    yr <- years[i]
    iso <- all_events$SOURCE[i]
    val <- attr_conflict$intensity_level[attr_conflict$side_a == iso & attr_conflict$year == yr]
    if (length(val) == 0) val <- 0
    all_events$ego_conflict[i] <- val[1]
  }

  attr_hr <- attr_hr[!is.na(attr_hr$iso),]
  med_hr <- tapply(attr_hr$hr_score, attr_hr$year, median, na.rm=TRUE)
  for (i in seq_len(nrow(all_events))) {
    yr <- years[i]
    iso <- all_events$SOURCE[i]
    val <- attr_hr$hr_score[attr_hr$iso == iso & attr_hr$year == yr]
    if (length(val) == 0) val <- med_hr[as.character(yr)]
    all_events$ego_hr[i] <- val[1]
  }
  
  attr_hr_ciri <- attr_hr_ciri[!is.na(attr_hr_ciri$iso),]
  med_hr_ciri <- tapply(attr_hr_ciri$hr_score, attr_hr_ciri$year, median, na.rm=TRUE)
  for (i in seq_len(nrow(all_events))) {
    yr <- years[i]
    iso <- all_events$SOURCE[i]
    val <- attr_hr_ciri$hr_score[attr_hr_ciri$iso == iso & attr_hr_ciri$year == yr]
    if (length(val) == 0) val <- med_hr_ciri[as.character(yr)]
    all_events$ego_hr_ciri[i] <- val[1]
  }

  med_staff <- tapply(attr_staff$mission_size_ordin, attr_staff$year, median, na.rm=TRUE)
  for (i in seq_len(nrow(all_events))) {
    yr <- years[i]
    iso <- all_events$SOURCE[i]
    val <- attr_staff$mission_size_ordin[attr_staff$iso == iso & attr_staff$year == yr]
    if (length(val) == 0) val <- med_staff[as.character(yr)]
    all_events$ego_staff[i] <- val[1]
  }

  med_group <- tapply(attr_group$group_any, attr_group$year, median, na.rm=TRUE)
  for (i in seq_len(nrow(all_events))) {
    yr <- years[i]
    iso <- all_events$SOURCE[i]
    val <- attr_group$group_any[attr_group$iso == iso & attr_group$year == yr]
    if (length(val) == 0) val <- med_group[as.character(yr)]
    all_events$ego_group[i] <- val[1]
  }

  med_regime <- tapply(attr_regime$bmr_dem, attr_regime$year, median, na.rm=TRUE)
  for (i in seq_len(nrow(all_events))) {
    yr <- years[i]
    iso <- all_events$SOURCE[i]
    val <- attr_regime$bmr_dem[attr_regime$ccodealp == iso & attr_regime$year == yr]
    if (length(val) == 0) val <- med_regime[as.character(yr)]
    all_events$ego_dem[i] <- val[1]
  }

  dir.create(dirname(out_csv_path), showWarnings = FALSE, recursive = TRUE)
  write.csv(all_events, file = out_csv_path, row.names = FALSE)
  invisible(out_csv_path)
}
