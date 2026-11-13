# Add ego attribute effects to eventlist for the initiation REM
# Same country-level covariates as the sponsoring analysis, plus
# ego_condemned (fixed indicator: country was condemned in 2012-2020,
# from data/drafts/country_condemnation_variables_final.csv).

source(file.path(Find(function(d) file.exists(file.path(d, "paths.R")),
                      c(".", "..", "../..")), "paths.R"), chdir = TRUE)

# Load eventlist
all_events <- read.csv(CREATION_EVENTLIST, stringsAsFactors = FALSE)
n <- nrow(all_events)

# Load country-level attribute files (same sources as sponsoring)
attr_gdp <- read.csv(GDP_CSV, stringsAsFactors = FALSE)
attr_conflict <- read.csv(CONFLICT_CSV, stringsAsFactors = FALSE)
attr_hr <- read.csv(HR_FARISS_CSV, stringsAsFactors = FALSE)
attr_hr_ciri <- read.csv(HR_CIRI_CSV, stringsAsFactors = FALSE)
attr_staff <- read.csv(STAFF_CSV, stringsAsFactors = FALSE)
attr_group <- read.csv(GROUP_CSV, stringsAsFactors = FALSE)
attr_regime <- read.csv(REGIME_CSV, stringsAsFactors = FALSE)
attr_condemn <- read.csv(CONDEMNATION_CSV, stringsAsFactors = FALSE)

# Add new attribute columns
all_events$ego_log_gdp <- NA
all_events$ego_conflict <- NA
all_events$ego_hr <- NA
all_events$ego_hr_ciri <- NA
all_events$ego_staff <- NA
all_events$ego_group <- NA
all_events$ego_dem <- NA
all_events$ego_condemned <- NA

get_year <- function(sec) as.integer(format(as.POSIXct(sec, origin="1970-01-01", tz="UTC"), "%Y"))
years <- get_year(all_events$TIME)

# GDP: use gdp_log
med_gdp <- tapply(attr_gdp$gdp_log, attr_gdp$year, median, na.rm=TRUE)
for (i in seq_len(n)) {
  if (i %% max(1, floor(n/100)) == 0) cat(sprintf("%d%% done (gdp)\n", round(100*i/n)))
  yr <- years[i]
  iso <- all_events$SOURCE[i]
  val <- attr_gdp$gdp_log[attr_gdp$iso == iso & attr_gdp$year == yr]
  if (length(val) == 0 || is.na(val[1])) val <- med_gdp[as.character(yr)]
  all_events$ego_log_gdp[i] <- val[1]
}

# Conflict: use intensity_level (0 if not present for country-year)
for (i in seq_len(n)) {
  if (i %% max(1, floor(n/100)) == 0) cat(sprintf("%d%% done (conflict)\n", round(100*i/n)))
  yr <- years[i]
  iso <- all_events$SOURCE[i]
  val <- attr_conflict$intensity_level[attr_conflict$side_a == iso & attr_conflict$year == yr]
  if (length(val) == 0 || is.na(val[1])) val <- 0
  all_events$ego_conflict[i] <- val[1]
}

# HR (Fariss): use hr_score
attr_hr <- attr_hr[!is.na(attr_hr$iso),]
med_hr <- tapply(attr_hr$hr_score, attr_hr$year, median, na.rm=TRUE)
for (i in seq_len(n)) {
  if (i %% max(1, floor(n/100)) == 0) cat(sprintf("%d%% done (hr)\n", round(100*i/n)))
  yr <- years[i]
  iso <- all_events$SOURCE[i]
  val <- attr_hr$hr_score[attr_hr$iso == iso & attr_hr$year == yr]
  if (length(val) == 0 || is.na(val[1])) val <- med_hr[as.character(yr)]
  all_events$ego_hr[i] <- val[1]
}

# HR alternative (CIRI): use hr_score
attr_hr_ciri <- attr_hr_ciri[!is.na(attr_hr_ciri$iso),]
med_hr_ciri <- tapply(attr_hr_ciri$hr_score, attr_hr_ciri$year, median, na.rm=TRUE)
for (i in seq_len(n)) {
  if (i %% max(1, floor(n/100)) == 0) cat(sprintf("%d%% done (hr ciri)\n", round(100*i/n)))
  yr <- years[i]
  iso <- all_events$SOURCE[i]
  val <- attr_hr_ciri$hr_score[attr_hr_ciri$iso == iso & attr_hr_ciri$year == yr]
  if (length(val) == 0 || is.na(val[1])) val <- med_hr_ciri[as.character(yr)]
  all_events$ego_hr_ciri[i] <- val[1]
}

# Staff: use mission_size_ordin
med_staff <- tapply(attr_staff$mission_size_ordin, attr_staff$year, median, na.rm=TRUE)
for (i in seq_len(n)) {
  if (i %% max(1, floor(n/100)) == 0) cat(sprintf("%d%% done (staff)\n", round(100*i/n)))
  yr <- years[i]
  iso <- all_events$SOURCE[i]
  val <- attr_staff$mission_size_ordin[attr_staff$iso == iso & attr_staff$year == yr]
  if (length(val) == 0 || is.na(val[1])) val <- med_staff[as.character(yr)]
  all_events$ego_staff[i] <- val[1]
}

# Group: use group_any
med_group <- tapply(attr_group$group_any, attr_group$year, median, na.rm=TRUE)
for (i in seq_len(n)) {
  if (i %% max(1, floor(n/100)) == 0) cat(sprintf("%d%% done (group)\n", round(100*i/n)))
  yr <- years[i]
  iso <- all_events$SOURCE[i]
  val <- attr_group$group_any[attr_group$iso == iso & attr_group$year == yr]
  if (length(val) == 0 || is.na(val[1])) val <- med_group[as.character(yr)]
  all_events$ego_group[i] <- val[1]
}

# Regime: use bmr_dem
med_regime <- tapply(attr_regime$bmr_dem, attr_regime$year, median, na.rm=TRUE)
for (i in seq_len(n)) {
  if (i %% max(1, floor(n/100)) == 0) cat(sprintf("%d%% done (dem)\n", round(100*i/n)))
  yr <- years[i]
  iso <- all_events$SOURCE[i]
  val <- attr_regime$bmr_dem[attr_regime$ccodealp == iso & attr_regime$year == yr]
  if (length(val) == 0 || is.na(val[1])) val <- med_regime[as.character(yr)]
  all_events$ego_dem[i] <- val[1]
}

# Condemned: fixed indicator condemned_2012_2020, keyed by acronym (iso3).
# Countries absent from the file were never condemned -> 0.
for (i in seq_len(n)) {
  if (i %% max(1, floor(n/100)) == 0) cat(sprintf("%d%% done (condemned)\n", round(100*i/n)))
  iso <- all_events$SOURCE[i]
  val <- attr_condemn$condemned_2012_2020[attr_condemn$acronym == iso]
  if (length(val) == 0 || is.na(val[1])) val <- 0
  all_events$ego_condemned[i] <- val[1]
}

# NA sanity check before export
cat("Remaining NAs per column:\n")
print(colSums(is.na(all_events)))

write.csv(all_events, file = CREATION_EVENTLIST, row.names = FALSE)
cat("Added ego attribute effects (incl. ego_condemned) and exported to eventlist_for_rem.csv\n")
