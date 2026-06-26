# ============================================================
# Stage 3.2 — Extract ICM grid values from NetCDF and join to evaluation table
# Hatvani I.G. & Kern Z. — Mind the gap (HESS, 2026)
# https://github.com/istvan60/SPbI
#
# Input : interp_point_eval.rds          (from Stage 3)
#         minimal_SPbI_input_bundle.rds   (from Stage 0, for station coords)
#         Monthly_194709_202403.nc         (IAEA GNIP ICM NetCDF)
# Output: combined_point_eval_with_ICM.rds
#         ICM_NetCDF_join_match_report.csv
# Next  : STAGE 3.3. final plotting.R
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ncdf4)
})

# ------------------------------------------------------------
# USER SETTINGS
# ------------------------------------------------------------

base_dir <- "~/Downloads/tst/SPbI_threshold_extracted/extracted_for_fig2"

interp_eval_file <- file.path(base_dir, "interp_point_eval.rds")

# NetCDF is currently in the same folder according to your test.
nc_file <- file.path(base_dir, "Monthly_194709_202403.nc")

# Stage 0 bundle. Usually one folder above extracted_for_fig2.
station_bundle_rds <- file.path(dirname(base_dir), "minimal_SPbI_input_bundle.rds")

icm_method_name <- "ICM-grid"

# NetCDF variable names from your ncdf4 inspection.
nc_var_d18O <- "d18Op"
nc_var_d2H  <- "d2Hp"
nc_var_dex  <- "dexc"

nc_lon_var  <- "longitude"
nc_lat_var  <- "latitude"
nc_time_dim <- "time"

# If TRUE, save large CSV copies as well as RDS.
# CSV is useful for inspection but can be slow/large.
write_large_csv <- TRUE

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Cannot find file: ", path, call. = FALSE)
  }
}

normalise_site <- function(x) {
  gsub("\\s+", " ", trimws(as.character(x)))
}

first_existing <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) == 0) NA_character_ else hit[1]
}

as_month_idate <- function(x) {
  x <- as.IDate(x)
  as.IDate(sprintf(
    "%04d-%02d-01",
    as.integer(format(x, "%Y")),
    as.integer(format(x, "%m"))
  ))
}

add_months_to_origin <- function(origin_date, month_offsets) {
  origin_date <- as.IDate(origin_date)
  y0 <- as.integer(format(origin_date, "%Y"))
  m0 <- as.integer(format(origin_date, "%m"))
  total_months <- y0 * 12L + (m0 - 1L) + as.integer(month_offsets)
  yy <- total_months %/% 12L
  mm <- total_months %% 12L + 1L
  as.IDate(sprintf("%04d-%02d-01", yy, mm))
}

parse_months_since_origin <- function(units_string) {
  # Expected: "months since 1947-09-01"
  m <- regexec("months since ([0-9]{4}-[0-9]{2}-[0-9]{2})", units_string)
  hit <- regmatches(units_string, m)[[1]]
  if (length(hit) < 2) {
    stop("Could not parse NetCDF time units: ", units_string, call. = FALSE)
  }
  as.IDate(hit[2])
}

nearest_lon_index <- function(site_lon, grid_lon) {
  site_lon <- as.numeric(site_lon)
  grid_lon <- as.numeric(grid_lon)
  
  # If grid is 0..360, convert negative station longitudes to 0..360.
  if (min(grid_lon, na.rm = TRUE) >= 0 && max(grid_lon, na.rm = TRUE) > 180) {
    site_lon_adj <- site_lon %% 360
    d <- abs(grid_lon - site_lon_adj)
    d <- pmin(d, abs(grid_lon - (site_lon_adj + 360)), abs(grid_lon - (site_lon_adj - 360)))
  } else {
    # If grid is -180..180, convert station longitudes >180 if needed.
    site_lon_adj <- ifelse(site_lon > 180, ((site_lon + 180) %% 360) - 180, site_lon)
    d <- abs(grid_lon - site_lon_adj)
  }
  which.min(d)
}

nearest_lat_index <- function(site_lat, grid_lat) {
  which.min(abs(as.numeric(grid_lat) - as.numeric(site_lat)))
}

replace_fill_with_na <- function(values, missval) {
  values <- as.numeric(values)
  if (!is.null(missval) && length(missval) > 0 && is.finite(missval[1])) {
    values[abs(values - as.numeric(missval[1])) < 1e-6] <- NA_real_
  }
  values[!is.finite(values)] <- NA_real_
  values
}

# ------------------------------------------------------------
# File checks
# ------------------------------------------------------------

stop_if_missing(interp_eval_file)
stop_if_missing(station_bundle_rds)
stop_if_missing(nc_file)

dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)

message("Base directory: ", normalizePath(base_dir))
message("Interpolation evaluation file: ", normalizePath(interp_eval_file))
message("Station bundle file: ", normalizePath(station_bundle_rds))
message("NetCDF file: ", normalizePath(nc_file))

# ------------------------------------------------------------
# Load current-run interpolation point evaluation table
# ------------------------------------------------------------

interp_eval <- readRDS(interp_eval_file)
data.table::setDT(interp_eval)

needed_interp <- c("Site", "date", "X", "X_pct", "boot", "method", "var", "orig", "pred", "diff")
missing_interp <- setdiff(needed_interp, names(interp_eval))
if (length(missing_interp) > 0) {
  stop(
    "interp_point_eval.rds is missing column(s): ",
    paste(missing_interp, collapse = ", "),
    "\nAvailable columns: ", paste(names(interp_eval), collapse = ", "),
    call. = FALSE
  )
}

interp_eval[, Site := normalise_site(Site)]
interp_eval[, date := as_month_idate(date)]
interp_eval[, var := as.character(var)]
interp_eval[, method := as.character(method)]
interp_eval[, X_pct := as.character(X_pct)]

if (!"source" %in% names(interp_eval)) {
  interp_eval[, source := "masked_interpolation"]
}

message("Interpolation rows loaded: ", format(nrow(interp_eval), big.mark = ","))
message("Interpolation variables: ", paste(sort(unique(interp_eval$var)), collapse = ", "))
message("Interpolation methods: ", paste(sort(unique(interp_eval$method)), collapse = ", "))

# Unique withheld observations. This avoids attaching ICM once per interpolation method.
# The X and boot columns are intentionally kept because ICM must be evaluated against
# the exact withheld station-month set for each masking fraction and bootstrap replicate.
obs_keys <- unique(interp_eval[, .(Site, date, X, X_pct, boot, var, orig)])
message("Unique withheld station-month observations including X/boot/var: ", format(nrow(obs_keys), big.mark = ","))
message("Unique Site-date-var combinations to extract from NetCDF: ",
        format(nrow(unique(obs_keys[, .(Site, date, var)])), big.mark = ","))

# ------------------------------------------------------------
# Load station coordinates from Stage 0 bundle
# ------------------------------------------------------------

bundle <- readRDS(station_bundle_rds)
if (!"stations_min" %in% names(bundle)) {
  stop("Bundle does not contain stations_min. Available names: ", paste(names(bundle), collapse = ", "), call. = FALSE)
}

stations <- as.data.table(bundle$stations_min)

site_col <- first_existing(names(stations), c("Site", "site", "Station", "station", "Name", "name"))
lon_col  <- first_existing(names(stations), c("Longitude", "longitude", "lon", "Lon", "LON", "x", "X"))
lat_col  <- first_existing(names(stations), c("Latitude", "latitude", "lat", "Lat", "LAT", "y", "Y"))

if (is.na(site_col) || is.na(lon_col) || is.na(lat_col)) {
  stop(
    "Could not identify Site/Longitude/Latitude columns in stations_min.\n",
    "Available columns: ", paste(names(stations), collapse = ", "),
    call. = FALSE
  )
}

setnames(stations, site_col, "Site")
setnames(stations, lon_col,  "station_lon")
setnames(stations, lat_col,  "station_lat")

stations <- stations[, .(
  Site = normalise_site(Site),
  station_lon = as.numeric(station_lon),
  station_lat = as.numeric(station_lat)
)]
stations <- unique(stations[is.finite(station_lon) & is.finite(station_lat)], by = "Site")

needed_sites <- sort(unique(obs_keys$Site))
missing_sites <- setdiff(needed_sites, stations$Site)
if (length(missing_sites) > 0) {
  fwrite(data.table(Site = missing_sites), file.path(base_dir, "ICM_missing_station_coordinates.csv"))
  stop(
    "Missing coordinates for ", length(missing_sites), " site(s). ",
    "Wrote ICM_missing_station_coordinates.csv. First missing sites: ",
    paste(head(missing_sites, 20), collapse = ", "),
    call. = FALSE
  )
}

stations_use <- stations[Site %in% needed_sites]
message("Station coordinates available for all needed sites: ", nrow(stations_use))

# ------------------------------------------------------------
# Open NetCDF and read coordinate/time axes
# ------------------------------------------------------------

nc <- ncdf4::nc_open(nc_file)
on.exit(ncdf4::nc_close(nc), add = TRUE)

nc_vars <- names(nc$var)
required_nc_vars <- c(nc_var_d18O, nc_var_d2H, nc_lon_var, nc_lat_var)
missing_nc_vars <- setdiff(required_nc_vars, nc_vars)
if (length(missing_nc_vars) > 0) {
  stop(
    "NetCDF is missing required variable(s): ", paste(missing_nc_vars, collapse = ", "),
    "\nAvailable variables: ", paste(nc_vars, collapse = ", "),
    call. = FALSE
  )
}

has_dex_nc <- nc_var_dex %in% nc_vars

lon <- as.numeric(ncdf4::ncvar_get(nc, nc_lon_var))
lat <- as.numeric(ncdf4::ncvar_get(nc, nc_lat_var))
time_vals <- as.integer(nc$dim[[nc_time_dim]]$vals)
time_units <- nc$dim[[nc_time_dim]]$units
origin_date <- parse_months_since_origin(time_units)
time_dates <- add_months_to_origin(origin_date, time_vals)

message("NetCDF lon range: ", min(lon, na.rm = TRUE), " to ", max(lon, na.rm = TRUE), " (n=", length(lon), ")")
message("NetCDF lat range: ", min(lat, na.rm = TRUE), " to ", max(lat, na.rm = TRUE), " (n=", length(lat), ")")
message("NetCDF time range: ", min(time_dates), " to ", max(time_dates), " (n=", length(time_dates), ")")
message("NetCDF d-excess variable present: ", has_dex_nc)

# ------------------------------------------------------------
# Assign each station to nearest ICM grid cell
# ------------------------------------------------------------

station_grid <- copy(stations_use)
station_grid[, grid_x_index := vapply(station_lon, nearest_lon_index, integer(1), grid_lon = lon)]
station_grid[, grid_y_index := vapply(station_lat, nearest_lat_index, integer(1), grid_lat = lat)]
station_grid[, grid_lon := lon[grid_x_index]]
station_grid[, grid_lat := lat[grid_y_index]]
station_grid[, abs_lon_diff := abs(grid_lon - ifelse(min(lon) >= 0 && max(lon) > 180, station_lon %% 360, station_lon))]
station_grid[, abs_lat_diff := abs(grid_lat - station_lat)]

fwrite(station_grid, file.path(base_dir, "ICM_station_to_nearest_grid_cell.csv"))
message("Wrote station-grid lookup: ICM_station_to_nearest_grid_cell.csv")

# ------------------------------------------------------------
# Extract ICM time series for needed stations and variables
# ------------------------------------------------------------

needed_vars <- sort(unique(obs_keys$var))
message("Variables requested by current run: ", paste(needed_vars, collapse = ", "))

var_map <- data.table(
  var = c("d18O", "d2H"),
  nc_var = c(nc_var_d18O, nc_var_d2H)
)

if ("d-excess" %in% needed_vars && has_dex_nc) {
  var_map <- rbind(var_map, data.table(var = "d-excess", nc_var = nc_var_dex))
}

var_map <- var_map[var %in% needed_vars]
if (nrow(var_map) == 0) {
  stop("None of the variables in interp_eval can be extracted from the NetCDF. interp vars: ",
       paste(needed_vars, collapse = ", "), call. = FALSE)
}

extract_one_site_var <- function(site_row, var_label, nc_var_name) {
  ix <- as.integer(site_row$grid_x_index)
  iy <- as.integer(site_row$grid_y_index)
  
  vals <- ncdf4::ncvar_get(
    nc,
    nc_var_name,
    start = c(ix, iy, 1),
    count = c(1, 1, -1)
  )
  
  vals <- replace_fill_with_na(vals, nc$var[[nc_var_name]]$missval)
  
  data.table(
    Site = site_row$Site,
    date = time_dates,
    var = var_label,
    pred = vals,
    grid_x_index = ix,
    grid_y_index = iy,
    grid_lon = as.numeric(site_row$grid_lon),
    grid_lat = as.numeric(site_row$grid_lat),
    station_lon = as.numeric(site_row$station_lon),
    station_lat = as.numeric(site_row$station_lat)
  )
}

message("Extracting ICM time series from NetCDF ...")
extract_list <- vector("list", nrow(station_grid) * nrow(var_map))
k <- 1L
for (i in seq_len(nrow(station_grid))) {
  site_row <- station_grid[i]
  for (j in seq_len(nrow(var_map))) {
    extract_list[[k]] <- extract_one_site_var(site_row, var_map$var[j], var_map$nc_var[j])
    k <- k + 1L
  }
}

icm_pred_long <- rbindlist(extract_list, use.names = TRUE, fill = TRUE)
icm_pred_long <- icm_pred_long[is.finite(pred)]

# If d-excess is requested but no NetCDF dexc variable is available, calculate from d2H - 8*d18O.
if ("d-excess" %in% needed_vars && !("d-excess" %in% unique(icm_pred_long$var))) {
  message("Calculating ICM d-excess as d2H - 8*d18O because dexc was not extracted.")
  d18 <- icm_pred_long[var == "d18O", .(Site, date, pred_d18O = pred, grid_x_index, grid_y_index, grid_lon, grid_lat, station_lon, station_lat)]
  d2h <- icm_pred_long[var == "d2H",  .(Site, date, pred_d2H  = pred)]
  dex <- merge(d18, d2h, by = c("Site", "date"), all = FALSE)
  dex[, `:=`(
    var = "d-excess",
    pred = pred_d2H - 8 * pred_d18O
  )]
  dex <- dex[, .(Site, date, var, pred, grid_x_index, grid_y_index, grid_lon, grid_lat, station_lon, station_lat)]
  icm_pred_long <- rbindlist(list(icm_pred_long, dex), use.names = TRUE, fill = TRUE)
}

setkey(icm_pred_long, Site, date, var)
message("ICM prediction rows extracted: ", format(nrow(icm_pred_long), big.mark = ","))

# Save extracted station-grid monthly values for auditability.
saveRDS(icm_pred_long, file.path(base_dir, "ICM_grid_monthly_values_at_stations.rds"), compress = FALSE)
if (isTRUE(write_large_csv)) {
  fwrite(icm_pred_long, file.path(base_dir, "ICM_grid_monthly_values_at_stations.csv"))
}

# ------------------------------------------------------------
# Join ICM values to exact current withheld observations
# ------------------------------------------------------------

icm_eval <- merge(
  obs_keys,
  icm_pred_long,
  by = c("Site", "date", "var"),
  all = FALSE,
  allow.cartesian = TRUE
)

icm_eval[, `:=`(
  source = "ICM_database_grid",
  method = icm_method_name,
  pred = as.numeric(pred),
  orig = as.numeric(orig)
)]
icm_eval[, diff := pred - orig]
icm_eval <- icm_eval[is.finite(orig) & is.finite(pred) & is.finite(diff)]

message("ICM evaluation rows matched to current withheld observations: ", format(nrow(icm_eval), big.mark = ","))

# ------------------------------------------------------------
# Matching diagnostics
# ------------------------------------------------------------

match_report <- merge(
  obs_keys[, .(n_withheld = .N), by = var],
  icm_eval[, .(n_matched = .N), by = var],
  by = "var",
  all = TRUE
)
match_report[is.na(n_withheld), n_withheld := 0L]
match_report[is.na(n_matched), n_matched := 0L]
match_report[, pct_matched := fifelse(n_withheld > 0, 100 * n_matched / n_withheld, NA_real_)]

unique_match_report <- merge(
  unique(obs_keys[, .(Site, date, var)])[, .(n_unique_withheld = .N), by = var],
  unique(icm_eval[, .(Site, date, var)])[, .(n_unique_matched = .N), by = var],
  by = "var",
  all = TRUE
)
unique_match_report[is.na(n_unique_withheld), n_unique_withheld := 0L]
unique_match_report[is.na(n_unique_matched), n_unique_matched := 0L]
unique_match_report[, pct_unique_matched := fifelse(
  n_unique_withheld > 0,
  100 * n_unique_matched / n_unique_withheld,
  NA_real_
)]

match_report_full <- merge(match_report, unique_match_report, by = "var", all = TRUE)
fwrite(match_report_full, file.path(base_dir, "ICM_NetCDF_join_match_report.csv"))
print(match_report_full)


matched_keys <- unique(icm_eval[, .(Site, date, var, X, boot)])
setkey(obs_keys, Site, date, var, X, boot)
setkey(matched_keys, Site, date, var, X, boot)
unmatched_keys <- obs_keys[!matched_keys]
if (nrow(unmatched_keys) > 0) {
  fwrite(unmatched_keys, file.path(base_dir, "ICM_unmatched_current_withheld_keys.csv"))
  warning("Some withheld keys were not matched to ICM. See ICM_unmatched_current_withheld_keys.csv")
}

# Save ICM-only point-level evaluation for auditability.
saveRDS(icm_eval, file.path(base_dir, "ICM_grid_point_eval_current_run.rds"), compress = FALSE)
if (isTRUE(write_large_csv)) {
  fwrite(icm_eval, file.path(base_dir, "ICM_grid_point_eval_current_run.csv"))
}

# ------------------------------------------------------------
# Combine interpolation + ICM as another method
# ------------------------------------------------------------

common_cols <- union(names(interp_eval), names(icm_eval))
for (cc in setdiff(common_cols, names(interp_eval))) interp_eval[, (cc) := NA]
for (cc in setdiff(common_cols, names(icm_eval)))    icm_eval[,    (cc) := NA]
setcolorder(interp_eval, common_cols)
setcolorder(icm_eval, common_cols)

combined_eval <- rbindlist(list(interp_eval, icm_eval), use.names = TRUE, fill = TRUE)
combined_eval <- combined_eval[is.finite(orig) & is.finite(pred) & is.finite(diff)]

message("Combined evaluation rows: ", format(nrow(combined_eval), big.mark = ","))
message("Combined sources: ", paste(sort(unique(combined_eval$source)), collapse = ", "))
message("Combined methods: ", paste(sort(unique(combined_eval$method)), collapse = ", "))

# ------------------------------------------------------------
# Recalculate summaries
# ------------------------------------------------------------

perf_by_boot_all <- combined_eval[
  ,
  .(
    MAD  = mean(abs(diff), na.rm = TRUE),
    RMSE = sqrt(mean(diff^2, na.rm = TRUE)),
    bias = mean(diff, na.rm = TRUE),
    n    = .N
  ),
  by = .(source, var, X, X_pct, boot, method)
]

perf_summary_all <- perf_by_boot_all[
  ,
  .(
    MAD_median  = median(MAD, na.rm = TRUE),
    RMSE_median = median(RMSE, na.rm = TRUE),
    bias_median = median(bias, na.rm = TRUE),
    MAD_mean    = mean(MAD, na.rm = TRUE),
    RMSE_mean   = mean(RMSE, na.rm = TRUE),
    n_points    = sum(n, na.rm = TRUE),
    n_boot      = .N
  ),
  by = .(source, var, X, X_pct, method)
]

ba_stats_all <- combined_eval[
  ,
  .(
    bias = mean(diff, na.rm = TRUE),
    n_points = .N
  ),
  by = .(source, var, X, X_pct, method)
]

# ------------------------------------------------------------
# Save combined outputs
# ------------------------------------------------------------

saveRDS(combined_eval,      file.path(base_dir, "combined_point_eval_with_ICM.rds"), compress = FALSE)
saveRDS(perf_by_boot_all,   file.path(base_dir, "combined_perf_by_boot_with_ICM.rds"), compress = FALSE)
saveRDS(perf_summary_all,   file.path(base_dir, "combined_perf_summary_with_ICM.rds"), compress = FALSE)
saveRDS(ba_stats_all,       file.path(base_dir, "combined_ba_stats_with_ICM.rds"), compress = FALSE)

if (isTRUE(write_large_csv)) {
  fwrite(combined_eval,      file.path(base_dir, "combined_point_eval_with_ICM.csv"))
  fwrite(perf_by_boot_all,   file.path(base_dir, "combined_perf_by_boot_with_ICM.csv"))
  fwrite(perf_summary_all,   file.path(base_dir, "combined_perf_summary_with_ICM.csv"))
  fwrite(ba_stats_all,       file.path(base_dir, "combined_ba_stats_with_ICM.csv"))
} else {
  fwrite(perf_summary_all,   file.path(base_dir, "combined_perf_summary_with_ICM.csv"))
  fwrite(ba_stats_all,       file.path(base_dir, "combined_ba_stats_with_ICM.csv"))
}

message("Done. Files saved in: ", normalizePath(base_dir))
message("Next run: plot combined BA/performance script using combined_point_eval_with_ICM.rds")

