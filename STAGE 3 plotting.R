# ============================================================
# Stage 3 — Export point-level evaluation table from loaded session
# Hatvani I.G. & Kern Z. — Mind the gap (HESS, 2026)
# https://github.com/istvan60/SPbI
#
# Run inside the RStudio session where all_imputed is loaded.
# Does NOT reload the .RData — works with the object in memory.
#
# Input : all_imputed  (object in the Global Environment)
# Output: interp_point_eval.rds
#         interp_perf_summary.csv
# Next  : STAGE 3.2 ICM extraction.R
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------
# USER SETTINGS
# -----------------------------
out_dir <- "/Users/hatvaniistvan/Downloads/tst/SPbI_threshold_extracted/extracted_for_fig2"

target_object <- "all_imputed"

# Optional filters. Leave NULL to keep all.
keep_X <- NULL
# keep_X <- c(0.02, 0.04, 0.08, 0.16, 0.32)

keep_methods <- NULL
# keep_methods <- c("Kalman", "Linear", "LOCF", "Moving-average", "Spline", "Stine", "Sinusoidal", "SPbI")

# Writing the full point-level CSV can be very slow/large. RDS is enough for the next scripts.
write_point_csv <- FALSE

# Set TRUE only after checking that the export worked.
remove_large_object_after_export <- FALSE

# -----------------------------
# Helpers
# -----------------------------
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

first_existing <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) == 0) NA_character_ else hit[1]
}

as_month_date <- function(x) {
  if (inherits(x, "IDate")) x <- as.Date(x)
  if (!inherits(x, "Date")) x <- as.Date(x)
  as.IDate(sprintf(
    "%04d-%02d-01",
    as.integer(format(x, "%Y")),
    as.integer(format(x, "%m"))
  ))
}

stop_with_names <- function(msg, DT) {
  stop(paste0(
    msg, "\n\nAvailable columns are:\n",
    paste(names(DT), collapse = ", ")
  ), call. = FALSE)
}

# -----------------------------
# Use already-loaded all_imputed
# -----------------------------
if (!exists(target_object, envir = .GlobalEnv, inherits = FALSE)) {
  stop(
    "Object '", target_object, "' is not present in the Global Environment.\n",
    "Open the RStudio session where the huge RData is already loaded, then run this script.",
    call. = FALSE
  )
}

message("Using already-loaded object: ", target_object)

# Avoid loading/copying the .RData again.
# setDT modifies by reference and is usually the most memory-efficient option here.
DT <- get(target_object, envir = .GlobalEnv)
setDT(DT)

message("Rows in all_imputed: ", format(nrow(DT), big.mark = ","))
message("Columns in all_imputed: ", paste(names(DT), collapse = ", "))

# -----------------------------
# Harmonise / create required columns
# -----------------------------
if (!"d_ex_orig" %in% names(DT)) {
  if ("d_ex_obs" %in% names(DT)) {
    setnames(DT, "d_ex_obs", "d_ex_orig")
  } else if (all(c("H2_orig", "O18_orig") %in% names(DT))) {
    DT[, d_ex_orig := H2_orig - 8 * O18_orig]
  } else {
    warning("Could not create d_ex_orig: missing d_ex_obs and/or H2_orig/O18_orig.")
  }
}

if (!"d_ex_imp" %in% names(DT)) {
  if (all(c("H2_imp", "O18_imp") %in% names(DT))) {
    DT[, d_ex_imp := H2_imp - 8 * O18_imp]
  } else {
    warning("Could not create d_ex_imp: missing H2_imp and/or O18_imp.")
  }
}

required_basic <- c("method", "X", "O18_orig", "O18_imp", "H2_orig", "H2_imp")
missing_basic <- setdiff(required_basic, names(DT))
if (length(missing_basic) > 0) {
  stop_with_names(
    paste0("Missing required column(s): ", paste(missing_basic, collapse = ", ")),
    DT
  )
}

site_col <- first_existing(
  names(DT),
  c("Site", "site", "Station", "station", "station_name", "Station_name",
    "station_id", "Station_ID", "site_id", "Site_ID", "Name", "name")
)
if (is.na(site_col)) {
  stop_with_names(
    "Could not identify a station/site column. I need this to join the ICM grid values.",
    DT
  )
}
if (site_col != "Site") setnames(DT, site_col, "Site")

date_col <- first_existing(
  names(DT),
  c("date", "Date", "sample_date", "Sample_date", "month_date",
    "Month_date", "time", "Time", "datetime", "DateTime")
)
year_col  <- first_existing(names(DT), c("year", "Year", "YEAR"))
month_col <- first_existing(names(DT), c("month", "Month", "MONTH"))

if (!is.na(date_col)) {
  if (date_col != "date") setnames(DT, date_col, "date")
  DT[, date := as_month_date(date)]
} else if (!is.na(year_col) && !is.na(month_col)) {
  DT[, date := as.IDate(sprintf(
    "%04d-%02d-01",
    as.integer(get(year_col)),
    as.integer(get(month_col))
  ))]
} else {
  stop_with_names(
    "Could not identify date/month information. I need either a date column or year + month columns.",
    DT
  )
}

boot_col <- first_existing(
  names(DT),
  c("boot", "bootstrap", "Bootstrap", "replicate", "rep", "iter", "iteration")
)
if (is.na(boot_col)) {
  DT[, boot := NA_integer_]
} else if (boot_col != "boot") {
  setnames(DT, boot_col, "boot")
}

if (!is.null(keep_X)) {
  DT <- DT[X %in% keep_X]
}
if (!is.null(keep_methods)) {
  DT <- DT[method %in% keep_methods]
}

x_levels <- sort(unique(as.numeric(DT$X)))
DT[, X_pct := factor(
  paste0(round(as.numeric(X) * 100, 0), "%"),
  levels = paste0(round(x_levels * 100, 0), "%")
)]

# -----------------------------
# Build point-level long table
# -----------------------------
make_eval_long <- function(DT, var_name, orig_col, pred_col) {
  if (!all(c(orig_col, pred_col) %in% names(DT))) return(NULL)
  
  out <- DT[, .(
    source = "masked_interpolation",
    var    = var_name,
    Site   = as.character(Site),
    date   = date,
    X      = as.numeric(X),
    X_pct  = as.character(X_pct),
    boot   = boot,
    method = as.character(method),
    orig   = as.numeric(get(orig_col)),
    pred   = as.numeric(get(pred_col))
  )]
  
  out[, diff := pred - orig]
  out[is.finite(orig) & is.finite(pred) & is.finite(diff)]
}

message("Creating point-level evaluation table...")
interp_eval <- rbindlist(
  list(
    make_eval_long(DT, "d18O",     "O18_orig",  "O18_imp"),
    make_eval_long(DT, "d2H",      "H2_orig",   "H2_imp"),
    make_eval_long(DT, "d-excess", "d_ex_orig", "d_ex_imp")
  ),
  use.names = TRUE,
  fill = TRUE
)

message("Rows in extracted point-level table: ", format(nrow(interp_eval), big.mark = ","))

# -----------------------------
# Summary statistics
# First calculate per bootstrap replicate, then summarize across replicates.
# -----------------------------
perf_by_boot <- interp_eval[
  ,
  .(
    MAD  = mean(abs(diff), na.rm = TRUE),
    RMSE = sqrt(mean(diff^2, na.rm = TRUE)),
    bias = mean(diff, na.rm = TRUE),
    n    = .N
  ),
  by = .(source, var, X, X_pct, boot, method)
]

perf_summary <- perf_by_boot[
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

ba_stats <- interp_eval[
  ,
  .(
    bias = mean(diff, na.rm = TRUE),
    n_points = .N
  ),
  by = .(source, var, X, X_pct, method)
]

# -----------------------------
# Save small files
# -----------------------------
saveRDS(interp_eval,    file.path(out_dir, "interp_point_eval.rds"), compress = FALSE)
saveRDS(perf_by_boot,  file.path(out_dir, "interp_perf_by_boot.rds"), compress = FALSE)
saveRDS(perf_summary,  file.path(out_dir, "interp_perf_summary.rds"), compress = FALSE)
saveRDS(ba_stats,      file.path(out_dir, "interp_ba_stats.rds"), compress = FALSE)

fwrite(perf_summary, file.path(out_dir, "interp_perf_summary.csv"))
if (isTRUE(write_point_csv)) {
  fwrite(interp_eval, file.path(out_dir, "interp_point_eval.csv"))
}

message("Done. Saved extracted files to: ", normalizePath(out_dir))
message("Next: close/restart RStudio, then run 02_extract_ICM_and_join_to_interpolation.R and 03_plot_combined_BA_with_ICM.R")

if (isTRUE(remove_large_object_after_export)) {
  rm(list = target_object, envir = .GlobalEnv)
  gc()
  message("Removed ", target_object, " from the Global Environment and ran gc().")
}


