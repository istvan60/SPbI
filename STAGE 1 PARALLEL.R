# ============================================================
# Stage 1 — Spatial proximity-based imputation (SPbI) error extraction
# Hatvani I.G. & Kern Z. — Mind the gap (HESS, 2026)
# https://github.com/istvan60/SPbI
#
# Input : minimal_SPbI_input_bundle.rds  (from Stage 0)
# Output: spbi_nonoverlapping_band_errors.rds
#         paired_SPbI_vs_baselines_nonoverlapping_bands.rds
#         baseline_errors_site_X_boot_method_isotope.rds
# Next  : STAGE 2 CI detection.R
#
# Run this in a fresh RStudio session.
# ============================================================

library(data.table)
library(future.apply)

# -----------------------------
# 0) Settings
# -----------------------------

out_dir   <- "D:/HIG/teszt, hogy megy-e/SPbI_threshold_extracted"

# Start conservatively on Windows.
# Increase only if RAM usage is acceptable.
N_WORKERS <- 8

# Number of Site-X-boot combinations handled by one worker task.
# Smaller chunks = better recovery and smaller temporary results.
# Larger chunks = less scheduling overhead.
CHUNK_SIZE <- 75

# Allow future to export large objects.
# This is not the RAM limit; it only prevents future from refusing large globals.

# Avoid oversubscription: each worker should use one thread internally.
Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)


# If TRUE, chunk results are saved to disk and then combined.
# This is safer for long runs.
SAVE_CHUNKS <- TRUE

chunk_dir <- "~/Downloads/tst/SPbI_threshold_extracted/stage1_parallel_chunks"

dir.create(chunk_dir, showWarnings = FALSE, recursive = TRUE)

# Increase if future complains about exported object size.
# Adjust downward if your RAM is limited.
options(future.globals.maxSize = 20 * 1024^3)  # 20 GB

# Non-overlapping spatial bands
band_width_km <- 50
max_distance_km <- 400

band_edges <- data.table(
  band_low  = seq(0, max_distance_km - band_width_km, by = band_width_km),
  band_high = seq(band_width_km, max_distance_km, by = band_width_km)
)

# Altitude lapse corrections
lapse_d18O <- 1.2 / 1000
lapse_d2H  <- 7.9 / 1000

# -----------------------------
# 1) Load minimal bundle
# -----------------------------

bundle <- readRDS(
  file.path(out_dir, "minimal_SPbI_input_bundle.rds")
)

removed_dates   <- bundle$removed_dates
core_combos     <- bundle$core_combos
baseline_errors <- bundle$baseline_errors
df_O18_min      <- bundle$df_O18_min
df_H2_min       <- bundle$df_H2_min
stations_min    <- bundle$stations_min
dist_km         <- bundle$dist_km

setDT(removed_dates)
setDT(core_combos)
setDT(baseline_errors)
setDT(df_O18_min)
setDT(df_H2_min)
setDT(stations_min)

removed_dates[, Date := as.Date(Date)]
df_O18_min[, Date := as.Date(Date)]
df_H2_min[, Date := as.Date(Date)]

dist_mat <- as.matrix(dist_km)

if (is.null(rownames(dist_mat)) || is.null(colnames(dist_mat))) {
  stop("dist_km must have station names as rownames and colnames.")
}

setindex(removed_dates, Site, X, boot)
setindex(df_O18_min, Site, Date)
setindex(df_H2_min, Site, Date)
setindex(stations_min, Site)

message("Number of Site-X-boot combinations: ", nrow(core_combos))
message("Number of workers: ", N_WORKERS)
message("Chunk size: ", CHUNK_SIZE)

# -----------------------------
# 2) Split core_combos into chunks
# -----------------------------

chunk_id <- ceiling(seq_len(nrow(core_combos)) / CHUNK_SIZE)

combo_chunks <- split(
  core_combos,
  chunk_id
)

message("Number of chunks: ", length(combo_chunks))

# -----------------------------
# 3) Bandwise SPbI function for one chunk
# -----------------------------

calc_spatial_band_errors_chunk <- function(
    combo_chunk,
    iso_dt,
    value_col,
    isotope_name,
    lapse_per_m,
    removed_dates,
    stations_min,
    dist_mat,
    band_edges
) {
  
  setDT(combo_chunk)
  setDT(iso_dt)
  setDT(removed_dates)
  setDT(stations_min)
  
  out <- vector("list", nrow(combo_chunk) * nrow(band_edges))
  out_i <- 1L
  
  for (k in seq_len(nrow(combo_chunk))) {
    
    this_site <- combo_chunk$Site[k]
    this_X    <- combo_chunk$X[k]
    this_boot <- combo_chunk$boot[k]
    
    if (!this_site %in% rownames(dist_mat)) {
      next
    }
    
    alt0 <- stations_min[Site == this_site, Altitude][1]
    
    if (!is.finite(alt0)) {
      next
    }
    
    rmv_dates <- removed_dates[
      Site == this_site & X == this_X & boot == this_boot,
      Date
    ]
    
    rmv_dates <- unique(as.Date(rmv_dates))
    
    if (length(rmv_dates) == 0) {
      next
    }
    
    # Original values at removed dates
    orig_dt <- iso_dt[
      Site == this_site & Date %in% rmv_dates,
      .(Date, orig = get(value_col))
    ]
    
    if (nrow(orig_dt) > 0) {
      orig_dt <- orig_dt[
        ,
        .(orig = mean(orig, na.rm = TRUE)),
        by = Date
      ]
    }
    
    # Candidate neighbour observations at removed dates
    cand_all <- iso_dt[
      Date %in% rmv_dates & Site != this_site,
      .(
        Date,
        Site,
        value = get(value_col),
        Altitude
      )
    ]
    
    dist_vec <- dist_mat[this_site, ]
    
    for (bb in seq_len(nrow(band_edges))) {
      
      lo <- band_edges$band_low[bb]
      hi <- band_edges$band_high[bb]
      band_label <- paste0(lo, "-", hi, " km")
      
      nbrs <- names(dist_vec)[
        is.finite(dist_vec) &
          dist_vec > lo &
          dist_vec <= hi
      ]
      
      if (length(nbrs) == 0 || nrow(cand_all) == 0) {
        
        out[[out_i]] <- data.table(
          Site = this_site,
          X = this_X,
          boot = this_boot,
          isotope = isotope_name,
          band_low = lo,
          band_high = hi,
          band_label = band_label,
          MAD = NA_real_,
          RMSE = NA_real_,
          n_pairs = 0L,
          n_imputed_dates = 0L,
          mean_n_neighbour_values = NA_real_,
          mean_n_neighbour_sites = NA_real_
        )
        
        out_i <- out_i + 1L
        next
      }
      
      cand <- cand_all[
        Site %chin% nbrs & is.finite(value)
      ]
      
      if (nrow(cand) == 0) {
        
        out[[out_i]] <- data.table(
          Site = this_site,
          X = this_X,
          boot = this_boot,
          isotope = isotope_name,
          band_low = lo,
          band_high = hi,
          band_label = band_label,
          MAD = NA_real_,
          RMSE = NA_real_,
          n_pairs = 0L,
          n_imputed_dates = 0L,
          mean_n_neighbour_values = NA_real_,
          mean_n_neighbour_sites = NA_real_
        )
        
        out_i <- out_i + 1L
        next
      }
      
      cand[, corrected_value := value + (Altitude - alt0) * lapse_per_m]
      
      band_means <- cand[
        is.finite(corrected_value),
        .(
          imp = mean(corrected_value, na.rm = TRUE),
          n_neighbour_values = sum(is.finite(corrected_value)),
          n_neighbour_sites = uniqueN(Site)
        ),
        by = Date
      ]
      
      comp <- merge(
        data.table(Date = rmv_dates),
        orig_dt,
        by = "Date",
        all.x = TRUE
      )
      
      comp <- merge(
        comp,
        band_means,
        by = "Date",
        all.x = TRUE
      )
      
      diff <- comp$orig - comp$imp
      n_ok <- sum(is.finite(diff))
      
      out[[out_i]] <- data.table(
        Site = this_site,
        X = this_X,
        boot = this_boot,
        isotope = isotope_name,
        band_low = lo,
        band_high = hi,
        band_label = band_label,
        MAD = if (n_ok > 0) mean(abs(diff), na.rm = TRUE) else NA_real_,
        RMSE = if (n_ok > 0) sqrt(mean(diff^2, na.rm = TRUE)) else NA_real_,
        n_pairs = n_ok,
        n_imputed_dates = sum(is.finite(comp$imp)),
        mean_n_neighbour_values = if (nrow(band_means) > 0) {
          mean(band_means$n_neighbour_values, na.rm = TRUE)
        } else {
          NA_real_
        },
        mean_n_neighbour_sites = if (nrow(band_means) > 0) {
          mean(band_means$n_neighbour_sites, na.rm = TRUE)
        } else {
          NA_real_
        }
      )
      
      out_i <- out_i + 1L
    }
  }
  
  rbindlist(out, use.names = TRUE, fill = TRUE)
}

# -----------------------------
# 4) Parallel runner for one isotope
# -----------------------------

run_isotope_parallel <- function(
    isotope_name,
    iso_dt,
    value_col,
    lapse_per_m,
    combo_chunks,
    removed_dates,
    stations_min,
    dist_mat,
    band_edges,
    chunk_dir
) {
  
  message("Starting parallel extraction for ", isotope_name)
  
  files <- future_lapply(
    seq_along(combo_chunks),
    function(i) {
      
      message(
        isotope_name,
        " chunk ",
        i,
        " / ",
        length(combo_chunks),
        " started at ",
        format(Sys.time(), "%H:%M:%S")
      )
      
      res <- calc_spatial_band_errors_chunk(
        combo_chunk = combo_chunks[[i]],
        iso_dt = iso_dt,
        value_col = value_col,
        isotope_name = isotope_name,
        lapse_per_m = lapse_per_m,
        removed_dates = removed_dates,
        stations_min = stations_min,
        dist_mat = dist_mat,
        band_edges = band_edges
      )
      
      fn <- file.path(
        chunk_dir,
        paste0("spbi_", isotope_name, "_chunk_", sprintf("%04d", i), ".rds")
      )
      
      saveRDS(res, fn, compress = "gzip")
      
      rm(res)
      gc()
      
      message(
        isotope_name,
        " chunk ",
        i,
        " finished at ",
        format(Sys.time(), "%H:%M:%S")
      )
      
      fn
    },
    future.seed = TRUE
  )
  
  files <- unlist(files)
  
  message("Combining chunks for ", isotope_name)
  
  res_all <- rbindlist(
    lapply(files, readRDS),
    use.names = TRUE,
    fill = TRUE
  )
  
  res_all
}

# -----------------------------
# 5) Start parallel plan
# -----------------------------

# For Windows, use multisession.
# On Linux/macOS, multicore can be faster, but multisession is safer in RStudio.
plan(multisession, workers = N_WORKERS)

# -----------------------------
# 6) Run d18O and d2H
# -----------------------------

spbi_d18O <- run_isotope_parallel(
  isotope_name = "d18O",
  iso_dt = df_O18_min,
  value_col = "O18",
  lapse_per_m = lapse_d18O,
  combo_chunks = combo_chunks,
  removed_dates = removed_dates,
  stations_min = stations_min,
  dist_mat = dist_mat,
  band_edges = band_edges,
  chunk_dir = chunk_dir
)

saveRDS(
  spbi_d18O,
  file.path(out_dir, "spbi_nonoverlapping_band_errors_d18O.rds"),
  compress = "gzip"
)

fwrite(
  spbi_d18O,
  file.path(out_dir, "spbi_nonoverlapping_band_errors_d18O.csv.gz")
)

gc()

spbi_d2H <- run_isotope_parallel(
  isotope_name = "d2H",
  iso_dt = df_H2_min,
  value_col = "H2",
  lapse_per_m = lapse_d2H,
  combo_chunks = combo_chunks,
  removed_dates = removed_dates,
  stations_min = stations_min,
  dist_mat = dist_mat,
  band_edges = band_edges,
  chunk_dir = chunk_dir
)

saveRDS(
  spbi_d2H,
  file.path(out_dir, "spbi_nonoverlapping_band_errors_d2H.rds"),
  compress = "gzip"
)

fwrite(
  spbi_d2H,
  file.path(out_dir, "spbi_nonoverlapping_band_errors_d2H.csv.gz")
)

gc()

# Stop workers after calculation
plan(sequential)

# -----------------------------
# 7) Combine isotope results
# -----------------------------

spbi_band_errors <- rbindlist(
  list(spbi_d18O, spbi_d2H),
  use.names = TRUE,
  fill = TRUE
)

saveRDS(
  spbi_band_errors,
  file.path(out_dir, "spbi_nonoverlapping_band_errors.rds"),
  compress = "gzip"
)

fwrite(
  spbi_band_errors,
  file.path(out_dir, "spbi_nonoverlapping_band_errors.csv.gz")
)

# -----------------------------
# 8) Save baseline RDS if not already saved
# -----------------------------

saveRDS(
  baseline_errors,
  file.path(out_dir, "baseline_errors_site_X_boot_method_isotope.rds"),
  compress = "gzip"
)

# -----------------------------
# 9) Create paired SPbI vs baseline table
# -----------------------------

spbi_long <- melt(
  spbi_band_errors,
  id.vars = c(
    "Site", "X", "boot", "isotope",
    "band_low", "band_high", "band_label",
    "n_pairs", "n_imputed_dates",
    "mean_n_neighbour_values", "mean_n_neighbour_sites"
  ),
  measure.vars = c("MAD", "RMSE"),
  variable.name = "metric",
  value.name = "spbi_error"
)

baseline_long <- melt(
  baseline_errors,
  id.vars = c("Site", "X", "boot", "isotope", "method", "n_baseline"),
  measure.vars = c("MAD", "RMSE"),
  variable.name = "metric",
  value.name = "baseline_error"
)

paired_errors <- merge(
  spbi_long,
  baseline_long,
  by = c("Site", "X", "boot", "isotope", "metric"),
  allow.cartesian = TRUE
)

paired_errors[, error_difference := baseline_error - spbi_error]

paired_errors <- paired_errors[
  is.finite(spbi_error) &
    is.finite(baseline_error) &
    is.finite(error_difference)
]

setcolorder(
  paired_errors,
  c(
    "Site", "X", "boot", "isotope", "metric",
    "band_low", "band_high", "band_label",
    "method",
    "spbi_error", "baseline_error", "error_difference",
    "n_pairs", "n_baseline",
    "n_imputed_dates",
    "mean_n_neighbour_values",
    "mean_n_neighbour_sites"
  )
)

saveRDS(
  paired_errors,
  file.path(out_dir, "paired_SPbI_vs_baselines_nonoverlapping_bands.rds"),
  compress = "gzip"
)

fwrite(
  paired_errors,
  file.path(out_dir, "paired_SPbI_vs_baselines_nonoverlapping_bands.csv.gz")
)

message("Done.")
message("Main output:")
message(normalizePath(file.path(out_dir, "paired_SPbI_vs_baselines_nonoverlapping_bands.rds")))
