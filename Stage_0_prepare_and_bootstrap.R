# ============================================================
# Stage 0 — Data preparation, imputation bootstrap, bundle export
# Hatvani I.G. & Kern Z. — Mind the gap (HESS, 2026)
# https://github.com/istvan60/SPbI
#
# Input : inputation_test6.xlsx
# Output: minimal_SPbI_input_bundle.rds
# Next  : STAGE 1 PARALLEL.R
# ============================================================
# STAGE 0  —  Data preparation, imputation bootstrap, bundle export
#
# Purpose:
#   Reproduce the full analysis from the raw input file.
#   Run this script once from a fresh R session.
#   All subsequent stages (Phase I, Phase II, 03c) require only
#   the bundle saved at the end of this script.
#
# Input : inputation_test6.xlsx
# Output: minimal_SPbI_input_bundle.rds  (required by Phase I)
#         [optional] diagnostic plots (Bland–Altman, bar charts)
#
# Estimated runtime: ~30–120 min depending on CPU cores and
#   the dynamic n_boot values chosen automatically below.
# ============================================================

# NOTE: kgc depends on plyr, which masks dplyr verbs.
# All dplyr calls in this script use explicit dplyr:: prefixes to avoid conflicts.

# install.packages(c(
#   "readxl", "dplyr", "lubridate", "tidyr", "purrr",
#   "geosphere", "imputeTS", "minpack.lm",
#   "ggplot2", "patchwork", "scales",
#   "furrr", "parallelly", "data.table"
# ))

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(lubridate)
  library(tidyr)
  library(purrr)
  library(geosphere)
  library(imputeTS)
  library(minpack.lm)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(furrr)
  library(parallelly)
  library(data.table)
})

# ============================================================
# USER SETTINGS  —  edit these, leave everything else alone
# ============================================================

input_xlsx <- "/Users/hatvaniistvan/Library/CloudStorage/GoogleDrive-hatvaniig@gmail.com/My Drive/hatvani istvan gabor/cikkek/interpolation paper HESS/rev/inputation_test6.xlsx"

# Folder where all output files will be written.
out_dir <- "/Users/hatvaniistvan/Downloads/tst/SPbI_threshold_extracted"

# Masking fractions to test.
X_vals <- c(0.01, 0.02, 0.04, 0.08, 0.16, 0.32)

# Target total imputed rows per method per X value.
# Higher = more stable estimates, longer runtime.
N_target <- 200000L

# Minimum continuous months required to include a site.
min_continuous_months <- 84L

# First date to retain (earlier records discarded).
data_start <- as.Date("1973-01-01")

# Lapse-rate corrections (‰ per metre).
lapse_d18O <- 1.2 / 1000
lapse_d2H  <- 7.9 / 1000

# Spatial radius used inside the bootstrap for the "spatial" method (km).
spatial_radius_bootstrap_km <- 100

set.seed(123)

# ============================================================

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ────────────────────────────────────────────────────────────
# §1  Read and prepare data
# ────────────────────────────────────────────────────────────

message("§1  Reading input data ...")

remove_sites <- function(df, sites) df %>% filter(!Site %in% sites)

df_O18 <- read_excel(input_xlsx, sheet = "O18") %>%
  mutate(Date = ymd(Date)) %>%
  remove_sites("VIENNA") %>%
  filter(Date >= data_start)

df_H2 <- read_excel(input_xlsx, sheet = "H2") %>%
  mutate(Date = ymd(Date)) %>%
  remove_sites("VIENNA") %>%
  filter(Date >= data_start)

stations <- read_excel(input_xlsx, sheet = "coords") %>%
  remove_sites("VIENNA")

df_O18 <- df_O18 %>% left_join(stations %>% select(Site, Altitude), by = "Site")
df_H2  <- df_H2  %>% left_join(stations %>% select(Site, Altitude), by = "Site")

# Observed d-excess (needed only for site selection below)
df_de <- inner_join(
  df_O18 %>% select(Site, Date, O18),
  df_H2  %>% select(Site, Date, H2),
  by = c("Site", "Date")
) %>%
  mutate(d_ex_orig = H2 - 8 * O18)

# ────────────────────────────────────────────────────────────
# §2  Select focus sites (≥ min_continuous_months uninterrupted)
# ────────────────────────────────────────────────────────────

message("§2  Selecting long-run sites ...")

compute_longest_run <- function(dates) {
  all_months <- seq(
    floor_date(min(dates), "month"),
    floor_date(max(dates), "month"),
    by = "1 month"
  )
  present <- all_months %in% floor_date(dates, "month")
  r <- rle(present)
  if (any(r$values)) max(r$lengths[r$values]) else 0L
}

# Returns a Date vector of all months that fall inside ANY continuous run
# of >= min_months paired observations for a given site.
qualifying_dates <- function(dates, min_months) {
  all_months <- seq(
    floor_date(min(dates), "month"),
    floor_date(max(dates), "month"),
    by = "1 month"
  )
  present <- all_months %in% floor_date(dates, "month")
  r       <- rle(present)
  ends    <- cumsum(r$lengths)
  starts  <- c(1L, ends[-length(ends)] + 1L)
  good    <- which(r$values & r$lengths >= min_months)
  if (length(good) == 0L) return(as.Date(character(0)))
  do.call(c, lapply(good, function(i) all_months[starts[i]:ends[i]]))
}

long_sites <- df_O18 %>%
  dplyr::group_by(Site) %>%
  dplyr::summarise(longest = compute_longest_run(Date), .groups = "drop") %>%
  dplyr::filter(longest >= min_continuous_months) %>%
  dplyr::pull("Site")

message("  Sites retained: ", length(long_sites))

# For each focus site keep only dates inside qualifying continuous windows.
# If a site has multiple runs >= min_continuous_months, all are retained.
qualifying_dates_per_site <- df_O18 %>%
  semi_join(df_de, by = c("Site", "Date")) %>%
  dplyr::filter(Site %in% long_sites) %>%
  dplyr::group_by(Site) %>%
  dplyr::group_modify(~{
    keep <- qualifying_dates(.x$Date, min_continuous_months)
    dplyr::filter(.x, floor_date(Date, "month") %in% keep)
  }) %>%
  dplyr::ungroup() %>%
  dplyr::select(Site, Date)

message("  Qualifying site-months retained: ", nrow(qualifying_dates_per_site))

df_O18_trim <- df_O18 %>%
  semi_join(qualifying_dates_per_site, by = c("Site", "Date"))

df_H2_trim <- df_H2 %>%
  semi_join(qualifying_dates_per_site, by = c("Site", "Date"))

# ────────────────────────────────────────────────────────────
# §3  Distance matrix
# ────────────────────────────────────────────────────────────

message("§3  Building distance matrix ...")

coords_mat <- as.matrix(stations[, c("Longitude", "Latitude")])
dist_km    <- distm(coords_mat, fun = distHaversine) / 1000
rownames(dist_km) <- stations$Site
colnames(dist_km) <- stations$Site

# ────────────────────────────────────────────────────────────
# §4  [Optional] Station map
#     Requires: maps, metR, sf, rnaturalearth
#     Safe to skip — comment out if packages unavailable.
# ────────────────────────────────────────────────────────────

if (all(c("maps", "metR", "sf", "rnaturalearth") %in%
        rownames(installed.packages()))) {

  suppressPackageStartupMessages({
    library(maps); library(metR); library(sf); library(rnaturalearth)
  })

  message("§4  Drawing station map ...")

  sites_map <- stations %>%
    dplyr::mutate(InTrim = Site %in% long_sites) %>%
    tidyr::drop_na(Longitude, Latitude)

  xlim <- range(sites_map$Longitude, na.rm = TRUE) + c(-0.5, 0.5)
  ylim <- range(sites_map$Latitude,  na.rm = TRUE) + c(-0.5, 0.5)

  lon_w <- floor(xlim[1]);   lon_e <- ceiling(xlim[2])
  lat_s <- floor(ylim[1]);   lat_n <- ceiling(ylim[2])

  bg <- map_data("world", region = c(
    "Austria", "Hungary", "Slovenia", "Slovakia",
    "Czech Republic", "Germany", "Italy", "Croatia",
    "Switzerland", "Poland"
  ))

  topo <- tryCatch(
    metR::GetTopography(
      lon.west  = lon_w, lon.east  = lon_e,
      lat.south = lat_s, lat.north = lat_n,
      resolution = 1 / 10
    ),
    error = function(e) NULL
  )

  if (!is.null(topo)) {
    if ("h" %in% names(topo) && !"height" %in% names(topo))
      topo <- dplyr::rename(topo, height = h)
    if ("x" %in% names(topo) && !"lon" %in% names(topo))
      topo <- dplyr::rename(topo, lon = x)
    if ("y" %in% names(topo) && !"lat" %in% names(topo))
      topo <- dplyr::rename(topo, lat = y)

    sf::sf_use_s2(FALSE)
    bb       <- st_as_sfc(st_bbox(
      c(xmin = xlim[1], xmax = xlim[2], ymin = ylim[1], ymax = ylim[2]),
      crs = 4326
    ))
    land     <- ne_countries(scale = 50, returnclass = "sf")
    land_bb  <- suppressWarnings(st_crop(land, bb))
    land_u   <- suppressWarnings(st_union(land_bb))
    water_sf <- suppressWarnings(st_difference(st_as_sf(bb), land_u))
    coords_w <- st_coordinates(st_cast(water_sf, "MULTIPOLYGON"))
    water_df <- as.data.frame(coords_w)[, c("X", "Y", "L1", "L2")]
    names(water_df) <- c("long", "lat", "g1", "g2")
    water_df$group  <- interaction(water_df$g1, water_df$g2, drop = TRUE)
    water_df        <- water_df[, c("long", "lat", "group")]

    p_map <- ggplot() +
      geom_raster(data = topo, aes(x = lon, y = lat, fill = height),
                  interpolate = TRUE) +
      metR::geom_relief(data = topo, aes(x = lon, y = lat, z = height),
                        alpha = 0.25) +
      geom_polygon(data = water_df, aes(long, lat, group = group),
                   fill = "#9ecae1", color = NA, inherit.aes = FALSE) +
      geom_polygon(data = bg, aes(long, lat, group = group),
                   fill = NA, color = "white", linewidth = 0.5) +
      geom_point(
        data = dplyr::filter(sites_map, !InTrim),
        aes(x = Longitude, y = Latitude, shape = "Stations"),
        size = 2.8, color = "black"
      ) +
      geom_point(
        data = dplyr::filter(sites_map, InTrim),
        aes(x = Longitude, y = Latitude, shape = "Focus sites"),
        size = 2.8, color = "black"
      ) +
      scale_shape_manual(name = NULL, values = c("Stations" = 1, "Focus sites" = 17)) +
      scale_fill_gradientn(
        colours = c("#8c2d04", "#31a354", "#a1d99b", "#e5f5e0"),
        name = "Elevation (m)"
      ) +
      coord_quickmap(xlim = xlim, ylim = ylim) +
      labs(x = "Longitude", y = "Latitude") +
      theme_minimal(base_size = 11) +
      theme(panel.grid = element_blank())

    print(p_map)
    ggsave(file.path(out_dir, "map_focus_sites.png"),
           p_map, width = 8, height = 6, dpi = 300)

  } else {
    message("  Topography download failed — map skipped.")
  }

} else {
  message("§4  Map skipped (maps / metR / sf / rnaturalearth not installed).")
}

# ────────────────────────────────────────────────────────────
# §5  Imputation helper functions
# ────────────────────────────────────────────────────────────

message("§5  Defining imputation functions ...")

safe_fun <- function(fun) {
  function(x) {
    out <- tryCatch(fun(x), error = function(e) rep(NA_real_, length(x)))
    out[!is.finite(out)] <- NA_real_
    out
  }
}

impute_funs <- list(
  # LOCF: leading NAs have no predecessor → stay NA (na_remaining = "keep")
  LOCF = function(x)
    imputeTS::na_locf(x, option = "locf", na_remaining = "keep"),

  # Linear/Spline/Stine: no extrapolation beyond first/last observation → NA.
  # Logic is inlined (not a helper) so furrr workers can find it.
  Linear = function(x) {
    out <- imputeTS::na_interpolation(x, option = "linear")
    obs <- which(!is.na(x))
    if (length(obs) > 0L) {
      if (min(obs) > 1L)           out[seq_len(min(obs) - 1L)]          <- NA_real_
      if (max(obs) < length(out))  out[seq(max(obs) + 1L, length(out))] <- NA_real_
    }
    out
  },
  Spline = function(x) {
    out <- imputeTS::na_interpolation(x, option = "spline")
    obs <- which(!is.na(x))
    if (length(obs) > 0L) {
      if (min(obs) > 1L)           out[seq_len(min(obs) - 1L)]          <- NA_real_
      if (max(obs) < length(out))  out[seq(max(obs) + 1L, length(out))] <- NA_real_
    }
    out
  },
  Stine = function(x) {
    out <- imputeTS::na_interpolation(x, option = "stine")
    obs <- which(!is.na(x))
    if (length(obs) > 0L) {
      if (min(obs) > 1L)           out[seq_len(min(obs) - 1L)]          <- NA_real_
      if (max(obs) < length(out))  out[seq(max(obs) + 1L, length(out))] <- NA_real_
    }
    out
  },

  # Kalman: numerical failures already caught by safe_fun wrapper
  Kalman = imputeTS::na_kalman,

  # Moving average k=5: requires exactly 5 non-NA neighbours on each side;
  # boundary positions and positions adjacent to other masked months → NA
  `Moving-average` = function(x) {
    k   <- 5L
    n   <- length(x)
    out <- x
    for (i in which(is.na(x))) {
      li <- i - k; ri <- i + k
      if (li >= 1L && ri <= n) {
        nbrs <- x[c(li:(i - 1L), (i + 1L):ri)]
        if (!any(is.na(nbrs))) out[i] <- mean(nbrs)
      }
    }
    out
  }
)
impute_funs_safe <- purrr::map(impute_funs, safe_fun)

impute_sin_nls <- function(x, dates) {
  start_m <- floor_date(min(dates, na.rm = TRUE), "month")
  t <- as.numeric(
    interval(start_m, floor_date(dates, "month")) / months(1)
  ) + 1
  obs <- !is.na(x)
  df  <- data.frame(t = t[obs], y = x[obs])
  if (nrow(df) < 3) return(rep(NA_real_, length(x)))
  A0 <- mean(df$y); B0 <- diff(range(df$y)) / 2; phi0 <- 0
  fit <- try(
    nlsLM(
      y ~ A + B * sin(2 * pi * t / 12 + phi),
      data  = df,
      start = list(A = A0, B = B0, phi = phi0),
      lower = c(-Inf, -Inf, -2 * pi),
      upper = c(Inf, Inf, 2 * pi),
      control = nls.lm.control(maxiter = 200)
    ),
    silent = TRUE
  )
  if (!inherits(fit, "nls")) return(rep(NA_real_, length(x)))
  co   <- coef(fit)
  yhat <- co["A"] + co["B"] * sin(2 * pi * t / 12 + co["phi"])
  yhat[!is.finite(yhat)] <- NA_real_
  yhat
}

safe_sin_nls <- function(x, dates) {
  out <- tryCatch(
    impute_sin_nls(x, dates),
    error = function(e) rep(NA_real_, length(x))
  )
  out[!is.finite(out)] <- NA_real_
  out
}

# Spatial imputation: mean of altitude-corrected neighbours within radius.
# Closes over df_O18 and df_H2 (full records) so neighbours contribute any
# available month, not just their own qualifying-window months.
spatial_impute_O18 <- function(x, dates, site) {
  alt0 <- stations$Altitude[stations$Site == site]
  out  <- x
  for (j in which(is.na(x))) {
    nbrs <- stations$Site[
      which(dist_km[site, ] > 0 & dist_km[site, ] <= spatial_radius_bootstrap_km)
    ]
    tmp <- df_O18 %>%
      filter(Site %in% nbrs, Date == dates[j]) %>%
      mutate(O18c = O18 + (Altitude - alt0) * lapse_d18O)
    if (any(!is.na(tmp$O18c))) out[j] <- mean(tmp$O18c, na.rm = TRUE)
  }
  out
}

spatial_impute_H2 <- function(x, dates, site) {
  alt0 <- stations$Altitude[stations$Site == site]
  out  <- x
  for (j in which(is.na(x))) {
    nbrs <- stations$Site[
      which(dist_km[site, ] > 0 & dist_km[site, ] <= spatial_radius_bootstrap_km)
    ]
    tmp <- df_H2 %>%
      filter(Site %in% nbrs, Date == dates[j]) %>%
      mutate(H2c = H2 + (Altitude - alt0) * lapse_d2H)
    if (any(!is.na(tmp$H2c))) out[j] <- mean(tmp$H2c, na.rm = TRUE)
  }
  out
}

# ────────────────────────────────────────────────────────────
# §6  Bootstrap imputation loop (parallel)
# ────────────────────────────────────────────────────────────

message("§6  Setting up bootstrap loop ...")

methods <- c(names(impute_funs), "Sinusoidal", "SPbI")
sites   <- unique(df_O18_trim$Site)

# Per-site cache: qualifying-window data only for the target site.
# SPbI neighbour lookup uses full records (df_O18 / df_H2) separately.
site_cache <- setNames(vector("list", length(sites)), sites)
for (s in sites) {
  tmpO <- df_O18_trim %>% filter(Site == s) %>% arrange(Date)
  tmpH <- df_H2_trim  %>% filter(Site == s) %>% arrange(Date)
  site_cache[[s]] <- list(dates = tmpO$Date, yO = tmpO$O18, yH = tmpH$H2)
}

# Dynamic n_boot: scale iterations so total imputed rows ≈ N_target per method/X
site_lengths     <- vapply(site_cache, function(sc) length(sc$yO), integer(1))
masked_per_boot  <- sapply(X_vals, function(x) sum(ceiling(x * site_lengths)))
n_boot_by_X      <- pmax(1L, ceiling(N_target / masked_per_boot))

boot_plan <- tibble(
  X                       = X_vals,
  n_boot                  = n_boot_by_X,
  approx_rows_per_method  = masked_per_boot * n_boot_by_X
)
message("  Boot plan:")
print(boot_plan)

# Full job grid
job <- map2_dfr(
  X_vals, n_boot_by_X,
  ~ expand_grid(Site = sites, X = .x, boot = seq_len(.y), method = methods)
)

message("  Total jobs: ", nrow(job))

# Single-job runner
do_one <- function(Site, X, boot, method) {
  sc     <- site_cache[[Site]]
  yO     <- sc$yO; yH <- sc$yH; dates <- sc$dates

  seed_val <- as.integer(
    (sum(utf8ToInt(Site)) + round(1e6 * X) * 10007 + boot) %% .Machine$integer.max
  )
  set.seed(seed_val)

  n_rm <- ceiling(X * length(yO))
  idx  <- sample.int(length(yO), n_rm)
  
  yO_na <- yO; yO_na[idx] <- NA_real_
  yH_na <- yH; yH_na[idx] <- NA_real_
  
  n_neighbors <- if (method == "SPbI") {
    length(stations$Site[
      which(dist_km[Site, ] > 0 & dist_km[Site, ] <= spatial_radius_bootstrap_km)
    ])
  } else {
    NA_integer_
  }
  
  yO_imp <- switch(
    method,
    Sinusoidal = safe_sin_nls(yO_na, dates),
    SPbI       = spatial_impute_O18(yO_na, dates, Site),
    impute_funs_safe[[method]](yO_na)
  )
  yH_imp <- switch(
    method,
    Sinusoidal = safe_sin_nls(yH_na, dates),
    SPbI       = spatial_impute_H2(yH_na, dates, Site),
    impute_funs_safe[[method]](yH_na)
  )
  
  O18_imp_vals <- yO_imp[idx]
  H2_imp_vals  <- yH_imp[idx]

  tibble(
    Site        = Site,
    X           = X,
    boot        = boot,
    method      = method,
    n_neighbors = n_neighbors,
    Date        = dates[idx],
    O18_orig    = yO[idx],
    O18_imp     = O18_imp_vals,
    H2_orig     = yH[idx],
    H2_imp      = H2_imp_vals,
    fallback    = is.na(O18_imp_vals) | is.na(H2_imp_vals)
  ) %>%
    mutate(
      d_ex_orig = H2_orig - 8 * O18_orig,
      d_ex_imp  = H2_imp  - 8 * O18_imp
    )
}


workers <- max(1L, availableCores() - 3L)
message("  Using ", workers, " parallel workers ...")
plan(multisession, workers = workers)

all_imputed <- future_pmap_dfr(
  job,
  do_one,
  .options = furrr_options(seed = TRUE)
)

plan(sequential)
message("  all_imputed rows: ", nrow(all_imputed))

# ── Fallback summary ─────────────────────────────────────────
fallback_summary <- all_imputed %>%
  dplyr::group_by(method, X) %>%
  dplyr::summarise(
    n_total    = dplyr::n(),
    n_fallback = sum(fallback),
    rate_pct   = round(100 * mean(fallback), 3),
    .groups    = "drop"
  ) %>%
  dplyr::arrange(method, X)

message("\n── Fallback rates per method and masking fraction ──")
print(as.data.frame(fallback_summary), row.names = FALSE)

saveRDS(fallback_summary,
        file.path(out_dir, "fallback_summary.rds"))
fwrite(fallback_summary,
       file.path(out_dir, "fallback_summary.csv"))

# ────────────────────────────────────────────────────────────
# §7  [Optional] Diagnostic performance plots
#     (Bland–Altman and MAD/RMSE bar charts)
#     Safe to skip — outputs are informative but not required
#     by any downstream stage.
# ────────────────────────────────────────────────────────────

message("§7  Drawing diagnostic plots ...")

all_imputed <- all_imputed %>%
  mutate(
    X_pct = factor(
      paste0(round(X * 100, 0), "%"),
      levels = paste0(round(sort(unique(X)) * 100, 0), "%")
    )
  )

# ── §7.0  Qualifying periods + masked months figure ──────────────────────────

message("§7.0  Drawing qualifying periods and masked months figure ...")

# All qualifying continuous windows per site (for background shading)
qualifying_windows <- df_O18_trim %>%
  dplyr::group_by(Site) %>%
  dplyr::group_modify(~{
    dates     <- sort(.x$Date)
    all_months <- seq(floor_date(min(dates), "month"),
                      floor_date(max(dates), "month"), by = "1 month")
    present   <- all_months %in% floor_date(dates, "month")
    r         <- rle(present)
    ends      <- cumsum(r$lengths)
    starts    <- c(1L, ends[-length(ends)] + 1L)
    good      <- which(r$values & r$lengths >= min_continuous_months)
    if (length(good) == 0L) return(tibble())
    tibble(
      win_start = all_months[starts[good]],
      win_end   = all_months[ends[good]]
    )
  }) %>%
  dplyr::ungroup()

# Show masked months for a single representative replicate (X=8%, boot=1).
# This gives a sparse, readable view of what a typical masking looks like.
rep_X    <- 0.08
rep_boot <- 1L

# Full record extent per site (all paired observations)
full_extent <- dplyr::inner_join(
    df_O18_trim %>% dplyr::select(Site, Date),
    df_H2_trim  %>% dplyr::select(Site, Date),
    by = c("Site", "Date")
  ) %>%
  dplyr::group_by(Site) %>%
  dplyr::summarise(
    rec_start = floor_date(min(Date), "month"),
    rec_end   = floor_date(max(Date), "month"),
    .groups   = "drop"
  )

masked_example <- all_imputed %>%
  dplyr::filter(X == rep_X, boot == rep_boot) %>%
  dplyr::distinct(Site, Date) %>%
  dplyr::mutate(Month = floor_date(Date, "month"))

# Site order: by earliest qualifying window start, then name
site_order_fig <- qualifying_windows %>%
  dplyr::group_by(Site) %>%
  dplyr::summarise(first_win = min(win_start), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(first_win), Site) %>%
  dplyr::pull(Site)

full_extent$Site        <- factor(full_extent$Site,        levels = site_order_fig)
qualifying_windows$Site <- factor(qualifying_windows$Site, levels = site_order_fig)
masked_example$Site     <- factor(masked_example$Site,     levels = site_order_fig)

p_periods <- ggplot() +
  # Full paired record — thin grey line
  geom_segment(
    data = full_extent,
    aes(x = rec_start, xend = rec_end, y = Site, yend = Site),
    colour = "grey70", linewidth = 1.2, lineend = "butt"
  ) +
  # Qualifying continuous windows — blue filled band
  geom_segment(
    data = qualifying_windows,
    aes(x = win_start, xend = win_end, y = Site, yend = Site),
    colour = "#2166ac", linewidth = 5, alpha = 0.25, lineend = "butt"
  ) +
  geom_segment(
    data = qualifying_windows,
    aes(x = win_start, xend = win_end, y = Site, yend = Site),
    colour = "#2166ac", linewidth = 0.7, lineend = "butt"
  ) +
  # Masked months — red tick marks
  geom_tile(
    data = masked_example,
    aes(x = Month, y = Site),
    fill = "#b2182b", width = 25, height = 0.6
  ) +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y", expand = c(0.01, 0)) +
  scale_y_discrete(limits = site_order_fig) +
  labs(
    title    = paste0("Focus sites: data coverage, qualifying periods, and example masking (X=",
                      round(rep_X * 100), "%, replicate ", rep_boot, ")"),
    subtitle = "Grey line = full paired δ¹⁸O & δ²H record  |  Blue band = qualifying continuous run ≥84 months  |  Red = withheld months",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(size = 7),
    plot.title         = element_text(face = "bold", size = 10),
    plot.subtitle      = element_text(size = 8, color = "grey40")
  )

print(p_periods)
ggsave(
  file.path(out_dir, "qualifying_periods_and_masked_months.png"),
  p_periods, width = 14, height = 16, dpi = 150
)

# ─────────────────────────────────────────────────────────────────────────────

summarise_var <- function(df, var) {
  o <- paste0(var, "_orig"); i <- paste0(var, "_imp")
  df %>%
    select(X, X_pct, method, orig = all_of(o), imp = all_of(i)) %>%
    filter(is.finite(orig), is.finite(imp)) %>%
    mutate(err = orig - imp) %>%
    group_by(X, X_pct, method) %>%
    summarise(
      MAD  = mean(abs(err), na.rm = TRUE),
      RMSE = sqrt(mean(err^2, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(var = var, .before = 1)
}

perf_summary <- bind_rows(
  summarise_var(all_imputed, "O18"),
  summarise_var(all_imputed, "H2"),
  summarise_var(all_imputed, "d_ex")
)

make_ba_long <- function(df, var_label, orig_col, imp_col) {
  df %>%
    transmute(
      var   = var_label,
      method,
      X_pct,
      orig  = .data[[orig_col]],
      diff  = .data[[imp_col]] - .data[[orig_col]]
    ) %>%
    filter(is.finite(orig), is.finite(diff))
}

ba_long <- bind_rows(
  make_ba_long(all_imputed, "δ18O",    "O18_orig",  "O18_imp"),
  make_ba_long(all_imputed, "δ2H",     "H2_orig",   "H2_imp"),
  make_ba_long(all_imputed, "d-excess", "d_ex_orig", "d_ex_imp")
)

ba_bias <- ba_long %>%
  dplyr::group_by(
    .data[["var"]],
    .data[["X_pct"]],
    .data[["method"]]
  ) %>%
  dplyr::summarise(
    bias = mean(.data[["diff"]], na.rm = TRUE),
    .groups = "drop"
  )

make_ba_plot <- function(which_var, y_min = NULL, y_max = NULL) {
  dfv <- ba_long %>%
    dplyr::filter(.data[["var"]] == which_var)

  statv <- ba_bias %>%
    dplyr::filter(.data[["var"]] == which_var)

  p <- ggplot(dfv, aes(x = orig, y = diff)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
    geom_bin2d(bins = 60, alpha = 0.9) +
    scale_fill_viridis_c(name = "count", option = "C") +
    geom_hline(
      data = statv, aes(yintercept = bias),
      color = "steelblue", linewidth = 0.6
    ) +
    facet_grid(rows = vars(X_pct), cols = vars(method), scales = "free_x") +
    labs(
      title = paste0("Bland–Altman — ", which_var, " (imp − orig)"),
      x     = paste0(which_var, " observed"),
      y     = paste0(which_var, " (imp − orig)")
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title         = element_text(hjust = 0.5),
      strip.background   = element_rect(fill = "#F0F0F0", color = NA),
      panel.grid.minor   = element_blank()
    )

  if (!is.null(y_min) || !is.null(y_max))
    p <- p + coord_cartesian(ylim = c(y_min, y_max))

  p
}

# ── Y-axis limits per variable ───────────────────────────────────────────────
ba_ylim <- list(
  "δ18O"     = c(y_min = -25,  y_max =  25),
  "δ2H"      = c(y_min = -200, y_max = 200),
  "d-excess" = c(y_min = -20,  y_max =  20)
)
# ─────────────────────────────────────────────────────────────────────────────

p_ba_O18 <- make_ba_plot("δ18O",    ba_ylim[["δ18O"]][["y_min"]],    ba_ylim[["δ18O"]][["y_max"]])
p_ba_H2  <- make_ba_plot("δ2H",     ba_ylim[["δ2H"]][["y_min"]],     ba_ylim[["δ2H"]][["y_max"]])
p_ba_dex <- make_ba_plot("d-excess", ba_ylim[["d-excess"]][["y_min"]], ba_ylim[["d-excess"]][["y_max"]])

print(p_ba_O18); print(p_ba_H2); print(p_ba_dex)

ggsave(file.path(out_dir, "BA_d18O.png"),   p_ba_O18, width = 14, height = 10, dpi = 150)
ggsave(file.path(out_dir, "BA_d2H.png"),    p_ba_H2,  width = 14, height = 10, dpi = 150)
ggsave(file.path(out_dir, "BA_dexcess.png"),p_ba_dex, width = 14, height = 10, dpi = 150)

# ────────────────────────────────────────────────────────────
# §8  Build and save the minimal reproducibility bundle
# ────────────────────────────────────────────────────────────
# ────────────────────────────────────────────────────────────
# §8  Build and save the minimal reproducibility bundle
# ────────────────────────────────────────────────────────────

message("§8  Building bundle ...")

core_keep <- c("Linear", "Sinusoidal")

removed_dates_out <- all_imputed %>%
  dplyr::filter(.data$method %in% core_keep) %>%
  dplyr::distinct(.data$Site, .data$X, .data$boot, .data$Date)

core_combos_out <- removed_dates_out %>%
  dplyr::distinct(.data$Site, .data$X, .data$boot)

baseline_errors_out <- dplyr::bind_rows(
  
  # d18O baselines
  all_imputed %>%
    dplyr::filter(
      .data$method %in% core_keep,
      is.finite(.data$O18_orig),
      is.finite(.data$O18_imp)
    ) %>%
    dplyr::group_by(.data$Site, .data$X, .data$boot, .data$method) %>%
    dplyr::summarise(
      n_baseline = dplyr::n(),
      MAD  = mean(abs(.data$O18_orig - .data$O18_imp), na.rm = TRUE),
      RMSE = sqrt(mean((.data$O18_orig - .data$O18_imp)^2, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(isotope = "d18O"),
  
  # d2H baselines
  all_imputed %>%
    dplyr::filter(
      .data$method %in% core_keep,
      is.finite(.data$H2_orig),
      is.finite(.data$H2_imp)
    ) %>%
    dplyr::group_by(.data$Site, .data$X, .data$boot, .data$method) %>%
    dplyr::summarise(
      n_baseline = dplyr::n(),
      MAD  = mean(abs(.data$H2_orig - .data$H2_imp), na.rm = TRUE),
      RMSE = sqrt(mean((.data$H2_orig - .data$H2_imp)^2, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(isotope = "d2H")
  
) %>%
  dplyr::select(
    .data$Site,
    .data$X,
    .data$boot,
    .data$isotope,
    .data$method,
    .data$n_baseline,
    .data$MAD,
    .data$RMSE
  )

bundle <- list(
  removed_dates   = removed_dates_out,
  core_combos     = core_combos_out,
  baseline_errors = baseline_errors_out,
  df_O18_min      = df_O18_trim %>%
    dplyr::select(.data$Site, .data$Date, .data$O18, .data$Altitude),
  df_H2_min       = df_H2_trim %>%
    dplyr::select(.data$Site, .data$Date, .data$H2, .data$Altitude),
  stations_min    = stations %>%
    dplyr::filter(.data$Site %in% long_sites) %>%
    dplyr::select(.data$Site, .data$Altitude, .data$Longitude, .data$Latitude),
  dist_km         = dist_km
)

out_path <- file.path(out_dir, "minimal_SPbI_input_bundle.rds")

saveRDS(
  bundle,
  out_path,
  compress = "gzip"
)

message("\nBundle saved: ", normalizePath(out_path))
message("Component sizes:")

for (nm in names(bundle)) {
  message(
    "  ",
    formatC(nm, width = 18, flag = "-"),
    format(object.size(bundle[[nm]]), units = "MB")
  )
}

message("\nDone. Run STAGE I next.")

