# ============================================================
# Build example_input.xlsx — a synthetic input file that mirrors
# the structure the pipeline expects.
# Hatvani I.G. & Kern Z. — Mind the gap (HESS, 2026)
# https://github.com/istvan60/SPbI
#
# The real input (inputation_test6.xlsx) is compiled from the national
# monitoring networks cited in the manuscript and is not redistributed
# here. This script documents the required layout and produces a file
# the pipeline will run on end to end, so users can verify their own
# data is shaped correctly before substituting it.
#
# THE VALUES ARE FABRICATED. Do not use them for any scientific purpose.
# ============================================================

suppressPackageStartupMessages({
  library(writexl)
})

set.seed(42)

# Written to the working directory.
out_file <- "example_input.xlsx"

# ── Sites ────────────────────────────────────────────────────
# Altitude in m a.s.l.; coordinates in decimal degrees (WGS84).
# Group is carried through the real file but is NOT read by any script.
# Sites are placed within ~80 km of one another. SPbI averages neighbours
# inside spatial_radius_bootstrap_km (100 km by default), so a sparser
# network would make it fall back on most months and the example would not
# demonstrate the method. Altitudes span a realistic relief gradient so the
# lapse-rate correction has something to act on.
coords <- data.frame(
  Site      = c("Alpha", "Bravo", "Charlie", "Delta", "Echo", "Foxtrot"),
  Latitude  = c(47.40, 47.05, 46.70, 46.95, 47.55, 46.45),
  Longitude = c(15.20, 14.40, 14.90, 15.60, 14.10, 15.35),
  Altitude  = c(  120,  1220,   550,   890,  1980,   300),
  Group     = c("NET_A", "NET_A", "NET_B", "NET_A", "NET_A", "NET_B"),
  stringsAsFactors = FALSE
)

# ── Monthly series ───────────────────────────────────────────
# Dates are the 15th of each month. The pipeline floors to month, so the
# day component only needs to be consistent.
months_all <- seq(as.Date("1973-01-15"), as.Date("2024-12-15"), by = "1 month")

# Per-site record window, plus an explicitly protected continuous block.
#
# This matters: the pipeline keeps only sites with an uninterrupted run of
# >= min_continuous_months (84) of PAIRED d18O + d2H. Scattering gaps at
# random over a long record fragments every run below that threshold, so
# each qualifying site is given a protected stretch that is left gap-free,
# with gaps applied only outside it. Echo and Foxtrot have no protected
# block and are expected to drop out -- that is deliberate, so the example
# exercises the selection step rather than trivially passing it.
spec <- list(
  Alpha   = list(start = "1973-01-15", end = "2024-12-15", gap = 0.10,
                 keep_from = "1980-01-15", keep_to = "1999-12-15"),   # 240 mo
  Bravo   = list(start = "1980-01-15", end = "2024-12-15", gap = 0.12,
                 keep_from = "1990-01-15", keep_to = "2004-12-15"),   # 180 mo
  Charlie = list(start = "1976-01-15", end = "2019-12-15", gap = 0.10,
                 keep_from = "1985-01-15", keep_to = "1993-12-15"),   # 108 mo
  Delta   = list(start = "1995-01-15", end = "2024-12-15", gap = 0.12,
                 keep_from = "2005-01-15", keep_to = "2012-12-15"),   #  96 mo
  Echo    = list(start = "2001-01-15", end = "2024-12-15", gap = 0.15,
                 keep_from = NA,          keep_to = NA),              # drops out
  Foxtrot = list(start = "2015-01-15", end = "2024-12-15", gap = 0.15,
                 keep_from = NA,          keep_to = NA)               # drops out
)

# Synthetic isotope signal: seasonal cycle + altitude effect + noise,
# with d2H tied to d18O through a meteoric-water-line relationship so
# d-excess behaves plausibly.
make_site <- function(site) {
  s   <- spec[[site]]
  alt <- coords$Altitude[coords$Site == site]

  d <- months_all[months_all >= as.Date(s$start) & months_all <= as.Date(s$end)]

  # drop a random subset of months to create realistic gappiness, but never
  # inside the protected block (see `spec` above)
  keep <- runif(length(d)) > s$gap
  if (!is.na(s$keep_from)) {
    protected <- d >= as.Date(s$keep_from) & d <= as.Date(s$keep_to)
    keep <- keep | protected
  }
  d <- d[keep]

  mon   <- as.integer(format(d, "%m"))
  seas  <- -4.5 * cos(2 * pi * (mon - 1) / 12)          # summer-enriched
  base  <- -9.0 - 0.20 * (alt / 100)                    # isotopic lapse rate
  o18   <- round(base + seas + rnorm(length(d), 0, 1.4), 2)
  h2    <- round(8.0 * o18 + 10 + rnorm(length(d), 0, 4.0), 1)

  list(
    o18 = data.frame(Site = site, Date = d, O18 = o18,
                     Dexcess = round(h2 - 8 * o18, 1),
                     Group = coords$Group[coords$Site == site],
                     stringsAsFactors = FALSE),
    h2  = data.frame(Site = site, Date = d, H2 = h2,
                     Dexcess = round(h2 - 8 * o18, 1),
                     Group = coords$Group[coords$Site == site],
                     stringsAsFactors = FALSE)
  )
}

parts  <- lapply(coords$Site, make_site)
sheet_O18 <- do.call(rbind, lapply(parts, `[[`, "o18"))
sheet_H2  <- do.call(rbind, lapply(parts, `[[`, "h2"))

write_xlsx(
  list(O18 = sheet_O18, H2 = sheet_H2, coords = coords),
  path = out_file
)

cat("Wrote", out_file, "\n")
cat("  O18   :", nrow(sheet_O18), "rows\n")
cat("  H2    :", nrow(sheet_H2),  "rows\n")
cat("  coords:", nrow(coords),    "sites\n")
